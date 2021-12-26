#include "mpi.h"
#include <mpi-ext.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define  Max(a, b) ((a)>(b)?(a):(b))

#define  N   (1*64+2)
#define COORDINATOR_NUM 0

MPI_Comm mpi_comm_world_custom;  // represents a logical group of MPI processes
int error_occurred = 0;  // bool flag

float maxeps = 0.1e-7;
int itmax = 100;
float w = 0.5;
float eps, local_eps;

float A[N][N];
float tmp_A_row[N];

int rank, num_workers, rc;
int first_row, last_row, n_rows;

void relax();
void init();
void verify();
int master_job();


static void verbose_errhandler (MPI_Comm *comm, int *err, ...) {
//    char errstr[MPI_MAX_ERROR_STRING];
//    int size, num_failed, len;
//    MPI_Group group_failed;
//    int old_rank = rank;
//    error_occurred = 1;
//    MPI_Comm_size(mpi_comm_world_custom, &size);
//    MPIX_Comm_failure_ack(mpi_comm_world_custom);  // Acknowledge the current group of failed processes
//    MPIX_Comm_failure_get_acked(mpi_comm_world_custom, &group_failed);  // Get the group of acknowledged failures.
//    MPI_Group_size(group_failed, &num_failed);
//    MPI_Error_string(*err, errstr, &len);
//
//    MPIX_Comm_shrink(mpi_comm_world_custom, &mpi_comm_world_custom);
//    MPI_Comm_rank(mpi_comm_world_custom, &rank);
//    MPI_Comm_size(mpi_comm_world_custom, &size);
//    printf("New rank for process %d: %d\n", old_rank, rank);
//    MPI_Barrier(mpi_comm_world_custom);
//    int * ranks = malloc(sizeof(int)*size);
//    MPI_Gather(&old_rank, 1, MPI_INT, ranks, 1, MPI_INT, 0, mpi_comm_world_custom);
//    if (rank == COORDINATOR_NUM) {
//        int killed_proc_num;
//        for (int i = 0; i < size - 1; ++i) {
//            if (ranks[i + 1] - ranks[i] > 1) {
//                killed_proc_num = ranks[i] + 1;
//            }
//        }
//        printf("Killed proc: %d\n", killed_proc_num);
////        sprintf(killed_filename, "./logs/%d.txt", killed_proc_num);
//        int *read_array = calloc(splitted_size, sizeof(int));
////        read_from_file(killed_filename, read_array, splitted_size);
////        int *tmp_splitted_array = calloc(splitted_size, sizeof(int));
////        merge_serial(read_array, tmp_splitted_array, splitted_size);
////        save_into_file(new_killed_filename, read_array, splitted_size);
//    }
}


void initialize_glob_row_borders(int num_workers, int rank){
    // -2 because first and last row is zero
    n_rows = (N-2) / num_workers;
    first_row = n_rows * rank + 1;
    if (rank != num_workers-1){
        last_row = first_row + n_rows;
    } else {
        last_row = N-1;
    }
}

int main(int an, char **as) {

    // создаем группу процессов и область связи
    if ((rc = MPI_Init(&an, &as))) {
        printf("Ошибка запуска %d, выполнение остановлено\n", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
        return rc;
    }
    // получаем номер текущего процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // получаем общее кол-во процессов
    MPI_Comm_size(MPI_COMM_WORLD, &num_workers);
    mpi_comm_world_custom = MPI_COMM_WORLD;

    MPI_Errhandler errh;
    MPI_Comm_create_errhandler(verbose_errhandler, &errh);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, errh);
    MPI_Barrier(mpi_comm_world_custom);

    initialize_glob_row_borders(num_workers, rank);
    printf("rank %d: first_row %d, last_row %d\n", rank, first_row, last_row);

    struct timeval start, stop;
    double secs;
    if (!rank){
        gettimeofday(&start, NULL);
    }

    init();
    for (int it = 1; it <= itmax; it++) {
        eps = 0.;
        local_eps = 0.;
        relax();
        // if (!rank) {printf("it=%4i   eps=%f\n", it, eps);}
        if (eps < maxeps) break;
    }
    verify();

    if (!rank){
        gettimeofday(&stop, NULL);
        secs = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
        printf("time taken for thread=%d, N=%d: %f seconds\n", num_workers, N, secs);
    }

    /* Shut down MPI */
    MPI_Finalize();

    return 0;
}

/* matrix like, [NxN]
0 0 0 0 0
0 3 4 5 0
0 4 5 6 0
0 5 6 7 0
0 0 0 0 0 */
void init() {
    for (int i = 0; i <= N - 1; i++)
        for (int j = 0; j <= N - 1; j++){
            if (i == 0 || i == N - 1 || j == 0 || j == N - 1) A[i][j] = 0.;
            else A[i][j] = (1. + i + j);
        }
}


void relax() {
    MPI_Status status;
    MPI_Barrier(MPI_COMM_WORLD);

    int up_send_tag = 2, down_send_tag = 3;

    // меняются только нечётные
    for (int i = first_row; i < last_row; i++)
        for (int j = 1; j <= N - 2; j++){
            if ((i + j) % 2 == 1) {
                float b;
                b = w * ((A[i - 1][j] + A[i + 1][j] + A[i][j - 1] + A[i][j + 1]) / 4. - A[i][j]);
                local_eps = Max(fabs(b), local_eps);
                A[i][j] = A[i][j] + b;
            }}

    // int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    // IN buf	-	адрес начала расположения пересылаемых данных;
    // IN count	-	число пересылаемых элементов;
    // IN datatype	-	тип посылаемых элементов;
    // IN dest	-	номер процесса-получателя в группе, связанной с коммуникатором comm;
    // IN tag	-	идентификатор сообщения (аналог типа сообщения функций nread и nwrite PSE nCUBE2);
    // IN comm	-	коммуникатор области связи.
    if (rank!=0) {MPI_Send(A[first_row], N, MPI_FLOAT, rank-1, up_send_tag, MPI_COMM_WORLD);}
    if (rank!=num_workers-1){
        MPI_Recv(tmp_A_row, N, MPI_FLOAT, rank+1, up_send_tag, MPI_COMM_WORLD, &status);
        // (j + last_row) % 2 == 1   =>   j = 1+(last_row % 2)
        for (int j=1+(last_row % 2); j<=N-2; j+=2) {A[last_row][j] = tmp_A_row[j];}
    }
    if (rank!=num_workers-1) {MPI_Send(A[last_row-1], N, MPI_FLOAT, rank+1, down_send_tag, MPI_COMM_WORLD);}
    if (rank!=0){
        MPI_Recv(tmp_A_row, N, MPI_FLOAT, rank-1, down_send_tag, MPI_COMM_WORLD, &status);
        for (int j=1+((first_row-1) % 2); j<=N-2; j+=2) {A[first_row-1][j] = tmp_A_row[j];}
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // sleep(rank);
    // printf("\nrank %d, не чётный\n", rank);
    // for (int i = 0; i < N; ++i){
    //     for (int j = 0; j < N; ++j){
    //         printf("%f ", A[i][j]);
    //     }
    //     printf("\n");
    // }

    // меняются только чётные
    for (int i = first_row; i < last_row; i++)
        for (int j = 1; j <= N - 2; j++)
            if ((i + j) % 2 == 0) {
                float b;
                b = w * ((A[i - 1][j] + A[i + 1][j] + A[i][j - 1] + A[i][j + 1]) / 4. - A[i][j]);
                A[i][j] = A[i][j] + b;
            }

    if (rank!=0) {MPI_Send(A[first_row], N, MPI_FLOAT, rank-1, up_send_tag, MPI_COMM_WORLD);}
    if (rank!=num_workers-1){
        MPI_Recv(tmp_A_row, N, MPI_FLOAT, rank+1, up_send_tag, MPI_COMM_WORLD, &status);
        // (j + last_row) % 2 == 0   =>   j = (last_row % 2)
        for (int j=(last_row % 2); j<=N-2; j+=2) {A[last_row][j] = tmp_A_row[j];}
    }
    if (rank!=num_workers-1) {MPI_Send(A[last_row-1], N, MPI_FLOAT, rank+1, down_send_tag, MPI_COMM_WORLD);}
    if (rank!=0){
        MPI_Recv(tmp_A_row, N, MPI_FLOAT, rank-1, down_send_tag, MPI_COMM_WORLD, &status);
        for (int j=((first_row-1) % 2); j<=N-2; j+=2) {A[first_row-1][j] = tmp_A_row[j];}
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Reduce local_eps in each proc
    MPI_Allreduce(&local_eps, &eps, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

    // sleep(rank);
    // printf("\nrank %d, чётный\n", rank);
    // for (int i = 0; i < N; ++i){
    //     for (int j = 0; j < N; ++j){
    //         printf("%f ", A[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // MPI_Barrier(MPI_COMM_WORLD);
}


void verify() {
    float local_sum, sum;

    local_sum = 0.;
    int begin = first_row;
    int end = last_row;
    if (first_row == 1) {begin = 0;}
    if (last_row == N - 1) {end = N - 1;}

    for (int i = begin; i <= end; i++)
        for (int j = 0; j <= N - 1; j++){
            local_sum = local_sum + A[i][j] * (i + 1) * (j + 1) / (N * N);
        }

    MPI_Reduce(&local_sum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (!rank) {printf("  S = %f\n", sum);}
}

// mpicc -std=c99 -o run_1 redb_2d.c -lm
