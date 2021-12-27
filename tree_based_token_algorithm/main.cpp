#include <cstdio>
#include "mpi.h"
//#include <sys/stat.h>
#include <unistd.h>

#include <string>
#include <stdexcept>

const std::string file_name = "critical.txt";

void critical_section() {
    if (FILE * file = fopen(file_name.c_str(), "r")) {  // проверка наличия файла “critical.txt”
        throw std::invalid_argument("received negative value");  // TODO change exception
        //<сообщение об ошибке>;
    }
    FILE *file = fopen(file_name.c_str(), "w");//<создание файла “critical.txt”>;
    unsigned int microseconds = 10000;
    usleep(microseconds);//Sleep (<случайное время>);
    std::remove(file_name.c_str()); // уничтожение файла “critical.txt”
}

struct Node {
    int current_process;
    int parent_process;
    int left_process;
    int right_process;

    Node(int current_process,
         int parent_process,
         int left_process,
         int right_process) {
        this->current_process = current_process;
        this->parent_process = parent_process;
        this->left_process = left_process;
        this->right_process = right_process;
    }

    void call() {
        int number = 1;
        if (parent_process != -1) { // it is a root
            MPI_Recv(&number, 1, MPI_INT, parent_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        printf("process %d got the marker\n", current_process);
        critical_section();

        if (left_process != -1) {
            MPI_Send(&number, 1, MPI_INT, left_process, 0, MPI_COMM_WORLD);
            MPI_Recv(&number, 1, MPI_INT, left_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (right_process != -1) {
            MPI_Send(&number, 1, MPI_INT, right_process, 0, MPI_COMM_WORLD);
            MPI_Recv(&number, 1, MPI_INT, right_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (parent_process != -1) {
            MPI_Send(&number, 1, MPI_INT, parent_process, 0, MPI_COMM_WORLD);
        }
    }
};

int rank, num_workers, rc;

int main(int an, char **as) {
    // создаем группу процессов и область связи
    if ((rc = MPI_Init(&an, &as))) {
        printf("Ошибка запуска %d, выполнение остановлено\n", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
        return rc;
    }

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int parent_proc = world_rank / 2 - 1;  // -1 because we enumerate nodes starting from 0
    if (parent_proc < 0){
        parent_proc = -1;
    }
    int left_proc = world_rank * 2 - 1;
    if (left_proc >= world_size){
        left_proc = -1;
    }
    int right_proc = world_rank * 2 + 1 - 1;
    if (right_proc >= world_size){
        right_proc = -1;
    }
    Node node = Node(world_rank, parent_proc, left_proc, right_proc);
    node.call();


//    int number;
//    MPI_Request request;
//    MPI_Status status;
//    if (world_rank == 0) {
//        number = -1;
//        MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
//        critical_section();
//    } else if (world_rank == 1) {
//         MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        // void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status
////        MPI_Irecv(&number, 1, MPI_INT, MPI_ANY_TAG, 0, MPI_COMM_WORLD, &request);
//        // void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request
//        MPI_Wait(&request, &status);
//        critical_section();
//        printf("Process 1 received number %d from process 0\n",
//               number);
//    }
}
