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
//<завершение работы программы>;
//        fclose(file);
//        return true;
    }
    FILE *file = fopen(file_name.c_str(), "w");//<создание файла “critical.txt”>;
    unsigned int microseconds = 10000;
    usleep(microseconds);//Sleep (<случайное время>);
    printf("SLEEP");
    std::remove(file_name.c_str()); // уничтожение файла “critical.txt”
}

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

    int number;
    if (world_rank == 0) {
        number = -1;
        MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        critical_section();
    } else if (world_rank == 1) {
        MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process 1 received number %d from process 0\n",
               number);
    }
}
