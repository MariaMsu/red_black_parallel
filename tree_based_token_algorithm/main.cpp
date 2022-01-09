#include <cstdio>
#include "mpi.h"
#include <unistd.h>

#include <string>
#include <stdexcept>

const char* file_name = "critical.txt";

void critical_section() {
    if (FILE * file = fopen(file_name, "r")) {  // проверка наличия файла “critical.txt”
        throw std::runtime_error("two ore more process are in critical section simultaneously");
    }
    FILE *file = fopen(file_name, "w");  // создание файла “critical.txt”
    unsigned int microseconds = std::rand() % 1000000 ;
    usleep(microseconds);  // Sleep (<случайное время>);
    std::remove(file_name); // уничтожение файла “critical.txt”
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
        int number = 1;  // the very value is not important. We just send and receive anything
        if (parent_process != -1) {
            // if parent_process != -1 then it  is a root =>
            // the process already has the marker (according to the task, not to the tree-token algorithm )
            printf("process %d wait the marker from %d\n", current_process, parent_process);
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

    int parent_proc = (world_rank - 1) / 2;
    if (world_rank == 0) {  // it is the root
        parent_proc = -1;
    }
    int left_proc = world_rank * 2 + 1;
    if (left_proc >= world_size) {
        left_proc = -1;
    }
    int right_proc = world_rank * 2 + 2;
    if (right_proc >= world_size) {
        right_proc = -1;
    }
    Node node = Node(world_rank, parent_proc, left_proc, right_proc);
    node.call();

    /* Shut down MPI */
    MPI_Finalize();
    return 0;
}
