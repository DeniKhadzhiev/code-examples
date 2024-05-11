#include <iostream>
using namespace std;

void free_all(double *A, double *b, double *A_copy, double *b_copy, double *solution, double *x, int *variable_order)
{
    free(A);
    free(b);
    free(A_copy);
    free(b_copy);
    free(variable_order);
    free(x);
    free(solution);
}

// добавить структуры и их тоже очищать