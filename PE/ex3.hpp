#include <iostream>
#include <cmath>
#include <pthread.h>
using namespace std;

typedef struct{
    int n;
    int p;
    int cur_p;
    double *A_copy;
    double *b_copy;
    double *x;
    int *variable_order;
    int *thread_create_flag;
    int *control;

    void fill(int n, int p, int cur_p, double *A_copy, double *b_copy, double *x, int *variable_order, int *thread_create_flag, int *control)
    {
        this->n = n;
        this->p = p;
        this->cur_p = cur_p;
        this->A_copy = A_copy;
        this->b_copy = b_copy;
        this->x = x;
        this->variable_order = variable_order;
        this->thread_create_flag = thread_create_flag;
        this->control = control;
    }
} Args;

void synchronize(int total_threads);

void swap_rows(double* A, double *b, int row1, int row2, int n);

void swap_cols(double* A, int col1, int col2, int n);

template<typename T>
void swap_array_elements(int n, int el1, int el2, T *arr);

bool is_singular(double* A, int n);

int determine_row(int n, int pos);
 
int determine_col(int n, int pos);
 
double determine_max(double *A, int n, int k);

int determine_max_pos(double *A, int n, int k);

double find_max_in_row(double *A, int n, int row);

// void normalize_matrix(double *A, double *b, int n);

int normalize_matrix2(double *A, double *b, int n);

int solve_system(double *A, double *b, double *x, int n, int *variable_order);

int step(void *args, int i);

int reverse_step(void *args, int i);

void *solve_system_treaded(void *input_data);