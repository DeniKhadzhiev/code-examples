#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
//#include "Exceptions3.hpp"
using namespace std;

typedef struct {
    int n;
    int p;
    int cur_p;
    double *A;
    double *b;
    double *x;
    double *residual;
    int *thread_create_flag;

    void fill(int n, int p, int cur_p, double *A, double *b, double *x, double *residual, int *thread_create_flag)
    {
        this->n = n;
        this->p = p;
        this->cur_p = cur_p;
        this->A = A;
        this->b = b;
        this->x = x;
        this->residual = residual;
        this->thread_create_flag = thread_create_flag;
    }
} R_args;

double f(int n, int k, int i, int j);

void make_b(int n, double *A, double *x);

void make_solution(int n, double *x);

bool special_symbol(char a);

int read(int n, string &in, double *matr_arr);

void generate(int n, int k, double *A);

void print(int max_displayed, int n, int m, double *arr);

void print_expanded(int max_displayed, int n, int m, double *A, double *b);

double eucl_v_norm(double *x, int n);

double max_v_norm(double *x, int n);

//double eucl_m_norm(double *A, int rows, int cols)
//{
//    double sum = 0;
//    for(int i = 0; i < rows; i++)
//    {
//        for (int j = 0; j < cols; j++)
//        {
//            sum += A[i*cols + j];
//        }
//    }
//    return sqrt(sum);
//}

//double max_m_norm(double *A, int rows, int cols)
//{
//    double max = A[0];
//    for(int i = 0; i < rows; i++)
//    {
//        for (int j = 0; j < cols; j++)
//        {
//            if (A[i*rows + j] > max)
//                max = A[i*rows + j];
//        }
//    }
//    return max;
//}

double matrix_norm(int n, int m, double *A);

double residual_norm(int n,  double *A, double* b, double *x);

void *residual_norm_threaded(void *input_r_data);

void *calculateNorm(void* input_r_data);

double error_norm(int n, double *x, double *y);

void read_generate(int n, string in, int k, double *matr_arr);

//void check_data(int n, int m, int k);

bool check_symbols(char *a, int n);
