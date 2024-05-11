#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <pthread.h>
#include <unistd.h>
#include "ex3.hpp"
#include "M3.hpp"
#include "Exceptions3.hpp"
using namespace std;

double f(int n, int k, int i, int j)
{
    if (k == 1)
        return double(n - max(i+1, j+1) + 1);
    else if (k == 2)
        return double(max(i+1, j+1));
    else if (k == 3)
        return double(abs(i-j));
    else if (k == 4)
        return 1./(i+j+1);
    return 0;
}

void make_b(int n, double *A, double *x)
{
    for (int i = 0; i < n; i++)
    {
        double x_i = 0;
        for (int k = 0; k < n; k++)
        {
            if (k%2==0)
            {
                x_i += A[i*n + k];
            }
        }
        x[i] = x_i;
    }
}

void make_solution(int n, double *x)
{
    for (int i = 0; i < n; i++)
    {
        i % 2 == 0 ? x[i] = 1 : x[i] = 0;
    }
}

bool special_symbol(char a)
{
    string s("!@#$%^&*(),?/|:;[]{}~`=_");
    for (int i = 0; i < int(s.length()); i++)
    {
        if(a == s[i])
            return true;
    }
    return false;
}

int read(int n, string &in, double *matr_arr)
{
    if (n <= 0)
    {
        //free(matr_arr);
        return -1;
        //throw MyException("\n\nEXCEPTION: read: Not able to read matrix\n\n\n");
    }
    ifstream input;
    input.open(in);
    if (!input.is_open())
    {
        //free(matr_arr);
        return -2;
        //throw MyException("\n\nEXCEPTION: read: Unable to open file for reading matrix A\n\n\n");
    }
    int i = 0;
    char a[30];
    while(input >> a && i< n*n)
    {
        // if (fabs(atof(a) - 0) < 1e-200 && a[0] != '0')
        // {
        //     input.close();
        //     free(matr_arr);
        //     throw MyException("EXCEPTION: read: wrong format of file, found letters\n");
        // }
        for (int j = 0; j < 30; j++)
        {
            if ((isalpha(a[j]) && a[j] != 'e')|| special_symbol(a[j]))
            {
                input.close();
                //free(matr_arr);
                return -3;
                //throw MyException("\n\nEXCEPTION: read: wrong format of file, found letters or special symbols\n\n\n");
            }
        }
        matr_arr[i] = atof(a);
        i++;
    }
    if (i < n*n)
    {
        //free(matr_arr);
        return -4;
        //throw MyException("\n\nEXCEPTION: read: Not enough numbers\n\n\n");
    }
    return 0;
}

// void read(int n, string &in, double *matr_arr)
// {
//     if (n <= 0)
//     {
//         free(matr_arr);
//         throw MyException("EXCEPTION: read: Not able to read matrix\n");
//     }
//     ifstream input(in);
//     if (!input)
//     {
//         free(matr_arr);
//         throw MyException("EXCEPTION: read: Unable to open file for reading matrix A\n");
//     }
//     int i = 0;
//     char a[10];
//     while(input >> a && i< n*n)
//     {
//         if (fabs(atof(a) - 0) < 1e-200 && a[0] != '0')
//         {
//             input.close();
//             free(matr_arr);
//             throw MyException("EXCEPTION: read: wrong format of file\n");

//         }
//         matr_arr[i] = atof(a);
//         i++;
//     }
//     if (i < n*n)
//     {
//         free(matr_arr);
//         throw MyException("EXCEPTION: read: Not enough numbers\n");
//     }
// }



void generate(int n, int k, double *A)
{
    if (k == 0)
    {}
    else
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                double cur = f(n, k, i, j);
                A[i*n+ j] = cur;
            }
        }
    }
}

void print(int max_displayed, int n, int m, double *arr)
{
    int rows = min(max_displayed, n);
    int cols = min(max_displayed, m);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            cout << setprecision(3) << scientific << arr[i*m + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print_expanded(int max_displayed, int n, int m, double *A, double *b)
{
    int rows = min(max_displayed, n);
    int cols = min(max_displayed, m);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            cout << setprecision(3) << scientific << A[i*m + j] << " ";
        }
        cout << " | " << b[i] << endl;
    }
    cout << endl;
}

double eucl_v_norm(double *x, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += x[i]*x[i];
    }
    return sqrt(sum);
}

double max_v_norm(double *x, int n)
{
    double max = x[0];
    for (int i = 1; i < n; i++)
    {
        if (x[i] > max)
            max = x[i];
    }
    return max;
}

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

double matrix_norm(int n, int m, double *A) // n - rows, m - cols
{
    double max_sum = 0;
    for (int i = 0; i < n; i++)
    {
        double cur_sum = 0;
        for (int j = 0; j < m; j++)
        {
            cur_sum += abs(A[i*m + j]);
        }
        if (cur_sum > max_sum)
            max_sum = cur_sum; 
    }
    return max_sum;
}

double residual_norm(int n,  double *A, double* b, double *x)
{
    // cout << endl;
    // for (int i = 0; i < n; i++)
    // {
    //     cout << "x[" << i << "]: " << x[i] << endl;
    // }
    // for (int i = 0; i < n; i++)
    // {
    //     cout << "b[" << i << "]: " << b[i] << endl;
    // }

    double y;
    double *Ax_b = (double *)malloc(n*sizeof(double));
    for (int i = 0; i < n; i++)
    {
        Ax_b[i] = 0;
        for(int j = 0; j < n; j++)
        {
            Ax_b[i] += A[i*n + j] * x[j];
        }
        // cout << "Ax_b[" << i << "]: " << Ax_b[i] << endl;
        Ax_b[i] -= b[i];
        // cout << "Ax_b[" << i << "]: " << Ax_b[i] << endl;
    }
    //double Ax_b_norm = eucl_v_norm(Ax_b, n);
    //double b_norm = eucl_v_norm(b, n);
    double Ax_b_norm = matrix_norm(n, 1, Ax_b);
    double b_norm = matrix_norm(n, 1, b);
    //cout << "Ax - b norm: " << Ax_b_norm << endl;
    //cout << "b norm: " << b_norm << endl;
    y = Ax_b_norm/b_norm;
    free(Ax_b);
    return y;
}

void *residual_norm_threaded(void *input_r_data)      // previous written functions do not help :(
{
    cout << "START" << endl;
    R_args *r_args = (R_args *)input_r_data;
    int n = r_args->n;
    int p = r_args->p;
    int cur_p = r_args->cur_p;

    while(r_args->thread_create_flag[0] != 1)   // first thread is created -- cuz we use this after all of the threads are created
                                                // so we can either have all of have none
    {
        usleep(1 + cur_p / 5000);
        if (r_args->thread_create_flag[0] == -1)
        {
            cout << "\n\nEXCEPTION: did not create threads to use them\n\n\n"; 
            //exit();
        }
    }
    synchronize(p);

    //cout << "Kostyl" << endl;

    //int j = 0;
    double max_sum = 0, cur_sum;//, tmp;
    int min = n * cur_p / p;        // bogachev style
    int max = n * (cur_p + 1) / p;

    //synchronize(p);
    for (int i = min; i < max; i++)
    {
        //cout << cur_p << endl;
        cur_sum = 0;
        for (int j = 0; i < n; j++)
        {
            cur_sum += r_args->A[i*n + j] * r_args->x[j];
            //cout << cur_sum << endl;
        }
        cur_sum -= r_args->b[i];
        //cout << "cur_sum" << endl;
        if (cur_sum > max_sum)
            max_sum = cur_sum;
        //cout << "max_sum" << endl;
        synchronize(p);
    }

    r_args->residual[cur_p] = max_sum;

    synchronize(p);
    cout << "Max_sum and synchronize" << endl;

    double norm = r_args->residual[0];
    if (cur_p == 0)
    {
        for (int i = 1; i < p; i++)     // с нуля можно начать
        {
            if (r_args->residual[i] > norm)
            {
                norm = r_args->residual[i];
            }
        }
    }

    r_args->residual[0] = norm;
    cout << "END" << endl;
    return 0;
}

void* calculateNorm(void *input_r_data) {
    R_args *r_args = (R_args *)input_r_data;
    int n = r_args->n;
    int p = r_args->p;
    int cur_p = r_args->cur_p;

    while(r_args->thread_create_flag[0] != 1)   // first thread is created -- cuz we use this after all of the threads are created
                                                // so we can either have all of have none
    {
        usleep(1 + cur_p / 5000);
        if (r_args->thread_create_flag[0] == -1)
        {
            cout << "\n\nEXCEPTION: did not create threads to use them\n\n\n"; 
            //exit();
        }
    }
    synchronize(p);

    int min = n * cur_p / p;
    int max = n * (cur_p+1) / p;

    // Вычисление нормы невязки для каждого потока
    for (int i = min; i < max; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++) {
            sum += r_args->A[i*n + j] * r_args->x[j];
        }
        r_args->residual[cur_p] = std::max(r_args->residual[cur_p], std::abs(sum - r_args->b[i]));
        //synchronize(p);
    }
    //synchronize(p);
    pthread_exit(nullptr);
}

double error_norm(int n, double *x, double *y)
{
    double *er = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
    {
        er[i] = x[i] - y[i];
    }
    double eu_norm = matrix_norm(n, 1, er);
    free(er);
    return eu_norm;
}

void read_generate(int n, string in, int k, double *matr_arr)
{
    if (k == 0)
    {
        read(n, in, matr_arr);
    }
    else
    {
        generate(n, k, matr_arr);
    }
}

// void check_data(int n, int m, int k)
// {
//     if (n<=0 || m <= 0 || k < 0 || k > 4)
//         throw MyException("\n\nEXCEPTION: check_data: n, m or k are incorrect\n\n\n");
// }