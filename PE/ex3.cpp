#include <iostream>
#include <cmath>
#include <thread>
#include <pthread.h>
#include <unistd.h>
#include "ex3.hpp"
#include "M3.hpp"
#include "Exceptions3.hpp"
using namespace std;

// напиши параллельную программу для метода гаусса с выбором главного элемента по всей матрице на с++. Используй библиотеку pthread.h, матрицу задай одномерным массивом, отлавливай вырожденные матрицы

void synchronize(int p)
{
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;

	pthread_mutex_lock(&mutex);

	threads_in++;
        if (threads_in >= p)
	{
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} 
    else
    {
                while (threads_in < p)
        {
			pthread_cond_wait(&condvar_in,&mutex);
        }
    }

	threads_out++;
        if (threads_out >= p)
	{
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} 
    else
    {
                while (threads_out < p)
        {
			pthread_cond_wait(&condvar_out,&mutex);
        }
    }
	pthread_mutex_unlock(&mutex);
}

void swap_rows(double* A, double *b, int row1, int row2, int n)
{
    if ((row1 < 0 || row1 >= n) || (row2 < 0 || row2 >= n))
    {
        // free(A);
        // free(b);
        //   throw MyException("\n\nEXCEPTION: swap_rows: index out of range\n\n\n");
    }
    for (int i = 0; i < n; i++) {
        double tmp = A[row1 * n + i];
        A[row1 * n + i] = A[row2 * n + i];
        A[row2 * n + i] = tmp;
    }
    double tmp = b[row1];
    b[row1] = b[row2];
    b[row2] = tmp;
}
void swap_cols(double* A, int col1, int col2, int n)
{
    if ((col1 < 0 || col1 >= n) || (col2 < 0 || col2 >= n))
    {
        //free(A);
        //   throw MyException("\n\nEXCEPTION: swap_cols: index out of range\n\n\n");
    }
    for (int i = 0; i < n; i++) {
        double tmp = A[i * n + col1];
        A[i * n + col1] = A[i * n + col2];
        A[i * n + col2] = tmp;
    }
}

template<typename T>
void swap_array_elements(int n, int el1, int el2, T *arr)
{
    if (el1 < 0 || el1 >= n || el2 < 0 || el2 >= n)
    {

        //free(arr);
        //   throw MyException("\n\nEXCEPTION: swap_array_elements: index out of range\n\n\n");
    }
    T tmp = arr[el1];
    arr[el1] = arr[el2];
    arr[el2] = tmp;
}

bool is_singular(double* A, int n) {
    for (int i = 0; i < n; i++) {
        if (fabs(A[i * n + i]) < 1e-16) {
            return true;
        }
    }
    return false;
}

int determine_row(int n, int pos)
{
    return pos/n;
}

int determine_col(int n, int pos)
{
    return pos % n;
}

double determine_max(double *A, int n, int k)
{
    double max_val = A[k*n + k];
    for (int r = k; r < n; r++)
    {
        for (int q = k; q < n; q++)
        {
            if (fabs(A[r*n + q]) > fabs(max_val))
            {
                max_val = A[r*n + q];
            }
        }
    }
    return max_val;
}

int determine_max_pos(double *A, int n, int k)
{
    int max_pos = k*n + k;
    double max_val = A[k*n + k];
    for (int i = k; i < n; i++)
    {
        for (int j = k; j < n; j++)
        {
            if (fabs(A[i*n + j]) > fabs(max_val))
            {
                max_val = A[i*n + j];
                max_pos = i*n + j;
            }
        }
    }
    return max_pos;
}

double find_max_in_row(double *A, int n, int row)
{
    double max_val = A[row*n];
    for (int i = 1; i < n; i++)
    {
        if (fabs(A[row*n + i]) > fabs(max_val))
        {
            max_val = A[row*n + i];
        }
    }
    return max_val;
}

// void normalize_matrix(double *A, double *b, int n) {
//     for (int i = 0; i < n; i++) {
//         double max_val = fabs(A[i * n + i]);

//         for (int j = i + 1; j < n; j++) {
//             if (fabs(A[j * n + i]) > max_val) {
//                 max_val = fabs(A[j * n + i]);
//             }
//         }

//         if (max_val < 1e-300) {
//             // free(A);
//             // free(b);
//             //   throw MyException("\n\nEXCEPTION: normalize_matrix: A is singular\n\n\n");
//         }

//         for (int j = 0; j < n; j++) {
//             A[i * n + j] /= max_val;
//         }

//         b[i] /= max_val;
//     }
// }

int normalize_matrix2(double *A, double *b, int n)
{
    for (int i = 0; i < n; i++)
    {
        double max_val = find_max_in_row(A, n, i);
        if (fabs(max_val) < 1e-300)
        {
            // free(A);
            // free(b);
            return -5;
            // throw MyException("\n\nEXCEPTION: normalize_matrix: A is singular\n\n\n");
        }
        for (int j = 0; j < n; j++)
        {
            A[i*n + j] /= max_val;
        }
        b[i] /= max_val;
    }
    return 0;
}

int solve_system(double *A_copy, double *b_copy, double *x, int n, int *variable_order)
{
    // p = 0;
    //double *b_copy = (double *)malloc(n*sizeof(double));
    //    for (int i = 0; i < n; i++)
    //    {
    //        b_copy[i] = b[i];
    //    }
    //    //double *A_copy = (double *)malloc(n*n*sizeof(double));
    //    for (int i = 0; i < n*n; i++)
    //    {
    //        A_copy[i] = A[i];
    //    }
    int c1 = normalize_matrix2(A_copy, b_copy, n); // zero row catched here
    if (c1 != 0)
    {
        return c1;
    }
    //cout << "Normalized A, b:" << endl;
    //print_expanded(m, n, n, A_copy, b_copy);
    //int max_positions[n-1];
    //int *variable_order = (int *)malloc(n * sizeof(int));


    // for (int i = 0; i < n; i++)
    // {
    //     variable_order[i] = n-1-i;
    // }


    int cols_swap_count = 0;
    //cout << "1" << endl;
    for (int i = 0; i < n-1; i++)
    {
        int max_pos = determine_max_pos(A_copy, n, i);
        int max_row = determine_row(n, max_pos);
        int max_col = determine_col(n, max_pos);
        double max_val = determine_max(A_copy, n, i);
        //cout << "max_val : " << max_val << endl;
        //if (max_val < 1e-16)
//        {
//            throw MyException("\n\nEXCEPTION: solve_system: Matrix is singular\n\n\n");
//        }
//        if (max_val < 1e-300)
//        {
//            free(variable_order);
//            free(b_copy);
//            free(A_copy);
//            throw MyException("\n\nEXCEPTION: solve_system: Matrix is singular\n\n\n");
//        }
        //cout << "found max : " << max_val << endl;
        //cout << "max pos: " << max_pos << endl;
        //cout << "max_col: " << max_col << " i : " << i << endl;
        //max_positions[i] = max_col;

        if (max_row != i)
        {
            //cout << "Swapped rows: " << max_row << " " << i << endl;
            swap_rows(A_copy, b_copy, max_row, i, n);
            max_row = i;
        }
        if (max_col != i)
        {
            //cout << "Swapped cols: " << max_col << " " << i << endl;
            swap_cols(A_copy, max_col, i, n);
            swap_array_elements(n, max_col, i, variable_order);
            max_col = i;
            cols_swap_count++;
        }
        //cout << "1/1" << endl;
        for (int j = i+1; j < n; j++)
        {
//            {
//                free(variable_order);
//                free(b_copy);
//                free(A_copy);
//                throw MyException("\n\nEXCEPTION: solve_system: Matrix is singular\n\n\n");
//            }
            //cout << "max_val" << endl;
            double coef = -A_copy[n*j + i]/ max_val;
            //cout << "coef" << endl;
            for (int k = i; k < n; k++)
            {
                A_copy[j*n + k] += A_copy[i*n + k] * coef;
            }
            b_copy[j] += coef * b_copy[i];
            //cout << "1/1/" << j << endl;
        }

//        if (i > 0)
//        {
//            double avg = 0;
//            for (int k = 0; k < i; k++)
//            {
//                avg += A_copy[k*n + k];
//            }
//            avg /= i;
//            if (A_copy[i*n + i]/avg < 1e-15)
//            {
//                throw MyException("\n\nEXCEPTION: solve_system: matrix is singular\n\n\n");
//            }
//         //print_expanded(m, n, n, A_copy, b_copy);
//        }
//        if (abs(A_copy[i*n + i]) < 1e-15)
//        {
//            throw MyException("\n\nEXCEPTION: solve_system: matrix is singular\n\n\n");
//        }

//        for (int i = 0; i < n; i++)
//        {
//            for (int j = 0; j < n; j++)
//            {
//                cout << A_copy[i*n + j] << " ";
//            }
//            cout << "| " << b_copy[i] << endl;
//        }
    }
    if (is_singular(A_copy, n))
    {
        // free(A);
        // free(b);
        // free(x);
        //free(variable_order);
        //free(b_copy);
        //free(A_copy);


        return -6;
        //   throw MyException("\n\nEXCEPTION: solve_system: Matrix is singular 2\n\n\n");
    }
    double avg = 0;
    for (int i = 0; i < n-1; i++)
    {
        avg += A_copy[i*n + i];
    }
    avg /= (n-1);
    if (avg/A_copy[(n-1)*n + n-1] > 1e15)
    {
        return -7;
        //   throw MyException("\n\nEXCEPTION: solve_system: matrix is singular 3\n\n\n");
    }
    //cout << "2" << endl;
    //cout << "made triangular" << endl;
    //print_expanded(m, n, n, A, b);
//    if (is_singular(A_copy, n))
//    {
//        // free(A);
//        // free(b);
//        // free(x);
//        //free(variable_order);
//        //free(b_copy);
//        //free(A_copy);
//        throw MyException("\n\nEXCEPTION: solve_system: Matrix is singular\n\n\n");
//    }
    //cout << "3" << endl;
    //cout << endl << "REVERSE" << endl << endl;
    for (int i = n-1; i > 0; i--)      // reverse(matrix is triangular)
    {
        for (int j = i-1; j > -1; j--)
        {
            double coef = - A_copy[j * n + i] / A_copy[i*n + i];
            A_copy[j*n + i] += coef * A_copy[i*n + i];
            b_copy[j] += coef * b_copy[i];
        }
    }

    //cout << "4" << endl;
    //print_expanded(m, n, n, A, b);
    //cout << "Making E:" << endl;
    for (int i = 0; i < n; i++)
    {
        //cout << "b[" << i << "]:" << b[i] << " -> ";
        b_copy[i] /= A_copy[i*n + i];
        //cout << b[i]  << " divided by " << A[i*n + i] << endl;
        A_copy[i*n + i] /= A_copy[i*n + i];
    }
    //cout << "done reverse way" << endl;
    //print_expanded(m, n, n, A, b);
    // for (int i =  cols_swap_count - 1; i >= 0; i--)
    // {
    //     //cout << "trying to swap: " << n-1 - i << endl;
    //     swap_cols(A, i, variable_order[i], n);
    //     //cout << "swapped: " << n-1 - i << endl;
    // }
    // print_expanded(m, n, n, A, b);
    // cout << "swapped back" << endl;
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         if (fabs(A[j*n + i] - 1e+0) < 1e-6)
    //         {
    //             x[j] = b[i];
    //         }
    //     }
    // }
    //print_expanded(m, n, n, A, b);
    //cout << "variable order:\n";
    for (int i = 0; i < n; i++)
    {
        int cur_var = variable_order[i];
        if ((A_copy[i*n + cur_var] - 1) < 1e-16)
            x[n-1-cur_var] = b_copy[i];
        // reversing variable  order back
    }
    return 0;
//    free(variable_order);
//    free(b_copy);
//    free(A_copy);
}








void *solve_system_treaded(void *input_data)
{
    //cout << "START" << endl;
    Args *args = (Args *)input_data;
    int n = args->n;
    int p = args->p;
    int cur_p = args->cur_p;

    int min = n * cur_p / p;        // bogachev style
    int max = n * (cur_p + 1) / p;

    while(args->thread_create_flag[0] != 1)   // first thread is created -- cuz we use this after all of the threads are created
                                                // so we can either have all of have none
    {
        usleep(1 + cur_p / 5000);
        if (args->thread_create_flag[0] == -1)
        {
            cout << "\n\nEXCEPTION: solve_system: did not create threads to use them\n\n\n"; 
            args->control[cur_p] = -12;     /////////
            pthread_exit(nullptr);
        }
    }
    synchronize(p);



    // ---------------------- normalize_matrix() -----------------------------


    double max_value_n;
    for (int i = min; i < max; i++)     // zero rows are caught here
    {
        max_value_n = abs(args->A_copy[i*n + i]);
        for (int j = 0; j < n; j++)     ///////// (int j = i+1; j < n)
        {
            if (abs(args->A_copy[i*n + j]) > max_value_n)     ///////// A_copy[j*n + i]
                max_value_n = abs(args->A_copy[i*n + j]);
        }
        if (max_value_n < 1e-300)
        {
            cout << "\n\nEXCEPTION: normalize_matrix: matrix is singular\n\n\n";
            args->control[cur_p] = -5;
            pthread_exit(nullptr);
        }

        for (int j = 0; j < n; j++)     /// int j = i или i+1; 
        {
            args->A_copy[i*n + j] /= max_value_n;
        }
        args->b_copy[i] /= max_value_n;
        // synchronize(p);
    }
    synchronize(p);     //////////
//    cout << "NORMALIZED" << endl;
//    print_expanded(5, n, n, args->A_copy, args->b_copy);



    // -----------------------------------------------------------------------



    if (n == 1 && args->A_copy[0] < 1e-16)      // check for 1-dimensional matrix
    {
        cout << "\n\nEXCEPTION: solve_system: n = 1, matrix is singular\n\n\n";
        args->control[cur_p] = -13;
        pthread_exit(nullptr);
    }
    else if (n == 1 && args->A_copy[0] > 1e-16)
    {
        args->x[0] = args->b_copy[0];
    }
    
//    cout << "\nNORMALIZED:\n\n";
//    print_expanded(10, n, n, args->A_copy, args->b_copy);


    // --------------- determine_max_pos() + make tridiagonal ----------------


    for (int k = 0; k < n-1; k++)
    {
        int max_row = k;
        int max_col = k;
        //int max_pos = k*n + k;
        double max_value = args->A_copy[k*n + k];
        for (int i = k; i < n; i++)
        {
            for (int j = k; j < n; j++)
            {
                if (args->A_copy[i*n + j] > abs(max_value))
                {
                    max_value = args->A_copy[i*n + j];
                    max_row = i;
                    max_col = j;
                }
            }
        }
        if (max_row != k)
        {
            swap_rows(args->A_copy, args->b_copy, max_row, k, n);

            max_row = k;
        }
        if (max_col != k)
        {
            swap_cols(args->A_copy, max_col, k, n);
            swap_array_elements(n, max_col, k, args->variable_order);
            max_col = k;
        }
        double coef;
        for (int i = min; i < max; i++)
        {
            if (i > k)
            {
                coef = -args->A_copy[i*n + k] / max_value;
                for (int j = k; j < n; j++)
                {
                    args->A_copy[i*n + j] += args->A_copy[k*n + j] * coef;
                }
                args->b_copy[i] += coef * args->b_copy[k];
            }
        }
        synchronize(p);
    }
    //synchronize(p);
    //cout << "MADE TRIANGULAR" << endl;
    //print_expanded(5, n, n, args->A_copy, args->b_copy);
    //return 0;

    // ----------------------------- bad matrices ------------------------------------


    if (is_singular(args->A_copy, n))
    {
        cout << "\n\nEXCEPTION: solve_system: matrix is singular 2\n\n\n";
        args->control[cur_p] = -6;
        pthread_exit(nullptr);
    }

    double avg = 0;
    for (int i = 0; i < n-1; i++)
    {
        avg += args->A_copy[i*n + i];
    }
    avg /= (n-1);
    if (abs(avg / args->A_copy[(n-1)*n + n - 1]) > 1e15)
    {
        cout << "\n\nEXCEPTION: solve_system: matrix is singular 3\n\n\n";
        args->control[cur_p] = -7;
        pthread_exit(nullptr);
    }


    // --------------------------- reverse --------------------------------


    //cout << "REVERSE:" << endl;         // maybe needs to be paralellized
    for (int i = n-1; i  > 0; i--)
    {
        for (int j = i-1; j > -1; j--)
        {
            double coef = -args->A_copy[j*n + i] / args->A_copy[i*n + i];
            args->A_copy[j*n + i] += coef * args->A_copy[i*n + i];
            args->b_copy[j] += coef * args->b_copy[i]; 
        }
    }

    for (int i = 0; i < n; i++)
    {
        args->b_copy[i] /= args->A_copy[i*n + i];
        args->A_copy[i*n + i] /= args->A_copy[i*n + i];
        //args->A_copy[i*n + i] = 1;
    }
    for (int i = 0; i < n; i++)
    {
        int cur_var = args->variable_order[i];
        if ((args->A_copy[i*n + i] - 1) < 1e-16)
            args->x[n-1-cur_var] = args->b_copy[i];
    }
    return 0;
}