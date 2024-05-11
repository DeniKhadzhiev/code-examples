#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <cstring>

#include <pthread.h>
#include <unistd.h>

#include "M3.hpp"
#include "Exceptions3.hpp"
#include "ex3.hpp"
#include "get_time.hpp"

using namespace std;

void free_all(double *A, double *b, double *A_copy, double *b_copy, double *solution, double *x, 
                int *variable_order, pthread_t *threads, double *residual, int *control, Args *args, R_args *r_args, double *time)
{
    if (A) free(A);
    if (b) free(b);
    if (A_copy) free(A_copy);
    if (b_copy) free(b_copy);
    if (variable_order) free(variable_order);
    if (x) free(x);
    if (solution) free(solution);
    if (threads) free(threads);
    if (residual) free(residual);
    if (control) free(control);
    if (args) free(args);
    if (r_args) free(r_args);
    if (time) free(time);
}

int main(int argc, char* argv[])
{
    if (argc != 6 && argc != 5)
    {
        std::cout << "\n\nEXCEPTION: error in command line arguments -- insufficient. There are only " << argc << "\n\n\n";
        return -1;
    }
    for (int i = 1; i < 5; i++)
    {
        for (int j = 0; j < static_cast<int>(strlen(argv[i])); j++)
        {
            if (isalpha((char)argv[i][j]))
            {
                std::cout << "\n\nEXCEPTION: error in command line arguments -- there are characters\n\n\n";
                return -1;
            }
        }
    }
    int n = atoi(argv[1]);
    int p = atoi(argv[2]);
    int m = atoi(argv[3]);
    int k = atoi(argv[4]);

    if ((argc == 5 && k == 0) || (argc == 6 && k != 0))
    {
        std::cout << "\n\nEXCEPTION: error in command line arguments \n\n\n";
        return -1;
    }

    if (n <= 0 || k < 0 || k > 4 || m <= 0 || p < 0)
    {
        std::cout << "\n\nEXCEPTION: check_data: n, m, k or p are incorrect\n\n\n";
        return -1;
    }

    long int t_full;
    double *A = (double *)malloc(n*n*sizeof(double));
    double *x = (double *)malloc(n*sizeof(double));
    for (int i = 0; i < n; i++) { x[i] = 0; }

    double *b = (double *)malloc(n*sizeof(double));
    double *solution = (double *)malloc(n*sizeof(double));
    make_solution(n, solution);

    double *A_copy = (double *)malloc(n*n*sizeof(double));
    double *b_copy = (double *)malloc(n*sizeof(double));
    int *variable_order = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++)
    {
        variable_order[i] =  n-1-i;
    }

    pthread_t *threads = (pthread_t*)malloc(p * sizeof(pthread_t));
    for (int i = 0; i < p; i++)
    {
        threads[i] = i;
    }
    double *residual = (double *)malloc(p * sizeof(double));
    for (int i = 0; i < p; i++)
    {
        residual[i] = 0;
    }
    int *control = (int *)malloc(p*sizeof(int));
    for (int i = 0; i < p; i++)
    {
        control[i] = 0;
    }        // для ошибок в потоках
    Args *args = (Args *)malloc(p * sizeof(Args));
    R_args *r_args = (R_args *)malloc(p * sizeof(R_args));
    double *time = (double *)malloc(p*sizeof(double));
    std::cout << "Made vectors" << endl;
    if (!A || !x || !b || !solution || !A_copy || !b_copy || !variable_order || !threads || !residual || !control || !args || !r_args || !time)
    {
        std::cout << "EXCEPTION: memory error" << endl;
        free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
        return -1;
    }
    if (argc == 6)
    {
        string input(argv[5]);
        int c1 = read(n, input, A);
        if (c1 != 0)
        {
            MyException C1(c1);
            C1.GetMessage(c1);
            free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
            return c1;
        }
    }
    else
    {
        generate(n, k, A);
    }
    std::cout << "\n\nINITIAL MATRIX:" << endl;
    make_b(n, A, b);
    print_expanded(m, n, n, A, b);
    for (int i = 0; i < n; i++)
    {
        b_copy[i] = b[i];
    }
    for (int i = 0; i < n*n; i++)
    {
        A_copy[i] = A[i];
    }


    std::cout << "finished preparing, started solution" << endl;

    int thread_create_flag = 0;
    int flag1 = 0;
    int flag2;

    t_full = get_full_time();
    for (int i = 0; i < p; i++)
    {
        //args[i].fill(n, p, i, A_copy, b_copy, x, variable_order, &thread_create_flag, control);
         args[i].n = n;
         args[i].p = p;
         args[i].cur_p = i;
         args[i].A_copy = A_copy;
         args[i].b_copy = b_copy;
         args[i].x = x;
         args[i].variable_order = variable_order;
         args[i].thread_create_flag = &thread_create_flag;
         args[i].control = control;

        flag1 = pthread_create(&threads[i], nullptr, solve_system_treaded, &args[i]); // try to make thread, if create != 0 -- idite nahuy
        // put pthread_create into if
        if (flag1 != 0)
        {
            thread_create_flag = -1;
            std::cout << "\n\nEXCEPTION: cannot create thread: " << i << "\n\n\n";
            for (int j = 0; j < i; j++)
            {
                flag2 = pthread_join(threads[i], nullptr);
                if (flag2 != 0)
                {
                    std::cout << "\n\nEXCEPTION: cannot join threads " << j << "\n\n\n"; 
                    free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
                    return -9;
                }
            }
            free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
            return -8;
        }
    }

    thread_create_flag = 1; // successfully created all threads
    for (int i = 0; i < p; i++)
    {
        flag1 = pthread_join(threads[i], nullptr);  // not joining? fuck it
        if (flag1 != 0)
        {
            std::cout << "\n\nEXCEPTION: cannot join threads after joining:" << i << "\n\n\n";
            free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
            return -10;
        }
    }
    t_full = get_full_time() - t_full;
    //double total = time[0];
    for (int i = 0; i < n; i++)
    {
        std::cout << x[i] << " ";
    }
    std::cout << "printed x" << endl << endl << endl << endl;

    //print_expanded(m, n, n, A_copy, b_copy);
    std::cout << "printed" << endl;


//    if (A) free(A);
//    if (b) free(b);
//    if (A_copy) free(A_copy);
//    if (b_copy) free(b_copy);
//    if (variable_order) free(variable_order);
//    if (x) free(x);
//    if (solution) free(solution);
//    if (threads) free(threads);
//    if (residual) free(residual);
//    if (control) free(control);
//    if (args) free(args);
//    if (r_args) free(r_args);
//    if (time) free(time);
//    std::std::cout << endl << "SUCCESS" << endl << endl << endl;
//    return 0;


    double elapsed = time[0];
    elapsed = t_full;


    std::cout << "Made solution, starting residual" << endl;
    //make_solution(n, solution);     ////////////
    if (threads)        // new task(residual_norm) -- new threads
        free(threads);
    
    pthread_t *threads2 = (pthread_t *)malloc(p*sizeof(pthread_t));
    if(!threads2)
    {
        std::cout << "\n\nEXCEPTION: memory error for new threads in residual norm\n\n\n";
        free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
        return -11;
    }
    for (int i = 0; i < p; i++)
    {
        threads2[i] = i;
    }

    thread_create_flag = 0;

    for (int i = 0; i < p; i++)
    {
        r_args[i].fill(n, p, i, A, b, x, residual, &thread_create_flag);
        // r_args[i].n = n;
        // r_args[i].p = p;
        // r_args[i].cur_p = i;
        // r_args[i].A = A;
        // r_args[i].b = b;
        // r_args[i].x = x;

        flag1 = pthread_create(&threads2[i], nullptr, residual_norm_threaded, &r_args[i]);
        if (flag1 != 0)
        {
            thread_create_flag = -1;
            std::cout << "\n\n EXCEPTION: cannot create thread " << i << "\n\n\n";
            for (int j = 0; j < i; j++)
            {
                flag2 = pthread_join(threads2[i], nullptr);
                if (flag2 != 0)
                {
                    std::cout << "\n\n EXCEPTION: cannot join threads " << j << "\n\n\n";
                    free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
                    return -13;
                }
            }
            free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
            return -12;
        }
        // if (control[i] != 0)
        // {
        //     MyException C3(control[i]);
        //     C3.GetMessage(control[i]);
        //     return control[i];
        // }
    }
    thread_create_flag = 1;
    for (int i = 0; i < p; i++)
    {
        flag1 = pthread_join(threads2[i], nullptr);
        if (flag1 != 0)
        {
            std::cout << "\n\n EXCEPTION: cannot join threads after successfully creating them in residual norm: " << i << "\n\n\n";
            free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
            return -14;
        }
    }

    for (int i = 0; i < p; i++)     // checking for any errors in the method
    {
        if (control[i] != 0)
        {
            MyException e(control[i]);
            e.GetMessage(control[i]);
            free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
            return -1;
        }
    }

    //double residual_norm = residual[0];

    //std::cout << "Residual norm: " << residual_norm << endl <<  "Elapsed = " << total << endl << "k = " << k << endl << "m = " << m << endl << "n = " << n << endl << "p = " << p << endl;
    printf("%s : residual = %e, elapsed = %.2f s = %d n = %d p = %d\n", argv[0], residual, elapsed, k, n, m, p);



    // -------------------------------------------------------------------------------------------------


    // auto start = chrono::high_resolution_clock::now();
    // //double time = clock();
    // int c2 = solve_system(A_copy, b_copy, x, n, variable_order);
    // if (c2 != 0)
    // {
    //     MyException C2(c2);
    //     C2.GetMessage(c2);
    //     free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, args, r_args, time);
    //     return c2;
    // }
    // //time = (clock() - time) / CLOCKS_PER_SEC;
    // auto stop = chrono::high_resolution_clock::now();
    // auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    // //std::cout << endl << "TIME: " << time << endl << endl << endl;
    // std::cout << endl << "\nTIME: " << duration.count()/10 << endcout << "" << endl;l << endl << endl;
    // std::cout << "\nSolution:\n";
    // print(m, n, 1, x);
    // std::cout << endl;
    // // double *solution = (double *)malloc(n*sizeof(double));


    // --------------------------------------------------------------------------------------------------



//    std::cout << "\nError norm: ";
//    std::cout << error_norm(n, x, solution) << endl;
//    std::cout << "\nResidual norm: ";
//    std::cout << residual_norm(n, A, b, x) << endl;
    //free_all(A, b, A_copy, b_copy, solution, x, variable_order, threads, residual, control, args, r_args, time);
//    free(A);
//    free(b);
//    free(A_copy);
//    free(b_copy);
//    free(solution);
//    free(x);
//    free(variable_order);
//    free(threads);
//    free(residual);
//    free(control);
//    free(args);
//    free(r_args);
//    free(time);
//    for (int i = 0; i < p; i++)
//    {
//        residual[i] = 0;
//    }

    if (A){ std::cout << "A" << endl; free(A); }
    if (b){ std::cout << "b" << endl; free(b); }
    if (A_copy) { std::cout << "A_copy" << endl; free(A_copy); }
    if (b_copy) { std::cout << "b_copy" << endl; free(b_copy); }
    if (variable_order) { std::cout << "var_ord" << endl; free(variable_order); }
    if (x) { std::cout << "x" << endl; free(x); }
    if (solution) { std::cout << "solution" << endl; free(solution); }
    if (threads) { std::cout << "threads" << endl; free(threads); }
    if (residual) { free(residual); std::cout << "residual" << endl; }
    if (control) { std::cout << "control" << endl; free(control); }
    if (args) { std::cout << "args" << endl; free(args); }
    if (r_args) { std::cout << "r_args" << endl; free(r_args); }
    if (time) { std::cout << "time" << endl; free(time); }
    std::cout << endl << "SUCCESS" << endl << endl << endl;
    return 0;
}