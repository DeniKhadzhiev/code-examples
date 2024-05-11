#include <iostream>
#include <map>
#include <string>

using namespace std;

class MyException
{
private:
    int code;
    const string message;
public:
    MyException() : code(0), message("OK") {};
    MyException(int code_, const string message_) : code(code_), message(message_) {};
    MyException(int code_) : code(code_) {};
    ~MyException() {};

    void GetMessage(int code)
    {
        string messages[16] = {
            "\n\nOK\n\n\n",

            "/n/n EXCEPTION: read: Not able to read matrix /n/n/n",
            "/n/n EXCEPTION: read: Unable to open file for reading matrix A /n/n/n",
            "/n/n EXCEPTION: read: wrong format of file, found letters or special symbols /n/n/n",
            "/n/n EXCEPTION: read: Not enough numbers /n/n/n",

            "/n/n EXCEPTION: normalize_matrix: A is singular /n/n/n",

            "/n/n EXCEPTION: solve_system: Matrix is singular 2 /n/n/n",
            "/n/n EXCEPTION: solve_system: matrix is singular 3 /n/n/n",
            "/n/n EXCEPTION: cannot create thread: /n/n/n",
            "/n/n EXCEPTION: cannot join threads: /n/n/n",
            "/n/n EXCEPTION: cannot join threads after joining: /n/n/n",
            "/n/n EXCEPTION: memory error for new threads in residual norm /n/n/n",
            "/n/n EXCEPTION: solve_system: did not create threads to use them /n/n/n",
            "/n/n EXCEPTION: solve_system: n = 1, matrix is singular /n/n/n",
            "/n/n  /n/n/n",
            "/n/n  /n/n/n"
        };
        //map<int, string> messages //= {
            // {0, "\n\nOK\n\n\n"},

            // {-1, "/n/n EXCEPTION: read: Not able to read matrix /n/n/n"},
            // {-2, "/n/n EXCEPTION: read: Unable to open file for reading matrix A /n/n/n"},
            // {-3, "/n/n EXCEPTION: read: wrong format of file, found letters or special symbols /n/n/n"},
            // {-4, "/n/n EXCEPTION: read: Not enough numbers /n/n/n"},

            // {-5, "/n/n EXCEPTION: normalize_matrix: A is singular /n/n/n"},

            // {-6, "/n/n EXCEPTION: solve_system: Matrix is singular 2 /n/n/n"},
            // {-7, "/n/n EXCEPTION: solve_system: matrix is singular 3 /n/n/n"},
            // {-8, "/n/n  /n/n/n"},
            // {-9, "/n/n  /n/n/n"},
            // {-10, "/n/n  /n/n/n"},
            // {-11, "/n/n  /n/n/n"},
            // {-12, "/n/n  /n/n/n"},
            // {-13, "/n/n  /n/n/n"},
            // {-14, "/n/n  /n/n/n"},
            // {-15, "/n/n  /n/n/n"}


        // map[0] = "\n\nOK\n\n\n";
        // map[-1] = "\n\n  \n\n\n";
        // map[-2] = "\n\n  \n\n\n";
        // map[-3] = "\n\n  \n\n\n";
        // map[-4] = "\n\n  \n\n\n";
        // map[-5] = "\n\n  \n\n\n";
        // map[-6] = "\n\n  \n\n\n";
        // map[-7] = "\n\n  \n\n\n";
        // map[-8] = "\n\n  \n\n\n";
        // map[-9] = "\n\n  \n\n\n";
        // map[-10] = "\n\n  \n\n\n";
        // map[-11] = "\n\n  \n\n\n";
        // map[-12] = "\n\n  \n\n\n";
        // map[-13] = "\n\n  \n\n\n";
        // map[-14] = "\n\n  \n\n\n";


        // };

        // auto it = messages.find(code);
        std::cout << "CODE: " << code << "\nMESSAGE:   " << messages[-code];
    }
};