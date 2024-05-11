#include <iostream>
#include <sys/time.h> 
#include <sys/resource.h> 
#include <unistd.h>
#include "get_time.hpp"


long int get_time()
{
    struct rusage buf;
    getrusage(RUSAGE_SELF, &buf);

    return buf.ru_utime.tv_sec * 10 + buf.ru_utime.tv_usec / 10000;
}

long int get_full_time ()
{
    struct timeval buf;
    gettimeofday (&buf, 0);
    return buf. tv_sec * 10 + buf.tv_usec / 10000;

}