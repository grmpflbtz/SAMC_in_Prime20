// timer class to calculate remaining time
#ifndef _timer_HPP
#define _timer_HPP
#include <ctime>
#include <iostream>
#include <iomanip>
#include <string>
class Timer{
    private:
        time_t start_time;
        time_t program_start;
        time_t program_time;

    public:
        void PrintProgress( long curstep, long total );             // print progress and remaining time estimate
        void endProgram();
        std::string curRunTime();

        Timer();
        ~Timer(){};
};

struct Time{
    long hrs;
    long mins;
    long secs;
};

Time convertSecToTime( time_t sec );
std::string convertTimeToString( Time t );
std::string convertSecToString( time_t sec );

#endif