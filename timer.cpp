#include "timer.hpp"
// variable initialization
Timer::Timer(){
    program_start = time(NULL);
    program_time = 0;
}

void Timer::PrintProgress( long curstep, long total ) {
    if(curstep == 0) { start_time = time(NULL); 
    curstep++; }
    double done = (double)(curstep)/(double)total;
    if(curstep == 1) { curstep--; }
    time_t curtime = time(NULL);
    time_t elapsed = curtime-start_time;
    time_t remain = (elapsed/done*(1.-done));
    long sec;
    long hrs = remain/3600;
    sec = remain%3600;
    long min = sec/60;
    sec = sec%60;

    std::cout << "\rstep " << std::fixed << std::setw(10) << curstep << "/" << total << " (" << std::setw(4) << std::right << std::setprecision(2) << done*100 << "%) - Remaining time: " << convertSecToString(remain) << std::flush;
}

void Timer::endProgram() {
    program_time = time(NULL) - program_start;
}

std::string Timer::curRunTime() {
    return convertSecToString(time(NULL) - program_start);
}

Time convertSecToTime( time_t sec ) {
    Time t;
    t.hrs = static_cast<long>(sec)/3600;
    t.secs= static_cast<long>(sec)%3600;
    t.mins= t.secs/60;
    t.secs= t.secs%60;
    return t;
}

std::string convertTimeToString( Time t ) {
    std::stringstream out;
    out << std::setw(5) << std::right << t.hrs << " hrs " << std::setw(2) << std::right << t.mins << " min " << std::setw(2) << std::right << t.secs << " sec.";
    return out.str();
}

std::string convertSecToString( time_t sec ) {
    return convertTimeToString( convertSecToTime( sec ) );
}