#include "timer.hpp"
#include <map>
#include <string>
#include <memory>
#include <fstream>

extern "C" {
    typedef void* timer_handle;

    timer_handle mytimer_start(const char* name) {
        return new CTimer(std::string(name));
    }

    void mytimer_stop(timer_handle h) {
        delete static_cast<CTimer*>(h);
    }

    void mytimer_print() {
        CTimer::print_timing_results();
    }
    
    void mytimer_gather_and_print() {
        
        std::vector<TimerData> all_timings;
        CTimer::gather_and_print(0, all_timings); // Root is 0
    }
}
