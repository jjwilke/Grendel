#ifndef _gigide_gigtimer_h_
#define _gigide_gigtimer_h_

#include <string>
#include <map>
#include <iostream>
#include <utility>
#include <printstream.h>

namespace gigtimer {

class Timer {

    typedef std::string timestr;
    typedef double time_t;

    private:
        typedef std::pair<time_t, unsigned int> time_count_t;
        typedef std::map<timestr, time_count_t> active_map;
        typedef std::map<timestr, time_t> total_map;
        static active_map active_;
        static total_map totals_;
    
    public:
        
        static void start(const timestr& region);

        static void stop(const timestr& region);

        static time_t getTime();

        static void getDisplay(time_t total, unsigned int& hours, unsigned int& minutes, time_t& seconds);

        static void print(std::ostream& os = std::cout);
        
};

}

#endif

