
#ifndef _gigide_printstream_h_
#define _gigide_printstream_h_

#include <sstream>
#include <stdio.h>
#include <iostream>
#include <map>

namespace std {

std::string 
stream_printf(const char* fmt, ...);

template <
  class Type
> void
print_keys(
    std::map<std::string, Type>& keymap,
    std::stringstream& sstr
)
{
    typename std::map<std::string, Type>::iterator it;
    for (it = keymap.begin(); it != keymap.end(); it++)
        sstr << it->first << std::endl;
}

void
print_map(
    std::map<std::string, std::string>& keymap,
    std::stringstream& sstr
);

}

#endif
