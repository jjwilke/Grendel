
#ifndef _gigide_printstream_h_
#define _gigide_printstream_h_

#include <sstream>
#include <stdio.h>
#include <iostream>
#include <map>

namespace std {

string 
stream_printf(const char* fmt, ...);

template <
  class Type
> void
print_keys(
    map<string, Type>& keymap,
    ostream& os = cout
)
{
    typename map<string, Type>::iterator it;
    for (it = keymap.begin(); it != keymap.end(); it++)
        os << it->first << endl;
}

void
print_map(
    map<string, string>& keymap,
    ostream& os = cout
);

}

#endif
