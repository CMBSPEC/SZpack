//==================================================================================================
//
// global functions that are part of the main SZpack setup and runmodes
//
//==================================================================================================

#ifndef globalfunctions_H
#define globalfunctions_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

//==================================================================================================
// level of messages printed while running
//==================================================================================================
static int show_mess = 1;

//==================================================================================================
//
// Messages to Screen
//
//==================================================================================================
static inline void print_message(string mess)
{
    if(show_mess >= 1)
    {
        cout << "\n " << setfill('=') << setw(90) << "=" << endl;  
        cout << " || " << mess << endl;
        cout << " " << setfill('=') << setw(90) << "=" << endl << endl;   
    }
}

//==================================================================================================
#ifdef UsePybind11
#include <pybind11/pybind11.h>
#include <stdexcept>
namespace py = pybind11;
static inline void print_error(string mess)
{
    py::print("ERROR:",mess);
}

static inline void exit_error(string mess)
{
    throw std::runtime_error(mess);
    exit(1);
}

#else
static inline void print_error(string mess)
{
    cerr << "\n" << setfill('#') << setw(90) << "#" << endl;
    cerr << " # ERROR: " << mess << endl;
    cerr << " " << setfill('#') << setw(90) << "#" << endl << endl;
}

static inline void exit_error(string mess)
{
    print_error(mess);
    exit(1);
}
#endif

//==================================================================================================
//==================================================================================================
#endif
