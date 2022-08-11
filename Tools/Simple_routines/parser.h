//====================================================================================================================
//
// simple parser written by Boris Bolliet and Jens Chluba (Feb 2018). This was motivated by a parser from Class.
// Slightly modified by Elizabeth Lee (Nov 2018)
//
// Purpose: read parameters with given identifier from parameter file. It is assumed that the entry is written in the
//          form 'id-string = entry'. Lines starting with '#' and empty lines are omitted. For entries of the form
//          'id-string = entry entry2' the read functions currently only account for the first entry.
//
//====================================================================================================================
// 23.04.2018: overcame problem with variable names appearing after '#' with was not at beginning of line

#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <vector>
#include <map>
#include <regex>
#include "global_functions.h"

using namespace std;

struct file_content
{
    vector<string> lines;
    string filename;
};

//====================================================================================================================
void parser_read_file(string filename, struct file_content &pfc, bool show_lines=false);
void parser_free(struct file_content &pfc);

// generic read of types like int, long int, double etc
template <class T>
bool parser_read(const struct file_content &pfc, string var_id, T &val);

bool parser_read_int(const struct file_content &pfc, string var_id, int &val);
bool parser_read_double(const struct file_content &pfc, string var_id, double &val);
bool parser_read_string(const struct file_content &pfc, string var_id, string &val);

bool parser_read_3vector(const struct file_content &pfc, string var_id, vector<double> &val);

#endif

//====================================================================================================================
//====================================================================================================================
