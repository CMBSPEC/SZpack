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

//====================================================================================================================
// Standards
//====================================================================================================================
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include "parser.h"

using namespace std;

//====================================================================================================================
void parser_read_file(string filename, struct file_content &pfc, bool show_lines)
{
    ifstream ifile(filename.c_str());
    
    if(!ifile) {
        exit_error("parser_read_file :: Error opening parameter file " + filename + ".");
    }

    pfc.filename=filename;
    
    do {
        string line="";
        getline(ifile, line);
        // drop comment lines
        if(line[0]!='#' && line.length()>1)
        {
            pfc.lines.push_back(line);
            // erase things after possible second '#'
            size_t pos=line.find('#');
            if(pos==string::npos) continue;
            pfc.lines.back().erase(pos, line.length());
        }
    }
    while(!ifile.eof());
    
    ifile.close();
    
    if(show_lines) for(int l=0; l<(int)pfc.lines.size(); l++) cout << pfc.lines[l] << endl;
}

//====================================================================================================================
void parser_free(struct file_content &pfc)
{
    pfc.lines.clear();
}

//====================================================================================================================
template <class T>
bool parser_read(const struct file_content &pfc, string var_id, T &val)
{
    for(int l=0; l<(int)pfc.lines.size(); l++)
    {
        string str=pfc.lines[l];
        size_t pos=str.find(var_id);
        
        // continue if at end of string
        if(pos>=str.length()) continue;
     
        // if string is found, continue search from positions on
        pos=str.find("=", pos)+1;
        for(; pos<str.length(); pos++) if(str[pos]!=' ') break;
            
        // now at position of entry and need to convert it (everything after entry is omitted)
        istringstream iss(str.substr(pos));
        iss >> val;
        
        return true;
    }
    return false;
}

bool parser_read_int(const struct file_content &pfc, string var_id, int &val){
    return parser_read(pfc, var_id, val);
}
bool parser_read_double(const struct file_content &pfc, string var_id, double &val){
    return parser_read(pfc, var_id, val);
}
bool parser_read_string(const struct file_content &pfc, string var_id, string &val){
    return parser_read(pfc, var_id, val);
}

bool parser_read_3vector(const struct file_content &pfc, string var_id, vector<double> &val){
    if (val.size() < 3) return false;

    for(int l=0; l<(int)pfc.lines.size(); l++)
    {
        string str=pfc.lines[l];
        size_t pos=str.find(var_id);
        
        // continue if at end of string
        if(pos>=str.length()) continue;
     
        // if string is found, continue search from positions on
        pos=str.find("=", pos)+1;
        for(; pos<str.length(); pos++) if(str[pos]!=' ') break;
        
        // now at position of entry and need to convert it (everything after entry is omitted)
        istringstream iss(str.substr(pos));
        iss >> val[0] >> val[1] >> val[2];
        
        return true;
    }
    return false;
}

//====================================================================================================================
//====================================================================================================================
