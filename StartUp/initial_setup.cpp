//==================================================================================================
//
// simple setup routines.
//
//==================================================================================================
//
// Author: Elizabeth Lee
//       : Based on work by Jens Chluba
// first implementation: November 2018
// last modification: November 2018
//
//==================================================================================================

//==================================================================================================
void read_startup_data_parser(string filename, Parameters &inp, int wait=0)
{//Must be initialised if given
    //==============================================================================================
    // read values (if given) from file
    //==============================================================================================
    inp.SetParametersFromFile(filename);

    if(show_mess>=0) cout << endl;
    
    if(show_mess>=1)
    {//Show read in parameters if set, so that they can be checked by the user, theoretically
        if(filename!="") cout << " Using cosmological parameters corresponding to " << filename << endl;
    }

    if (!inp.CheckValues()){
        print_error("Using default Parameters.");
        inp = Parameters();
    }
            
    if(show_mess>=0)
    {//Let the user know where the file will be saved.
        cout << "\n output path: " << inp.outputPath << endl;
        cout << " file ending  : " << inp.fileEnding << endl;
        cout << " run mode : " << inp.rare.RunMode << endl;
    }

    if(show_mess>=0) cout << endl;
    
    if(wait==1) wait_f_r();

    return;
}

//==================================================================================================
//==================================================================================================
