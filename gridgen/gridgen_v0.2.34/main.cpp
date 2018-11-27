///////////////////////////////////////
// Name:        main.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     13/10/2015
// Last Modification: 31/07/2018 A.Pisarski
///////////////////////////////////////

#include "ReadConfig.h"
#include "OptionValue.h"
#include "num.h"
#include "ReadEphemeris.h"
#include "FisherRM.h"
#include "stringmanip.h"
#include "DensityS1.h"
#include "DensityS2.h"
#include "GridS1.h"
#include "GridS2.h"
#include "Manual.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <cstdlib>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>
//#include <sys/ioctl.h> //linux
//#include <windows.h>

using namespace std;

int method(double c0, unsigned int nfft, unsigned int data_length, unsigned int spindown=1);
bool if_exist(map<string, vector<string> >&, string);


int main(int argc, char* argv[])
{
    try{
        //cout << string(*argv) <<endl;

        vector<string> commands;
        for(int i=1; i<argc; i++)           // convert to lowercase string
        {
            string temp = string(*(argv+i));

            //#mb this ain't any good
            //std::transform(temp.begin(), temp.end(), temp.begin(), ::tolower);
            commands.push_back(temp);
        }

        vector<string> options_available;   // order is important (mast be the same as in long version)
        options_available.push_back("-h");  // help
        options_available.push_back("-ht"); // print help to text file.
        options_available.push_back("-c");  // covariance   == (minimal match)^2
        options_available.push_back("-m");  // minimal match
        options_available.push_back("-d");  // directory
        options_available.push_back("-a");  // algorithm (to choose)
        options_available.push_back("-i");  // initial time
        options_available.push_back("-n");  // number of Fourier bins (nfft length)
        options_available.push_back("-nd"); // number of data
        options_available.push_back("-na"); // number of loops (incrementation deep) in root finding algorithm
        options_available.push_back("-nr"); // number of loops in covering radius (deep hole) finding algorithm
        options_available.push_back("-p");  // print g - grid; f - fisher matrix; d - density of covering; a - all;
                                            // o - print in original (author) mode (without that option program will work,
                                            // with Michal Bejger's sets (result saved to the 'grid.bin' file); ao - this same as: dfgo.
        options_available.push_back("-ch"); // chop (set to zero if number is less than num::epsilon)
        options_available.push_back("-cv"); // convert grid from hyper-sphere to hyper-ellipsoid space
        options_available.push_back("-v");  // program version
        options_available.push_back("-s");  // spindown
        options_available.push_back("-aa"); // show information about author(s)
        vector<string> options_available_long;              // order is important (mast be the same as without _long)
        options_available_long.push_back("--help");         // help
        options_available_long.push_back("--helptxt");      // print help to text file.
        options_available_long.push_back("--covariance");   // covariance   == (minimal match)^2
        options_available_long.push_back("--match");        // minimal match
        options_available_long.push_back("--directory");    // directory
        options_available_long.push_back("--algorithm");    // algorithm (to choose)
        options_available_long.push_back("--initial");      // initial time
        options_available_long.push_back("--nfft");         // number of Fourier bins (nfft length)
        options_available_long.push_back("--ndata");        // number of data
        options_available_long.push_back("--nalpha");       // number of loops (incrementation deep) in root finding algorithm
        options_available_long.push_back("--nradius");      // number of loops in covering radius (deep hole) finding algorithm
        options_available_long.push_back("--print");        // print g - grid; f - fisher matrix; d - density of covering; a - all;
                                                            // o - print in original (author) mode (without that option program will work,
                                                            // with Michal Bejger's sets (result saved to the 'grid.bin' file); ao - this same as: dfgo.
        options_available_long.push_back("--chop");         // chop (set to zero if number is less than num::epsilon)
        options_available_long.push_back("--convert");      // convert grid from hyper-sphere to hyper-ellipsoid space
        options_available_long.push_back("--version");      // program version
        options_available_long.push_back("--spindown");     // spindown
        options_available_long.push_back("--author");       // show information about author(s)

        /// to mark position of encounter options apart negative values
        vector<unsigned int> position;
        for(unsigned int i=0; i<commands.size(); i++)
        {
            string current = commands[i];
            string found = string(1, current.at(0));
            bool not_number = !(current.find_first_of(".0123456789") != std::string::npos);
            if( (found=="-" && not_number) || (current.find("--")!=std::string::npos && not_number) )
                position.push_back(i);
        }

        sort(position.begin(), position.end());

        map<string, vector<string> > options_found;
        for(unsigned int i=0; i<position.size(); i++)
        {
            unsigned int b, e;
            if(position[i]!=commands.size())
            {
                b = position[i];
                if(i<position.size()-1)
                    e = position[i+1];
                else
                    e = commands.size();

                unsigned int k = b+1;
                if(k != e)
                {
                    for(unsigned int j = k; j<e; j++)
                        options_found[commands[b]].push_back(commands[j]);
                }
                else{
                    options_found[commands[b]].push_back(string(""));
                }
            }
        }

        /// if map of long options exist and short options not
        /// create map of short options only (options_found);
        /// this allow to work further with one version of options
        for(unsigned int i=0; i<options_available.size(); ++i)
        {
            if( !if_exist(options_found, options_available[i])
               &&  if_exist(options_found, options_available_long[i]) )
            {
                map<string, vector<string> >::const_iterator it = options_found.find(options_available_long[i]);
                vector<string>::const_iterator line_it = it->second.begin();
                if(*line_it!="")
                {
                    while (line_it != it->second.end())
                    {
                        options_found[options_available[i]].push_back(*line_it);
                        ++line_it;
                    }
                }
                else{
                    options_found[options_available[i]].push_back("");
                }


                options_found.erase(options_available_long[i]);
            }
        }

        /// Convert minimal match to covariance
        if( !if_exist(options_found, "-c") &&  if_exist(options_found, "-m") )
        {
            map<string, vector<string> >::const_iterator it = options_found.find("-m");
            vector<string>::const_iterator line_it = it->second.begin();

            if(*line_it!="")
            {
                while (line_it != it->second.end())
                {
                    double mm = stod(*line_it);
                    double c0 = pow(mm, 2);
                    string cov = to_string(c0);
                    options_found["-c"].push_back(cov);
                    ++line_it;
                }
            }
            else{
                options_found["-c"].push_back("");
            }

            options_found.erase("-m");
        }
        /// Remove unnecessary (now) command
        if( if_exist(options_found, "-m") )
            options_found.erase("-m");

        /// create defaults flags (user don't need provide them)
        for(unsigned int i=0; i<options_available.size(); ++i)
        {
            if( !if_exist(options_found, options_available[i]) )
            {
                if(options_available[i]!="-h" && options_available[i]!="-ht"
                   && options_available[i]!="-v" && options_available[i]!="-aa"
                    && options_available[i]!="-c" && options_available[i]!="-s") // && options_available[i]!="-m"
                        options_found[options_available[i]].push_back("");
            }

        }

        Manual manual("\t * Build 0.2.34 (alpha).                   *");

        double CovarianceMin, CovarianceMax, CovarianceStep=0.01;
        double InitialTimeMin=0.5, InitialTimeMax, InitialTimeStep;
        unsigned int NfftExpLength=20, DataLength=0, Spindown=1, dim = 4;
        int NAlpha=35, NRadius=20;
        string PathSSB, PathDet, PathDetSSB;
        string PathSave;
        string FilePatternGrid, FilePatternFM, FilePatternDoC;
        string SaveGrid="False", SaveFisherMatrix="False", SaveDensityOfCovering="False";
        string DataChop="False";
        string ChooseMethod="Automatic";
        string DataS2Use, DataS2SaveBest;
        string Convert="True";
        string Quiet;
        string SayHello;
        size_t sizeCo, sizeIt, sizeNEL, sizeDLn;
        string snel, sdln;
        InitialTimeMax=InitialTimeMin + num::epsilon();
        InitialTimeStep=1.0;

        string OriginalMode="False";
        string gridout;
        vector<double> M, Mn;

        if( commands.size() != 0) /// if flags in command line are provided
        {
            ///Help
            if( if_exist(options_found, "-h") )
                manual.print_help();

            ///Help print to text file
            if( if_exist(options_found, "-ht") )
                manual.save_help();

            ///Version
            if( if_exist(options_found, "-v") )
                manual.print_version();

            ///About author
            if( if_exist(options_found, "-aa") )
                manual.print_author();


            if( if_exist(options_found, "-c") )
            {
                /// Directory; PathSSB, PathDet, PathDetSSB
                if( if_exist(options_found, "-d") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-d");
                    vector<string>::const_iterator line_it = it->second.begin();

                    PathSave = *line_it;

                    if(*line_it=="")
                    {
                        PathSSB="rSSB.bin";
                        PathDet="rDet.bin";
                        PathDetSSB="DetSSB.bin";
                    }
                    else{
                        PathSSB=*line_it+"rSSB.bin";
                        PathDet=*line_it+"rDet.bin";
                        PathDetSSB=*line_it+"DetSSB.bin";
                    }
                }
                //cout << PathSSB << endl;

                /// Covariance
                {   //local scope
                    map<string, vector<string> >::const_iterator it = options_found.find("-c");
                    vector<string>::const_iterator line_it = it->second.begin();
                    vector<double> temp_d;
                    int k=0;

                    if(*line_it!="")
                    {
                        while (line_it != it->second.end())
                        {
                            double tt = stod(*line_it);
                            if(num::sign(tt)==-1)
                                tt*=-1;

                            temp_d.push_back(tt);
                            ++line_it; k++;
                        }
                    }
                    else{
                        CovarianceMin=0.75;
                        CovarianceMax=CovarianceMin;
                        CovarianceStep=0.01;
                    }

                    //cout << "k(c)= " << k << endl;
                    if(k>=2)
                    {
                        sort(temp_d.begin(), temp_d.begin()+1);
                        CovarianceMin=temp_d[0];
                        CovarianceMax=temp_d[1];

                        if(k==3)
                            CovarianceStep=temp_d[2];
                        else   // k==2
                            CovarianceStep=0.01;
                    }

                    if(k==1)
                    {
                        CovarianceMin=temp_d[0];
                        CovarianceMax=CovarianceMin;
                        CovarianceStep=1.0;
                    }

                    CovarianceMax+=num::epsilon();
                    /*
                    if(*(it->second.begin())!="")
                    {
                        cout << "CovarianceMin= " << CovarianceMin << endl;
                        cout << setw(13) << setprecision(12) << "CovarianceMax= " << CovarianceMax << endl;
                        cout << "CovarianceStep= " << CovarianceStep << endl;
                    }
                    */

                }

                ///InitialTime
                if( if_exist(options_found, "-i") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-i");
                    vector<string>::const_iterator line_it = it->second.begin();
                    vector<double> temp_d;
                    int k=0;

                    if(*line_it!="")
                    {

                        while (line_it != it->second.end())
                        {
                            double tt = stod(*line_it);

                            temp_d.push_back(tt);
                            ++line_it; k++;
                        }
                    }
                    else{
                        InitialTimeMin=0.5;
                        InitialTimeMax=InitialTimeMin;
                        InitialTimeStep=1.0;
                    }

                    //cout << "k(i)= " << k << endl;
                    if(k>=2)
                    {
                        sort(temp_d.begin(), temp_d.begin()+1);
                        InitialTimeMin=temp_d[0];
                        InitialTimeMax=temp_d[1];

                        if(k==3)
                        {
                            InitialTimeStep=temp_d[2];
                            if(num::sign(InitialTimeStep)==-1)
                                InitialTimeStep*=-1;
                        }
                        else{   // k==2
                            InitialTimeStep=InitialTimeMax-InitialTimeMin;
                            if(abs(InitialTimeStep)<num::epsilon())
                                InitialTimeStep=1.0;
                        }
                    }

                    if(k==1)
                    {
                        InitialTimeMin=temp_d[0];
                        InitialTimeMax=InitialTimeMin;
                        InitialTimeStep=1.0;
                    }
                    InitialTimeMax+=num::epsilon();
                    /*
                    if(*(it->second.begin())!="")
                    {
                        cout << "InitialTimeMin= " << InitialTimeMin << endl;
                        cout << "InitialTimeMax= " << InitialTimeMax << endl;
                        cout << "InitialTimeStep= " << InitialTimeStep << endl;
                    }
                    */
                }

                ///Spindown
                if( if_exist(options_found, "-s") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-s");
                    vector<string>::const_iterator line_it = it->second.begin();
                    if(*line_it=="")
                        Spindown=1;
                    else{
                        unsigned int tt = stoi(*line_it);
                        if(tt == 1 || tt == 2)
                            Spindown = tt;
                        else{
                            string error = "Allowed spindown values: {1, 2}.\n";
                            throw domain_error(error);
                        }
                    }
                    dim = Spindown +3;
                }

                ///ChooseMethod
                if( if_exist(options_found, "-a") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-a");
                    vector<string>::const_iterator line_it = it->second.begin();

                    if(*line_it=="" || *line_it=="a")
                        ChooseMethod="Automatic";
                    else
                        ChooseMethod=*line_it;
                }

                ///Nfft
                if( if_exist(options_found, "-n") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-n");
                    vector<string>::const_iterator line_it = it->second.begin();

                    if(*line_it=="")
                        NfftExpLength=20;
                    else
                    {
                        NfftExpLength=stoi(*line_it);
                        if(NfftExpLength<=0)
                            NfftExpLength=20;
                    }
                }

                ///DataLength
                if( if_exist(options_found, "-nd") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-nd");
                    vector<string>::const_iterator line_it = it->second.begin();

                    if(*line_it=="")
                        DataLength=0;
                    else{
                        int tt = stoi(*line_it);
                        if(tt >= 0)
                            DataLength = static_cast<unsigned int> (tt);
                        else{
                            string error = "Data length can be set on"+*line_it+"\n";
                            throw domain_error(error);
                        }
                    }
                }

                ///NAlpha
                if( if_exist(options_found, "-na") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-na");
                    vector<string>::const_iterator line_it = it->second.begin();
                    if(*line_it=="")
                        NAlpha=35;
                    else{
                        int tt = stoi(*line_it);
                        if(tt > 0)
                            NAlpha = static_cast<unsigned int> (tt);
                        else{
                            string error = "NAlpha can not be less or equal zero!\n";
                            throw domain_error(error);
                        }
                    }
                }

                ///NRadius
                if( if_exist(options_found, "-nr") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-nr");
                    vector<string>::const_iterator line_it = it->second.begin();

                    if(*line_it=="")
                        NRadius=20;
                    else{
                        int tt = stoi(*line_it);
                        if(tt > 0)
                            NRadius = static_cast<unsigned int> (tt);
                        else{
                            string error = "NRadius can not be less or equal zero!\n";
                            throw domain_error(error);
                        }
                    }
                }

                ///Convert
                if( if_exist(options_found, "-cv") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-cv");
                    vector<string>::const_iterator line_it = it->second.begin();

                    if(*line_it=="" || *line_it=="t")
                        Convert="True";
                    else{
                        if(*line_it=="f")
                            Convert= "False";
                        else{
                            string error = "Arguments for -cv (--convert) can be only: t (true) or f (false).\n";
                            throw domain_error(error);
                        }
                    }
                }

                ///DataChop
                if( if_exist(options_found, "-ch") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-ch");
                    vector<string>::const_iterator line_it = it->second.begin();

                    if(*line_it=="" || *line_it=="f")
                        DataChop="False";
                    else{
                        if(*line_it=="t")
                            DataChop= "True";
                        else{
                            string error = "Arguments for -ch (--chop) can be only: t (true) or f (false).\n";
                            throw domain_error(error);
                        }
                    }
                }

                ///Print; SaveGrid, SaveFisherMatrix, SaveDensityOfCovering;
                if( if_exist(options_found, "-p") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-p");
                    vector<string>::const_iterator line_it = it->second.begin();

                    bool no_match = true;
                    size_t ts_end = (*line_it).size();
                    if(ts_end > 4)
                        ts_end = 4;

                    if(*line_it==""){
                        SaveGrid="True"; no_match = false;
                    }
                    else{
                        string temp = string((*line_it).begin(),(*line_it).begin()+ts_end); //[,)

                        if(temp.find_first_of("o")!=string::npos){
                            OriginalMode="True"; no_match = false;
                        }

                        if(temp.find_first_of("a")!=string::npos){
                            SaveDensityOfCovering="True";
                            SaveFisherMatrix="True";
                            SaveGrid="True";
                            no_match = false;
                        }
                        else{
                            if(temp.find_first_of("d")!=string::npos){
                                SaveDensityOfCovering="True"; no_match = false;
                            }

                            if(temp.find_first_of("f")!=string::npos){
                                SaveFisherMatrix="True"; no_match = false;
                            }

                            if(temp.find_first_of("g")!=string::npos){
                                SaveGrid="True"; no_match = false;
                            }
                        }
                    }

                    if(no_match){
                        string error = "Arguments for -p (--print) can be only: d (density of covering),\n"
                        "f (Fisher reduced matrix), g (grid), a (all: density, Fisher matrix, grid),\n"
                        "t - print on terminal.\n";
                        throw domain_error(error);
                    }
                }

                DataS2Use = "True";
                DataS2SaveBest = "True";
/*
        for (map<string, vector<string> >::const_iterator it = options_found.begin(); it != options_found.end(); ++it)
        {
            cout << "[" << it->first << "] have option(s): ";

            vector<string>::const_iterator line_it = it->second.begin();
            cout << *line_it;

            ++line_it;
            while (line_it != it->second.end())
            {
                cout << ", " << *line_it;
                ++line_it;
            }

            cout << endl;
        }
*/
                bool s2use = false;
                if (DataS2Use == "True") s2use = true;
                bool s2save = false;
                if (DataS2SaveBest == "True") s2save = true;
                bool d_chop = false;
                if (DataChop == "True") d_chop = true;

                bool conv = false;
                if (Convert == "True") conv = true;

                unsigned int nfft = pow(2, NfftExpLength);

                //cout << "DataLength= " << DataLength << endl;
                if(DataLength==0)    // DataLength == 0 means DataLength should be equal to ephemeris length
                {
                    //#mb binary file grid bin
                    gridout = PathSave + "grid.bin";
                    std::ofstream gridbin(gridout, std::ios::binary);
                    if(OriginalMode == "False")
                        cout << "grid.bin will be saved in " << gridout << endl;

                    vector<string> paths_to_ephemeris;
                    paths_to_ephemeris.push_back(PathSSB);
                    paths_to_ephemeris.push_back(PathDet);
                    paths_to_ephemeris.push_back(PathDetSSB);

                    ReadEphemeris RE = ReadEphemeris(paths_to_ephemeris);
                    FisherRM const *const cp_FM = new FisherRM(RE.get_ephemeris1(), RE.get_ephemeris2(), Spindown);

                    GridS1 gs1 = GridS1(cp_FM, Spindown);
                    GridS2 gs2 = GridS2(cp_FM, NAlpha, NRadius);

                    unsigned int data_length = cp_FM->get_ephemeris_length();
                    //ofstream fsd;

                    if (SaveDensityOfCovering == "True")
                    {
                        if( cout.good() == true )
                        {
                            if(OriginalMode == "True")
                            	cout << "Covariance, Density of covering\n";
                        }
                        else{
                            string error = "Can't open stream!\n";
                            throw runtime_error(error);
                        }
                    }

                    if (SaveGrid == "True" || SaveDensityOfCovering == "True")
                    for(double c0=CovarianceMin; c0<=CovarianceMax; c0+=CovarianceStep)
                    {
                        //cout << "[c0=" << c0 <<"]\n";
                        int ch_method = 0;
                        if(ChooseMethod == "s1")
                            ch_method = 1;

                        if(ChooseMethod == "s2")
                        {
                            switch(Spindown)
                            {
                                case 1:
                                    ch_method = 2;
                                    break;
                                case 2:
                                    cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                                    cerr << "'s1' will be applied instead." ;
                                    ch_method = 1;
                                    break;
                            }

                        }

                        if(ChooseMethod == "Automatic")
                            ch_method = method(c0, nfft, RE.get_length(), Spindown);

                        if(ch_method == 0)
                        {
                            string error = "Bad option value: " + ChooseMethod + ".\n"
                            "Available options: s1, s2, a; 's2' only for spindown = 1.\n";
                            throw domain_error(error);
                        }
                        ///cout << "test 1" << endl;
                        double density;
                        vector<double> sgrid;
                        if(ch_method == 1)
                        {
                            sgrid = gs1.grid_prim(c0, nfft, data_length);
                            density = gs1.DensityS1::density(sgrid);
                        }
                        ///cout << "test 2" << endl;
                        if(ch_method == 2)
                        {
                            sgrid = gs2.grid_prim(c0, nfft, s2use, s2save);
                            density = gs2.density(sgrid);
                        }

                        if(OriginalMode == "True"){
                            if (SaveDensityOfCovering == "True")
                            {
                                cout.fill('0');
                                cout.width(4);
                                cout << left << c0 << " " << setw(13) << setprecision(12)
                                << density << "\n";
                            }

                            if (SaveGrid == "True")
                            {
                                for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                                {
                                    vector<double> sgrid2;
                                    if( conv )
                                    {
                                        if( ch_method == 1)
                                           sgrid2 = gs1.convert(c0, xi, sgrid);

                                        if( ch_method == 2)
                                            sgrid2 = gs2.convert(c0, xi, sgrid);
                                    }
                                    else{
                                        sgrid2 = sgrid;
                                    }

                                    if(d_chop)
                                        num::chop(sgrid2);

                                    //fstream fs;
                                    //fs.open( PathSave + "\/" + name1, fstream::out | fstream::trunc);
                                    if( cout.good() == true )
                                    {   //fixed << left
                                        cout << setw(15) << setprecision(12) << left;
                                        cout << density << "\n";
                                        for(unsigned int i = 0; i<sgrid2.size(); i++)
                                        {
                                            cout << setw(15) << setprecision(11) << left << sgrid2[i];
                                            if( (i+1)%dim==0 )
                                                     cout << endl;
                                            else
                                                cout << "\t";
                                        }

                                        //        //fs.unget();
                                        //fs.close();
                                    }
                                    else
                                    {
                                        string error = "Can't open stream!\n";
                                        throw runtime_error(error);
                                    }
                                    //cout << name1 << endl;
                                }
                                cout << endl;
                            } // End of 'if(SaveGrid == "True")'
                        }
                        else{
                            //#mb
                            if (SaveDensityOfCovering == "True")
                            {
                                cout.fill('0');
                                cout.width(4);
                                cout << "Density of covering: " << left << c0 << " " << setw(13) << setprecision(12)
                                << density << "\n";
                            }


                                for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                                {
                                    vector<double> sgrid2;
                                    if( conv )
                                    {
                                        if( ch_method == 1)
                                           sgrid2 = gs1.convert(c0, xi, sgrid);

                                        if( ch_method == 2)
                                            sgrid2 = gs2.convert(c0, xi, sgrid);
                                    }
                                    else{
                                        sgrid2 = sgrid;
                                    }

                                    if(d_chop)
                                        num::chop(sgrid2);


                                    //#mb writing the fftpad to grid.bin
                                    int fftpad;
                                    fftpad = (log2(nfft) - ceil(log2(cp_FM->get_ephemeris_length())) + 1);

                                    cout << "fftpad: " << fftpad << endl ;

                                    gridbin.write((char*) &fftpad, sizeof(int));

                                    unsigned int datal = cp_FM->get_ephemeris_length();
                                    /// 'datal' is this same as 'data_length'
                                    M.resize(dim*dim);
                                    Mn.resize(dim*dim);
                                    for(unsigned int i = 0; i<sgrid2.size(); i++) {
                                        M[i] = sgrid2[i]/datal;
                                        Mn[i] = sgrid2[i];
                                    }

                                    // M[1] = M[1]/datal;
                                    // M[5] = M[5]/datal;
                                    // M[9] = M[9]/datal;
                                    // M[13] = M[13]/datal;
                                    for(unsigned int i=0; i<dim; i++)
                                        M[dim*i+1]/=datal;

                                    cout << "Normalized grid matrix:" << endl;
                                    for(unsigned int i = 0; i<M.size(); i++){
                                        cout << scientific << M[i];
                                            if( (i+1)%dim==0 )
                                                cout << "\n";
                                            else
                                                cout << "\t";
                                    }

                                    gridbin.write((char*)&M[0], M.size()*sizeof(double));

                                    if( cout.good() == true )
                                    {   //fixed << left
                                        cout << setw(15) << setprecision(12) << left;
                                        cout << density << "\n";
                                        for(unsigned int i = 0; i<sgrid2.size(); i++)
                                        {
                                            cout << setw(15) << setprecision(11) << left << sgrid2[i];
                                            if( (i+1)%dim==0 )
                                                     cout << endl;
                                            else
                                                cout << "\t";
                                        }

                                        //fs.unget();
                                        //fs.close();
                                    }
                                    else
                                    {
                                        string error = "Can't open stream!\n";
                                        throw runtime_error(error);
                                    }
                                    //cout << name1 << endl;
                                }
                                cout << endl;
                        } // End of else
                    } // End of 'for(double c0=CovarianceMin ...'

                    if(OriginalMode == "True"){
                        if (SaveFisherMatrix == "True")
                        {
                            for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                            {
                                vector<double> fm = cp_FM->postrmf(xi);

                                if( cout.good() == true )
                                {
                                    cout << scientific  << left;//<< setw(16) << setprecision(15) << fixed
                                    for(unsigned int i = 0; i<fm.size(); i++)
                                    {
                                        cout << fm[i];
                                        if( (i+1)%dim==0 )
                                        {
                                            cout << "\n";
                                        }
                                        else
                                        {
                                            cout << "\t";
                                        }
                                    }

                                    //fs.close();
                                }
                                else{
                                    string error = "Can't open stream!\n";
                                    throw runtime_error(error);
                                }
                            }
                            cout << endl;
                        }
                    }
                    else{
                        for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                        {
                            vector<double> fm = cp_FM->postrmf(xi);

                            if( cout.good() == true )
                            {
                                cout << setw(13) << setprecision(12) << fixed << left;

                                cout << "Fisher matrix:" << endl;
                                for(unsigned int i = 0; i<fm.size(); i++)
                                {

                                    //#mb Printout of Fisher matrix
                                    cout << scientific << fm[i];
                                    if( (i+1)%dim==0 )
                                    {
                                        cout << "\n";
                                    }
                                    else
                                    {
                                        cout << "\t";
                                    }
                                }

                                //fs.close();
                            }
                            else{
                                string error = "Can't open stream!\n";
                                throw runtime_error(error);
                            }
                        }   cout << endl;


                        //#mb Fisher matrix written to the grid.bin file
                        vector<double> fish_mat = cp_FM->postrmf(InitialTimeMin);
                        gridbin.write((char*)&fish_mat[0], fish_mat.size()*sizeof(double));

                        //#mb Mn matrix

                        cout << "Grid matrix:" << endl;
                        for(unsigned int i = 0; i<Mn.size(); i++) {
                            cout << scientific << Mn[i];
                            if( (i+1)%dim==0 )
                                cout << "\n";
                            else
                                cout << "\t";
                        }

                        gridbin.write((char*)&Mn[0], Mn.size()*sizeof(double));
                        gridbin.close();
                    }


                    delete cp_FM;

                } // End of 'if(DataLength==0)'
                else
                {

                    if (SaveDensityOfCovering == "True")
                    {
                        //ofstream fs;
                        //fs.open( PathSave + name, fstream::out | fstream::trunc);
                        if( cout.good() == true )
                        {
                            cout << "Covariance, Density, Algorithm\n";
                            for(double c0=CovarianceMin; c0<=CovarianceMax; c0+=CovarianceStep)
                            {
                                int ch_method = 0;
                                if(ChooseMethod == "s1")
                                    ch_method = 1;

                                if(ChooseMethod == "s2")
                                {
                                    switch(Spindown)
                                    {
                                        case 1:
                                            ch_method = 2;
                                            break;
                                        case 2:
                                            cerr << "Algorithm 's2' is implemented only for spindown = 1.\n";
                                            cerr << "'s1' will be applied instead." ;
                                            ch_method = 1;
                                            break;
                                    }
                                }

                                if(ChooseMethod == "Automatic")
                                    ch_method = method(c0, nfft, DataLength, Spindown);

                                if(ch_method == 0)
                                {
                                    string error = "Bad option value: " + ChooseMethod + ".\n"
                                    "Available options: s1, s2, a; 's2' only for spindown = 1.\n";
                                    throw domain_error(error);
                                }

                                int prec;
                                double density;
                                vector<double> sgrid;
                                if(ch_method == 1)
                                {
                                    DensityS1 ds1 = DensityS1(Spindown);
                                    density = ds1.density(c0, nfft, DataLength);
                                    prec = 6;
                                }

                                if(ch_method == 2)
                                {
                                    DensityS2 ds2 = DensityS2();
                                    density = ds2.density_approx(c0, nfft, DataLength);
                                    prec = 3;
                                }

                                if(d_chop)
                                    num::chop(sgrid);

                                //fs.fill('0');
                                cout.width(4);
                                cout << left << c0 << "\t" << setw(13) << setprecision(prec)
                                << density << " s" << ch_method << "\n";
                            }

                            cout << endl;
                            //fs.close();
                        }
                        else{
                            string error = "Can't open stream!\n";
                            throw runtime_error(error);
                        }
                    }

                }

            }
        }
        else{  /// if flag in command line NOT provided
            ReadConfig rc = ReadConfig();
            vector<OptionValue> vov = rc.get_values();

            CovarianceMin=vov[0].m_d;
            CovarianceMax=vov[1].m_d + num::epsilon();
            CovarianceStep=vov[2].m_d;
            InitialTimeMin=vov[3].m_d;
            InitialTimeMax=vov[4].m_d + num::epsilon();
            InitialTimeStep=vov[5].m_d;
            NfftExpLength=static_cast<unsigned int>(vov[6].m_i);
            DataLength=static_cast<unsigned int>(vov[7].m_i);
            NAlpha=vov[8].m_i;
            NRadius=vov[9].m_i;
            Spindown=static_cast<unsigned int>(vov[10].m_i);
            PathSSB=vov[11].m_s;
            PathDet=vov[12].m_s;
            PathDetSSB=vov[13].m_s;
            PathSave=vov[14].m_s;
            FilePatternGrid=vov[15].m_s;
            FilePatternFM=vov[16].m_s;
            FilePatternDoC=vov[17].m_s;
            SaveGrid=vov[18].m_s;
            SaveFisherMatrix=vov[19].m_s;
            SaveDensityOfCovering=vov[20].m_s;
            DataChop=vov[21].m_s;
            ChooseMethod=vov[22].m_s;
            DataS2Use=vov[23].m_s;
            DataS2SaveBest=vov[24].m_s;
            Convert=vov[25].m_s;
            Quiet=vov[26].m_s;
            SayHello=vov[27].m_s;

            sizeCo = vov[2].m_s.size(); /// Remove sizeC0 -> .m_s is a string !
            sizeIt = vov[5].m_s.size()+1;
            sizeNEL = to_string(NfftExpLength).size();
            sizeDLn = to_string(DataLength).size();
            snel = to_string(NfftExpLength).substr(0,sizeNEL);
            sdln = to_string(DataLength).substr(0,sizeDLn);
            dim = Spindown +3;

            bool s2use = false;
            if (DataS2Use == "True") s2use = true;
            bool s2save = false;
            if (DataS2SaveBest == "True") s2save = true;
            bool d_chop = false;
            if (DataChop == "True") d_chop = true;

            bool conv = false;
            if (Convert == "True") conv = true;
            bool hello = false;
            if (SayHello == "True") hello = true;

            unsigned int nfft = pow(2, NfftExpLength);
            if(hello)
            {
                //char oc = 64;
                manual.print_version();
                //system("pause");
            }

            if(DataLength==0)    // DataLength == 0 means DataLength should be equal to ephemeris length
            {
                vector<string> paths_to_ephemeris;
                paths_to_ephemeris.push_back(PathSSB);
                paths_to_ephemeris.push_back(PathDet);
                paths_to_ephemeris.push_back(PathDetSSB);

                ReadEphemeris RE = ReadEphemeris(paths_to_ephemeris);
                FisherRM const *const cp_FM = new FisherRM(RE.get_ephemeris1(), RE.get_ephemeris2(), Spindown);

                GridS1 gs1 = GridS1(cp_FM, Spindown);
                GridS2 gs2 = GridS2(cp_FM, NAlpha, NRadius);

                unsigned int data_length = cp_FM->get_ephemeris_length();
                size_t dL_Eph_size = to_string(data_length).size();
                string sd_Eph= to_string(data_length).substr(0,dL_Eph_size);
                string named = FilePatternDoC;
                stringmanip::sreplace(named, "%N", snel);
                stringmanip::sreplace(named, "%D", sd_Eph);
                ofstream fsd;

                if (SaveDensityOfCovering == "True")
                {
                    fsd.open( PathSave + "/" + named, fstream::out | fstream::trunc);
                    if( fsd.good() == true )
                    {
                        fsd << "Covariance, Density of covering\n";
                    }
                    else{
                        string error = "Can't open: " + named;
                        throw runtime_error(error);
                    }
                }

                if (SaveGrid == "True" || SaveDensityOfCovering == "True")
                for(double c0=CovarianceMin; c0<=CovarianceMax; c0+=CovarianceStep)
                {
                    //cout << "[c0=" << c0 <<"]\n";
                    int ch_method = 0;
                    if(ChooseMethod == "s1")
                        ch_method = 1;

                    if(ChooseMethod == "s2")
                    {
                        switch(Spindown)
                        {
                            case 1:
                                ch_method = 2;
                                break;
                            case 2:
                                cerr << "Algorithm 's2' is implemented only for spindown = 1.\n";
                                cerr << "'s1' will be applied instead." ;
                                ch_method = 1;
                                break;
                        }
                    }

                    if(ChooseMethod == "Automatic")
                        ch_method = method(c0, nfft, RE.get_length(), Spindown);

                    if(ch_method == 0)
                    {
                        string error = "Bad option value: " + ChooseMethod + ".\n"
                        "Available options: s1, s2, a; 's2' only for spindown = 1.\n";
                        throw domain_error(error);
                    }

                    double density;
                    vector<double> sgrid;
                    if(ch_method == 1)
                    {
                        sgrid = gs1.grid_prim(c0, nfft, data_length);
                        density = gs1.DensityS1::density(sgrid);
                    }

                    if(ch_method == 2)
                    {
                        sgrid = gs2.grid_prim(c0, nfft, s2use, s2save);
                        density = gs2.density(sgrid);
                    }

                    if (SaveDensityOfCovering == "True")
                    {
                        fsd.fill('0');
                        fsd.width(4);
                        fsd << left << c0 << " " << setw(13) << setprecision(12)
                        << density << "\n";
                    }

                        if (SaveGrid == "True")
                        {
                            string nameg = FilePatternGrid;
                            string sc0 = to_string(c0).substr(0,sizeCo); ///! sizeCo
                            stringmanip::sreplace(nameg, "%C", sc0);
                            stringmanip::sreplace(nameg, "%N", snel);

                            for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                            {
                                vector<double> sgrid2;
                                if( conv )
                                {
                                    if( ch_method == 1)
                                       sgrid2 = gs1.convert(c0, xi, sgrid);

                                    if( ch_method == 2)
                                        sgrid2 = gs2.convert(c0, xi, sgrid);
                                }
                                else{
                                    sgrid2 = sgrid;
                                }

                                if(d_chop)
                                    num::chop(sgrid2);

                                string sxi;
                                if(xi<0)
                                    sxi = to_string(xi).substr(0,sizeIt);
                                else
                                    sxi = ("+"+to_string(xi)).substr(0,sizeIt);

                                string name1 = nameg;
                                stringmanip::sreplace(name1, "%I", sxi);

                                fstream fs;
                                fs.open( PathSave + "/" + name1, fstream::out | fstream::trunc);
                                std::cerr << "Output: " << PathSave << std::endl;
                                if( fs.good() == true )
                                {   //fixed << left
                                    fs << setw(15) << setprecision(12) << left;
                                    fs << density << "\n";
                                    for(unsigned int i = 0; i<sgrid2.size(); i++)
                                    {
                                        fs << setw(15) << setprecision(12) << left << sgrid2[i];
                                        if( (i+1)%static_cast<int>(dim)==0 )
                                        {
                                            fs << "\n";
                                        }
                                        else
                                        {
                                            fs << "\t";
                                        }
                                    }

                                    //fs.unget();
                                    fs.close();
                                }
                                else
                                {
                                    string error = "Can't open: " + name1;
                                    throw runtime_error(error);
                                }
                                //cout << name1 << endl;
                            }
                        } // End of 'if(SaveGrid == "True")'
                } // End of 'for(double c0=CovarianceMin ...'

                if( fsd.good() == true )
                    fsd.close();

                if (SaveFisherMatrix == "True")
                for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                {
                    string name = FilePatternFM;
                    string sxi;
                    if(xi<0)
                        sxi = to_string(xi).substr(0,sizeIt);
                    else
                    {
                        if(xi<num::epsilon() && xi>-num::epsilon())
                        {
                            sxi = "0.0";
                        }
                        else
                        {
                            sxi = ("+"+to_string(xi)).substr(0,sizeIt);
                        }

                    }


                    stringmanip::sreplace(name, "%I", sxi);
                    vector<double> fm = cp_FM->postrmf(xi);

                    fstream fs;
                    fs.open( PathSave + "/" + name, fstream::out | fstream::trunc);
                    if( fs.good() == true )
                    {
                        fs << setw(16) << setprecision(15) << fixed << left;
                        for(unsigned int i = 0; i<fm.size(); i++)
                        {
                            fs << fm[i];
                            if( (i+1)%static_cast<int>(dim)==0 ){
                                fs << "\n";
                            }
                            else{
                                fs << "\t";
                            }
                        }

                        fs.close();
                    }
                    else{
                        string error = "Can't open: " + PathSave + "/" + name;
                        throw runtime_error(error);
                    }
                }

                delete cp_FM;

            } // End of 'if(DataLength==0)'
            else
            {

                if (SaveDensityOfCovering == "True")
                {
                    string name = FilePatternDoC;
                    stringmanip::sreplace(name, "%N", snel);
                    stringmanip::sreplace(name, "%D", sdln);

                    ofstream fs;
                    fs.open( PathSave + "/" + name, fstream::out | fstream::trunc);
                    if( fs.good() == true )
                    {
                        fs << "Covariance, Density, Algorithm\n";
                        for(double c0=CovarianceMin; c0<=CovarianceMax; c0+=CovarianceStep)
                        {
                            int ch_method = 0;
                            if(ChooseMethod == "s1")
                                ch_method = 1;

                            if(ChooseMethod == "s2")
                            {
                                switch(Spindown)
                                {
                                    case 1:
                                        ch_method = 2;
                                        break;
                                    case 2:
                                        cerr << "Algorithm 's2' is implemented only for spindown = 1.\n";
                                        cerr << "'s1' will be applied instead." ;
                                        ch_method = 1;
                                        break;
                                }
                            }

                            if(ChooseMethod == "Automatic")
                                ch_method = method(c0, nfft, DataLength);

                            if(ch_method == 0)
                            {
                                string error = "Bad option value: " + ChooseMethod + ".\n"
                                "Available options: s1, s2, a; 's2' only for spindown = 1.\n";
                                throw domain_error(error);
                            }


                            int prec;
                            double density;
                            vector<double> sgrid;
                            if(ch_method == 1)
                            {
                                DensityS1 ds1 = DensityS1(Spindown);
                                density = ds1.density(c0, nfft, DataLength);
                                prec = 6;
                            }

                            if(ch_method == 2)
                            {
                                DensityS2 ds2 = DensityS2();
                                density = ds2.density_approx(c0, nfft, DataLength);
                                prec = 3;
                            }

                            if(d_chop)
                                num::chop(sgrid);

                            //fs.fill('0');
                            fs.width(4);
                            fs << left << c0 << "\t" << setw(13) << setprecision(prec)
                            << density << " s" << ch_method << "\n";
                        }

                        fs.close();
                    }
                    else{
                        string error = "Can't open: " + name;
                        throw runtime_error(error);
                    }
                }

            }
        }


    }   // End of 'try'
    catch (const runtime_error& e){ cout << e.what() << endl; }
    catch (const domain_error& e) { cout << e.what() << endl; }
    catch (...) { cout << "Error unknown type"; }

    return 0;
}

int method(double c0, unsigned int nfft, unsigned int data_length, unsigned int spindown)
{
    int s = 0;
    switch(spindown){
        case 1:
            {
                double sq2= sqrt(2);
                double sq2_min = sq2 - num::epsilon();
                double sq2_max = sq2 + num::epsilon();

                double x = num::delta_omega_zero_prim(c0, nfft, data_length);

                if(x<=sqrt(3))
                {
                    if(x<1.64){
                        if( x>sq2_min && x<sq2_max ) //exception grid S1 == grid S2 = grid A4*
                            s=1;
                        else
                            s=2;
                    }
                    else
                        s=1;
                }
                else{
                    if(x<1.81)
                        s=2;
                    else
                        s=1;
                }
            }
            break;

        case 2:
            s=1;
            break;

    }

    return s;
}

bool if_exist(map<string, vector<string> >& map_name, string command)
{
    return map_name.find(command) != map_name.end();
}


