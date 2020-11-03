///////////////////////////////////////
// Name:        main.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     13/10/2015
// Modification:14/06/2020 A.Pisarski
///////////////////////////////////////

#include "ReadConfig.h"
#include "OptionValue.h"
#include "num.h"
#include "stringmanip.h"
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
/// All-Sky searches:
#include "DensityS1.h"
#include "DensityS2.h"
#include "GridS1.h"
#include "GridS2.h"
#include "FisherRM.h"
#include "ReadData.h"
#include "ReadEphemeris.h"
/// Directed searches: (from line 1814)
#include "DensityS1DS.h"
#include "GridS1DS.h"
#include "GridS2DS.h"
#include "FisherRMDS.h"
/// For creating directory
#include <sys/types.h>
#include <sys/stat.h>


using namespace std;

int method(double c0, unsigned int nfft, unsigned int data_length, unsigned int spindown=1);
int methodDS(double c0, unsigned int nfft, unsigned int data_length, unsigned int spindown=1);
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
        options_available.push_back("-i");  // initial time (for equations where ti: <-1/2 To, 1/2 To>. To - observation time).
        options_available.push_back("-I");  // initial time (for equations where ti: <0, To>. To - observation time).
        options_available.push_back("-n");  // number of Fourier bins (nfft length)
        options_available.push_back("-nd"); // number of data
        options_available.push_back("-na"); // number of loops (incrementation deep) in root finding algorithm
        options_available.push_back("-nr"); // number of loops in covering radius (deep hole) finding algorithm
        options_available.push_back("-p");  // print in original (author) mode: g - grid; f - fisher matrix; d - density of covering; a - all;
                                            // (without any of this option program will work,
                                            // with Michal Bejger's sets (result saved to the 'grid.bin' file);
        options_available.push_back("-ch"); // chop (set to zero if number is less than num::epsilon)
        options_available.push_back("-cv"); // convert grid from hyper-sphere to hyper-ellipsoid space
        options_available.push_back("-v");  // program version
        options_available.push_back("-s");  // spindown
        options_available.push_back("-aa"); // show information about author(s)
        options_available.push_back("-t");  // type of searches: a - all sky, d - directed.
        options_available.push_back("-u");  // variance estimator can be unbiased, or biased (default set: unbiased == true)

        vector<string> options_available_long;              // order is important (mast be the same as without _long)
        options_available_long.push_back("--help");         // help
        options_available_long.push_back("--helptxt");      // print help to text file.
        options_available_long.push_back("--covariance");   // covariance   == (minimal match)^2
        options_available_long.push_back("--match");        // minimal match
        options_available_long.push_back("--directory");    // directory
        options_available_long.push_back("--algorithm");    // algorithm (to choose)
        options_available_long.push_back("--initial");      // initial time (for equations where ti: <-1/2 To, 1/2 To>. To - observation time).
        options_available_long.push_back("--Initial");      // initial time (for equations where ti: <0, To>. To - observation time).
        options_available_long.push_back("--nfft");         // number of Fourier bins (nfft length)
        options_available_long.push_back("--ndata");        // number of data
        options_available_long.push_back("--nalpha");       // number of loops (incrementation deep) in root finding algorithm
        options_available_long.push_back("--nradius");      // number of loops in covering radius (deep hole) finding algorithm
        options_available_long.push_back("--print");        // print in original (author) mode: g - grid; f - fisher matrix; d - density of covering; a - all;
                                                            // (without any of this option program will work,
                                                            // with Michal Bejger's sets (result saved to the 'grid.bin' file);
        options_available_long.push_back("--chop");         // chop (set to zero if number is less than num::epsilon)
        options_available_long.push_back("--convert");      // convert grid from hyper-sphere to hyper-ellipsoid space
        options_available_long.push_back("--version");      // program version
        options_available_long.push_back("--spindown");     // spindown
        options_available_long.push_back("--author");       // show information about author(s)
        options_available_long.push_back("--type");         // type of searches: a - all sky, d - directed.
        options_available_long.push_back("--unbiased");     // variance estimator can be unbiased, or biased (default set: unbiased == true)

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
        /// Remove unnecessary (now) command (in case when flags "-m"
        /// and "-c" have been provided together, only flag "-c" will be accepted):
        if( if_exist(options_found, "-m") ) options_found.erase("-m");


        /// Convert initial time from notation "-I" to "-i"
        if( !if_exist(options_found, "-i") &&  if_exist(options_found, "-I") )
        {
            map<string, vector<string> >::const_iterator it = options_found.find("-I");
            vector<string>::const_iterator line_it = it->second.begin();

            if(*line_it!="")
            {
                while (line_it != it->second.end())
                {
                    double ti = stod(*line_it) - 0.5;
                    string t_init = to_string(ti);
                    options_found["-i"].push_back(t_init);
                    ++line_it;
                }
            }
            else{
                options_found["-i"].push_back("");
            }

            options_found.erase("-I");
        }
        /// Remove unnecessary (now) command (in case when flags "-I"
        /// and "-i" have been provided together, only flag "-i" will be accepted):
        if( if_exist(options_found, "-I") ) options_found.erase("-I");

        /// Create defaults flags (user don't need to provide them).
        /// Adds flags with default or temporary default flag option "",
        for(unsigned int i=0; i<options_available.size(); ++i)
        {
            if( !if_exist(options_found, options_available[i]) )
            {   /// ... except flags which don't need or don't have flag option:
                if(options_available[i]!="-h" && options_available[i]!="-ht"
                   && options_available[i]!="-v" && options_available[i]!="-aa"
                   && options_available[i]!="-c" && options_available[i]!="-s"
                   && options_available[i]!="-t" && options_available[i]!="-u") // && options_available[i]!="-m"
                        options_found[options_available[i]].push_back("");
            }

        }

        Manual manual("\t * Build 0.3.07 (alpha).                   *");

        double CovarianceMin, CovarianceMax, CovarianceStep=0.01;
        double InitialTimeMin=0.5, InitialTimeMax, InitialTimeStep;
        unsigned int NfftExpLength=20, DataLength=0, Spindown=1, dim = 4, DataLengthDS = 344656, mode = 1;
        /// DataLengthDS - Default data length for directed searches.
        /// mode = {1, 2} only for grid 2 (S2) for directed searches and first spindown.
        /// mode = 2 need to be used with caution.
        int NAlpha=35, NRadius=20;
        string PathSSB, PathDet, PathDetSSB;
        string Path, SegmentNo, DetectorH="", DetectorL="", DetectorV="", Band, DataFile;
        unsigned int PathK; /// Number of strings in pathname.
        vector<string> paths;
        string PathSave;
        string FilePatternGrid, FilePatternFM, FilePatternDoC;
        string SaveGrid="False", SaveFisherMatrix="False", SaveDensityOfCovering="False";
        string DataChop="False";
        string ChooseMethod="Automatic"; //s1, s2
        string TypeSearch="All-Sky";    /// v0.3.02
        string DataS2Use, DataS2SaveBest;
        string Convert="True";
        string ObserveTimeSymmetric = "True";   /// v.0.3.02 (False == initial time - ti: <0, To>. True == ti: <-1/2 To, 1/2 To>).
        string Quiet;					        /// ObserveTimeSymmetric: only for options read from file 'gg.ini'.
        string SayHello;
        bool unbiased = true;
        size_t sizeCo, sizeIt, sizeNEL, sizeDLn;
        string snel, sdln;
        InitialTimeMax=InitialTimeMin + num::epsilon();
        InitialTimeStep=1.0;

        string gridout;
        vector<double> M, Mn;

        if( commands.size() != 0) /// if flags in command line are provided
        {
            ///Help
            if( if_exist(options_found, "-h") )
                manual.print_help();

            ///Print a help to text file
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
                if( if_exist(options_found, "-t") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-t");
                    vector<string>::const_iterator line_it = it->second.begin();

                    if(*line_it=="" || *line_it=="a")
                        TypeSearch="All-Sky";
                    if(*line_it=="d")
                        TypeSearch="Directed";

                    if(*line_it!="" && *line_it!="a" && *line_it!="d"){
                        string error = "Allowed only types of searches: a \"All-Sky\" and d \"Directed\".\n";
                        throw domain_error(error);
                    }
                }

                /// Directory: Path, PathSave. SegmentNo, Detector#, Band, DataFile
                string lastDataName = "";
                if( if_exist(options_found, "-d") && TypeSearch=="All-Sky" )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-d");
                    vector<string>::const_iterator line_it = it->second.begin();
                    vector<string> temp;
                    PathK=0;

                    if(*line_it=="")
                    {
                        PathSave = "";
                    }
                    else{
                        while (line_it != it->second.end()){
                            temp.push_back(*line_it);
                            ++line_it; PathK++;
                        }

                        if(PathK == 1){
                            PathSave = temp[0];
                        }
                        if(PathK >=3 && PathK <= 7){
                            Path=temp[0];
                            SegmentNo=temp[1];

                            if (std::find(temp.begin()+2, temp.end(), "H1") != temp.end()){
                                DetectorH = "H1";
                            }

                            if (std::find(temp.begin()+2, temp.end(), "L1") != temp.end()){
                                DetectorL = "L1";
                            }

                            if (std::find(temp.begin()+2, temp.end(), "V1") != temp.end()){
                                DetectorV = "V1";
                            }

                            if(PathK >=5 && PathK <= 7){
                                string last = temp[temp.size()-1];
                                string dataFileName2find = "xdats";
                                if(last.find(dataFileName2find)!=std::string::npos){
                                    Band = temp[temp.size()-2];
                                    if(last.find(".bin")!=std::string::npos){
                                        size_t si = last.find("_");
                                        if(si>=6){
                                            lastDataName = last.substr(si-1, 1);// to add last letter xdats"C" to output
                                            // for test:
                                            //cout << lastDataName << endl;
                                        }

                                        DataFile = last;
                                    }
                                    else{
                                        size_t si = last.size();
                                        if(si>=6){
                                            lastDataName = last.substr(si-1, 1);
                                            // for test:
                                            //cout << lastDataName << endl;
                                        }

                                        DataFile = last + "_" + SegmentNo + "_" + Band + ".bin";
                                    }
                                }
                                else{
                                    Band = last;
                                    DataFile = dataFileName2find + "_" + SegmentNo + "_" + Band + ".bin";
                                }
                            }
                        }
                        if(PathK != 1 && (PathK <3 || PathK > 7)){
                            string error1 = "(1): -d path_to_ephemeris/.\n";
                            string error2 = "(2): -d path/ SegmentNo detectorX band.\n";
                            string error3 = "(3): -d path/ SegmentNo detectorX detectorY band.\n";
                            string error4 =" (4): -d path/ SegmentNo detectorX detectorY detectorZ band.\n detector# == H1 or L1 or V1.\n";
                            string error = error1 + error2 + error3 + error4;
                            throw domain_error(error);
                        }
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
                    ///cout << "InitialTimeMin= " << InitialTimeMin << endl;
                    ///cout << "InitialTimeMax= " << InitialTimeMax << endl;
                    ///cout << "InitialTimeStep= " << InitialTimeStep << endl;
                }



                /// Need to set a spindown and dimension for type search: "Directed"
                /// but spindown is not provided

                if( TypeSearch=="Directed" && !if_exist(options_found, "-s") ) {
                    Spindown=1;
                    dim = Spindown + 1;
                }

                ///Spindown
                if( if_exist(options_found, "-s") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-s");
                    vector<string>::const_iterator line_it = it->second.begin();
                    if( TypeSearch=="All-Sky") {
                        if(*line_it=="")
                            Spindown=1;
                        else{
                            unsigned int tt = stoi(*line_it);
                            if(tt == 1 || tt == 2)
                                Spindown = tt;
                            else{
                                string error = "Allowed spindown values for \"All-Sky\" searches: {1, 2}.\n";
                                throw domain_error(error);
                            }
                        }
                        dim = Spindown + 3;
                    }

                    if( TypeSearch=="Directed") {
                        if(*line_it=="")
                            Spindown=1;
                        else{
                            unsigned int tt = stoi(*line_it);
                            if(tt == 1 || tt == 2 || tt == 3)
                                Spindown = tt;
                            else{
                                string error = "Allowed spindown values for \"Directed\" searches: {1, 2, 3}.\n";
                                throw domain_error(error);
                            }
                        }
                        dim = Spindown + 1;
                    }
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
                            string error = "Data length can not be set on: "+*line_it+"\n";
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
                bool no_match = true;
                if( if_exist(options_found, "-p") )
                {
                    //std::cout << "test --- " << std::endl;
                    map<string, vector<string> >::const_iterator it = options_found.find("-p");
                    vector<string>::const_iterator line_it = it->second.begin();

                    size_t ts_end = (*line_it).size();
                    if(ts_end > 4)
                        ts_end = 4;

                    if(*line_it==""){
                        SaveGrid=="False"; SaveFisherMatrix=="False"; SaveDensityOfCovering=="False";
                        no_match = false;/// !!!
                    }
                    else{
                        string temp = string((*line_it).begin(),(*line_it).begin()+ts_end); //[,)

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
                        "f (Fisher reduced matrix), g (grid), a (all: density, Fisher matrix, grid).\n";
                        //"t - print on terminal.\n";
                        throw domain_error(error);
                    }

                    /* // for testing
                    cout << "test *1*" << no_match << endl;
                    cout << "SaveDensityOfCovering = " << SaveDensityOfCovering << endl;
                    cout << "SaveFisherMatrix = " << SaveFisherMatrix << endl;
                    cout << "SaveGrid = " << SaveGrid << endl;
                    */
                }
                else{
                    SaveGrid=="False"; SaveFisherMatrix=="False"; SaveDensityOfCovering=="False";
                    no_match = false;
                }
                //cout<< no_match << endl;


                ///Unbiased
                if( if_exist(options_found, "-u") )
                {
                    map<string, vector<string> >::const_iterator it = options_found.find("-u");
                    vector<string>::const_iterator line_it = it->second.begin();

                    if(*line_it=="" || *line_it=="t")
                        unbiased=true;
                    else{
                        if(*line_it=="f")
                            unbiased= false;
                        else{
                            string error = "Arguments for -u (--unbiased) can be only: t (true) or f (false).\n";
                            throw domain_error(error);
                        }
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
        */      bool s2use = false;
                if (DataS2Use == "True") s2use = true;
                bool s2save = false;
                if (DataS2SaveBest == "True") s2save = true;
                bool d_chop = false;
                if (DataChop == "True") d_chop = true;

                bool conv = false;
                if (Convert == "True") conv = true;

                unsigned int nfft = pow(2, NfftExpLength);

                // for test
                ///cerr << "ChooseMethod: " << ChooseMethod << endl;

                //cout << "DataLength= " << DataLength << endl;
                if(TypeSearch == "All-Sky" && Spindown == 1){
                    /// 01 All-Sky. Spindown = 1. Options read from command line.
                    if(DataLength==0)    // DataLength == 0 means DataLength should be equal to ephemeris length
                    {
                        map<string, string> gridouts;
                        vector<string> detectors;
                        unsigned int detectors_size;
                        if(DetectorH == "H1"){
                            detectors.push_back("H1");
                        }
                        if(DetectorL == "L1"){
                            detectors.push_back("L1");
                        }
                        if(DetectorV == "V1"){
                            detectors.push_back("V1");
                        }

                        detectors_size = detectors.size();

                        switch(detectors_size){
                            /*case 0:
                                gridouts["default"]=PathSave;
                                break;*/
                            case 1:
                                gridouts[detectors[0]]=Path + SegmentNo + "/" + detectors[0]+ "/";
                                break;
                            case 2:
                                gridouts[detectors[0]]=Path + SegmentNo + "/" + detectors[0]+ "/";
                                gridouts[detectors[1]]=Path + SegmentNo + "/" + detectors[1]+ "/";
                                break;
                            case 3:
                                gridouts[detectors[0]]=Path + SegmentNo + "/" + detectors[0]+ "/";
                                gridouts[detectors[1]]=Path + SegmentNo + "/" + detectors[1]+ "/";
                                gridouts[detectors[2]]=Path + SegmentNo + "/" + detectors[2]+ "/";
                                break;
                        }

                        FisherRM const * cp_FM;

                        vector<double> sigma4deta;
                        vector<string> paths_to_ephemeris;
                        string full_path_to_data;
                        vector<vector<double> > detectors_ephemeris;

                        if(detectors_size >= 2){
                            ///#mb binary file grid bin
                            gridout = Path + SegmentNo + "/grids/";
                            struct stat info;
                            int check = 0;
                            if( stat( gridout.c_str(), &info ) != 0){
                                check = mkdir( gridout.c_str(), 0777);
                            }
                            if(check){
                                std::string error = "Unable to create directory: " + gridout;
                            }
                            gridout+="grid_" + SegmentNo + "_" + Band + "_";

                            //string path_to_data = Path + SegmentNo + "/";
                            for (map<string, string>::const_iterator go = gridouts.begin(); go != gridouts.end(); ++go )
                            {
                                gridout+=go->first;
                                //full_path_to_data = path_to_data + go->first + "/" + DataFile;
                                full_path_to_data = go->second + DataFile;

                                PathSSB = go->second + "rSSB.bin";
                                PathDet = go->second + "rDet.bin";
                                PathDetSSB = go->second + "DetSSB.bin";

                                //cout << PathSSB << endl;

                                paths_to_ephemeris.push_back(PathSSB);
                                paths_to_ephemeris.push_back(PathDet);
                                paths_to_ephemeris.push_back(PathDetSSB);

                                ReadEphemeris RE = ReadEphemeris(paths_to_ephemeris);
                                //ephemeris_length = RE.get_length();
                                detectors_ephemeris.push_back(RE.get_ephemeris1());
                                detectors_ephemeris.push_back(RE.get_ephemeris2());

                                ReadData RD = ReadData(full_path_to_data);
                                double average = num::avg(RD.get_data());
                                double variance = num::var(RD.get_data(), average, unbiased);
                                cout.precision(15);
                                cout << "avg[" << go->first << "] = " << average;
                                cout << ", \u03C3" << "[" << go->first << "] = " << sqrt(variance);
                                cout << ", \u03C3" << "\u00B2" << "[" << go->first << "] = " << variance;
                                cout.precision(6); // default set
                                if(unbiased){
                                    cout << " (unbiased = true).\n";
                                }
                                else{
                                    cout << " (unbiased = false).\n";
                                }
                                sigma4deta.push_back(variance);
                            }
                            cp_FM = new FisherRM(detectors_ephemeris, sigma4deta, Spindown);

                            if(lastDataName!=""){
                                gridout+=lastDataName;
                            }
                            ///#mb binary file grid bin
                            gridout+=".bin";
                        }
                        else{
                            ///#mb binary file grid bin
                            if(detectors_size == 1){
                                gridout = Path + SegmentNo + "/" + detectors[0] + "/";
                                PathSSB = gridout + "rSSB.bin";
                                PathDet = gridout + "rDet.bin";
                                PathDetSSB = gridout + "DetSSB.bin";
                                gridout+="grid.bin";
                            }
                            else{
                                gridout = PathSave + "grid.bin";

                                PathSSB = PathSave + "rSSB.bin";
                                PathDet = PathSave + "rDet.bin";
                                PathDetSSB = PathSave + "DetSSB.bin";
                            }

                            paths_to_ephemeris.push_back(PathSSB);
                            paths_to_ephemeris.push_back(PathDet);
                            paths_to_ephemeris.push_back(PathDetSSB);

                            ReadEphemeris RE = ReadEphemeris(paths_to_ephemeris);
                            //ephemeris_length = RE.get_length();
                            cp_FM = new FisherRM(RE.get_ephemeris1(), RE.get_ephemeris2(), Spindown);
                        }

                        std::ofstream gridbin(gridout, std::ios::binary);
                        if(gridbin.fail())
                        {
                            std::string error="Can not find: " + gridout;
                            throw std::runtime_error(error);
                        }
                        if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False")
                            cout << "grid.bin will be saved in " << gridout << endl;

                        GridS1 gs1 = GridS1(cp_FM, Spindown);
                        GridS2 gs2 = GridS2(cp_FM, NAlpha, NRadius); /// Only for first spindown!

                        unsigned int data_length = cp_FM->get_ephemeris_length();
                        //ofstream fsd;

                        if (SaveDensityOfCovering == "True")
                        {
                            if( cout.good() == true ){
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
                                ch_method = 2;

                            if(ChooseMethod == "Automatic")
                                ch_method = method(c0, nfft, data_length, Spindown);

                            if(ch_method == 0)
                            {
                                string error = "Bad option value: " + ChooseMethod + ".\n"
                                "Available options: s1, s2, a.\n";
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

                            if (SaveDensityOfCovering == "True")
                            {
                                cout.fill('0');
                                cout.width(4);
                                cout << left << c0 << ", " << setw(13) << setprecision(12)
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
                        } // End of 'for(double c0=CovarianceMin ...'

                            /// Default Michal Bejger's options
                            if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False"){
                                for(double c0=CovarianceMin; c0<=CovarianceMax; c0+=CovarianceStep)
                                {
                                    int ch_method = 0;
                                    if(ChooseMethod == "s1")
                                        ch_method = 1;

                                    if(ChooseMethod == "s2")
                                        ch_method = 2;

                                    if(ChooseMethod == "Automatic")
                                        ch_method = method(c0, nfft, data_length, Spindown);

                                    if(ch_method == 0)
                                    {
                                        string error = "Bad option value: " + ChooseMethod + ".\n"
                                        "Available options: s1, s2, a.\n";
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
                                    //#mb
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
                                        else{
                                            string error = "Can't open stream!\n";
                                            throw runtime_error(error);
                                        }
                                        //cout << name1 << endl;
                                    }
                                    cout << endl;
                                } // End of 'for(double c0=CovarianceMin ...'
                            }

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

                            /// Default Michal Bejger's options
                            if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False"){
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
                                }
                                cout << endl;

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
                    else{
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
                                        ch_method = 2;

                                    if(ChooseMethod == "Automatic")
                                        ch_method = method(c0, nfft, DataLength, Spindown);

                                    if(ch_method == 0)
                                    {
                                        string error = "Bad option value: " + ChooseMethod + ".\n"
                                        "Available options: s1, s2, a.\n";
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
                if(TypeSearch == "All-Sky" && Spindown == 2){
                    /// 02 All-Sky. Spindown = 2. Options read from command line.
                    if(DataLength==0)    // DataLength == 0 means DataLength should be equal to ephemeris length
                    {
                        map<string, string> gridouts;
                        vector<string> detectors;
                        unsigned int detectors_size;
                        if(DetectorH == "H1"){
                            detectors.push_back("H1");
                        }
                        if(DetectorL == "L1"){
                            detectors.push_back("L1");
                        }
                        if(DetectorV == "V1"){
                            detectors.push_back("V1");
                        }

                        detectors_size = detectors.size();

                        switch(detectors_size){
                            /*case 0:
                                gridouts["default"]=PathSave;
                                break;*/
                            case 1:
                                gridouts[detectors[0]]=Path + SegmentNo + "/" + detectors[0]+ "/";
                                break;
                            case 2:
                                gridouts[detectors[0]]=Path + SegmentNo + "/" + detectors[0]+ "/";
                                gridouts[detectors[1]]=Path + SegmentNo + "/" + detectors[1]+ "/";
                                break;
                            case 3:
                                gridouts[detectors[0]]=Path + SegmentNo + "/" + detectors[0]+ "/";
                                gridouts[detectors[1]]=Path + SegmentNo + "/" + detectors[1]+ "/";
                                gridouts[detectors[2]]=Path + SegmentNo + "/" + detectors[2]+ "/";
                                break;
                        }

                        FisherRM const * cp_FM;

                        vector<double> sigma4deta;
                        vector<string> paths_to_ephemeris;
                        string full_path_to_data;
                        vector<vector<double> > detectors_ephemeris;
                        //unsigned int ephemeris_length = -1;

                        if(detectors_size >= 2){
                            ///#mb binary file grid bin
                            gridout = Path + SegmentNo + "/grids/";
                            struct stat info;
                            int check = 0;
                            if( stat( gridout.c_str(), &info ) != 0){
                                check = mkdir( gridout.c_str(), 0777);
                            }
                            if(check){
                                std::string error = "Unable to create directory: " + gridout;
                            }
                            gridout+="grid_" + SegmentNo + "_" + Band + "_";
                            //gridout = Path + SegmentNo + "/grids/grid_" + SegmentNo + "_" + Band + "_";
                            //string path_to_data = Path + SegmentNo + "/";
                            for (map<string, string>::const_iterator go = gridouts.begin(); go != gridouts.end(); ++go )
                            {
                                gridout+=go->first;
                                //full_path_to_data = path_to_data + go->first + "/" + DataFile;
                                full_path_to_data = go->second + DataFile;

                                PathSSB = go->second + "rSSB.bin";
                                PathDet = go->second + "rDet.bin";
                                PathDetSSB = go->second + "DetSSB.bin";

                                paths_to_ephemeris.push_back(PathSSB);
                                paths_to_ephemeris.push_back(PathDet);
                                paths_to_ephemeris.push_back(PathDetSSB);

                                ReadEphemeris RE = ReadEphemeris(paths_to_ephemeris);
                                //ephemeris_length = RE.get_ephemeris_length();
                                detectors_ephemeris.push_back(RE.get_ephemeris1());
                                detectors_ephemeris.push_back(RE.get_ephemeris2());

                                ReadData RD = ReadData(full_path_to_data);
                                double average = num::avg(RD.get_data());
                                double variance = num::var(RD.get_data(), average, unbiased);
                                cout.precision(15);
                                cout << "avg[" << go->first << "] = " << average;
                                cout << ", \u03C3" << "[" << go->first << "] = " << sqrt(variance);
                                cout << ", \u03C3" << "\u00B2" << "[" << go->first << "] = " << variance;
                                cout.precision(6); // default set
                                if(unbiased){
                                    cout << " (unbiased = true).\n";
                                }
                                else{
                                    cout << " (unbiased = false).\n";
                                }
                                sigma4deta.push_back(variance);
                            }
                            cp_FM = new FisherRM(detectors_ephemeris, sigma4deta, Spindown);

                            ///#mb binary file grid bin
                            if(lastDataName!=""){
                                gridout+=lastDataName;
                            }
                            gridout+=".bin";
                        }
                        else{
                            ///#mb binary file grid bin
                            if(detectors_size == 1){
                                gridout = Path + SegmentNo + "/" + detectors[0] + "/";
                                PathSSB = gridout + "rSSB.bin";
                                PathDet = gridout + "rDet.bin";
                                PathDetSSB = gridout + "DetSSB.bin";
                                gridout+="grid.bin";
                            }
                            else{
                                gridout = PathSave + "grid.bin";

                                PathSSB = PathSave + "rSSB.bin";
                                PathDet = PathSave + "rDet.bin";
                                PathDetSSB = PathSave + "DetSSB.bin";
                            }

                            paths_to_ephemeris.push_back(PathSSB);
                            paths_to_ephemeris.push_back(PathDet);
                            paths_to_ephemeris.push_back(PathDetSSB);

                            ReadEphemeris RE = ReadEphemeris(paths_to_ephemeris);
                            //ephemeris_length = RE.get_ephemeris_length();
                            cp_FM = new FisherRM(RE.get_ephemeris1(), RE.get_ephemeris2(), Spindown);
                        }

                        std::ofstream gridbin(gridout, std::ios::binary);
                        if(gridbin.fail())
                        {
                            std::string error="Can not find: " + gridout;
                            throw std::runtime_error(error);
                        }
                        if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False")
                            cout << "grid.bin will be saved in " << gridout << endl;

                        GridS1 gs1 = GridS1(cp_FM, Spindown);

                        unsigned int data_length = cp_FM->get_ephemeris_length();
                        //ofstream fsd;

                        if (SaveDensityOfCovering == "True")
                        {
                            if( cout.good() == true ){
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

                            if(ChooseMethod == "s2"){
                                cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                                cerr << "'s1' will be applied instead.";
                                ch_method = 1;
                            }

                            if(ChooseMethod == "Automatic")
                                ch_method = 1;

                            if(ch_method == 0)
                            {
                                string error = "Bad option value: " + ChooseMethod + ".\n"
                                "Available options: s1, a (s2 only for spindown = 1).\n";
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

                            if (SaveDensityOfCovering == "True")
                            {
                                cout.fill('0');
                                cout.width(4);
                                cout << left << c0 << ", " << setw(13) << setprecision(12)
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
                        } // End of 'for(double c0=CovarianceMin ...'


                        /// Default Michal Bejger's options
                        if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False"){
                            for(double c0=CovarianceMin; c0<=CovarianceMax; c0+=CovarianceStep)
                            {
                                int ch_method = 0;
                                if(ChooseMethod == "s1")
                                    ch_method = 1;

                                if(ChooseMethod == "s2"){
                                    cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                                    cerr << "'s1' will be applied instead.";
                                    ch_method = 1;
                                }

                                if(ChooseMethod == "Automatic")
                                    ch_method = 1;

                                if(ch_method == 0)
                                {
                                    string error = "Bad option value: " + ChooseMethod + ".\n"
                                    "Available options: s1, a (s2 only for spindown = 1).\n";
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
                                //#mb
                                for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                                {
                                    vector<double> sgrid2;
                                    if( conv )
                                    {
                                        if( ch_method == 1)
                                           sgrid2 = gs1.convert(c0, xi, sgrid);
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
                                    else{
                                        string error = "Can't open stream!\n";
                                        throw runtime_error(error);
                                    }
                                    //cout << name1 << endl;
                                }
                                cout << endl;
                            } // End of 'for(double c0=CovarianceMin ...'
                        } // End of if

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

                        /// Default Michal Bejger's options
                        if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False"){
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
                            }
                            cout << endl;

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
                    else{
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
                                        cerr << "Algorithm 's2' is implemented only for spindown = 1.\n";
                                        cerr << "'s1' will be applied instead." ;
                                        ch_method = 1;
                                    }

                                    if(ChooseMethod == "Automatic")
                                        ch_method = 1;

                                    if(ch_method == 0)
                                    {
                                        string error = "Bad option value: " + ChooseMethod + ".\n"
                                        "Available options: s1, a (s2 only for spindown = 1).\n";
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
                // Code for directed search:
                if(TypeSearch == "Directed" && Spindown == 1){
                    /// 03 Directed. Spindown = 1. Options read from command line.
                    unsigned int data_length;
                    if(DataLength==0){
                        data_length = DataLengthDS; /// Default length of data 'DataLengthDS'. DataLengthDS =
                    }
                    else{
                        data_length = DataLength;
                    }

                    //#mb binary file grid bin
                    gridout = PathSave + "grid.bin";
                    std::ofstream gridbin(gridout, std::ios::binary);
                    if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False")
                        cout << "grid.bin will be saved in " << gridout << endl;

                    FisherRMDS const *const cp_FMDS = new FisherRMDS(Spindown);

                    GridS1DS gs1DS = GridS1DS(cp_FMDS, Spindown);
                    GridS2DS gs2DS = GridS2DS(cp_FMDS);

                    if (SaveDensityOfCovering == "True")
                    {
                        if( cout.good() == true ){
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
                            ch_method = 2;

                        if(ChooseMethod == "Automatic")
                            ch_method = methodDS(c0, nfft, data_length, Spindown);

                        if(ch_method == 0)
                        {
                            string error = "Bad option value: " + ChooseMethod + ".\n"
                            "Available options: s1, s2, a.\n";
                            throw domain_error(error);
                        }
                        ///cout << "test 1" << endl;
                        double densityDS;
                        vector<double> sgridDS;
                        if(ch_method == 1)
                        {
                            sgridDS = gs1DS.grid_prim(c0, nfft, data_length);
                            densityDS = gs1DS.DensityS1DS::density(sgridDS);
                        }
                        ///cout << "test 2" << endl;
                        if(ch_method == 2)
                        {
                            sgridDS = gs2DS.grid_prim(c0, nfft, data_length, mode);
                            densityDS = gs2DS.density(sgridDS);
                        }

                        if (SaveDensityOfCovering == "True")
                        {
                            cout.fill('0');
                            cout.width(4);
                            cout << left << c0 << " " << setw(13) << setprecision(12)
                            << densityDS << "\n";
                        }

                        if (SaveGrid == "True")
                        {
                            for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                            {
                                vector<double> sgridDS2;
                                if( conv )
                                {
                                    if( ch_method == 1)
                                       sgridDS2 = gs1DS.convert(c0, xi, sgridDS);

                                    if( ch_method == 2)
                                        sgridDS2 = gs2DS.convert(c0, xi, sgridDS);
                                }
                                else{
                                    sgridDS2 = sgridDS;
                                }

                                if(d_chop)
                                    num::chop(sgridDS2);

                                //fstream fs;
                                //fs.open( PathSave + "\/" + name1, fstream::out | fstream::trunc);
                                if( cout.good() == true )
                                {   //fixed << left
                                    cout << setw(15) << setprecision(12) << left;
                                    cout << densityDS << "\n";
                                    for(unsigned int i = 0; i<sgridDS2.size(); i++)
                                    {
                                        cout << setw(15) << setprecision(11) << left << sgridDS2[i];
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
                    } // End of 'for(double c0=CovarianceMin ...'


                    /// Default Michal Bejger's options
                    if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False"){
                        for(double c0=CovarianceMin; c0<=CovarianceMax; c0+=CovarianceStep)
                        {
                            int ch_method = 0;
                            if(ChooseMethod == "s1")
                                ch_method = 1;

                            if(ChooseMethod == "s2")
                                ch_method = 2;

                            if(ChooseMethod == "Automatic")
                                ch_method = methodDS(c0, nfft, data_length, Spindown);

                            if(ch_method == 0)
                            {
                                string error = "Bad option value: " + ChooseMethod + ".\n"
                                "Available options: s1, s2, a.\n";
                                throw domain_error(error);
                            }
                            ///cout << "test 1" << endl;
                            double densityDS;
                            vector<double> sgridDS;
                            if(ch_method == 1)
                            {
                                sgridDS = gs1DS.grid_prim(c0, nfft, data_length);
                                densityDS = gs1DS.DensityS1DS::density(sgridDS);
                            }
                            ///cout << "test 2" << endl;
                            if(ch_method == 2)
                            {
                                sgridDS = gs2DS.grid_prim(c0, nfft, data_length, mode);
                                densityDS = gs2DS.density(sgridDS);
                            }
                            //#mb
                            for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                            {
                                vector<double> sgridDS2;
                                if( conv )
                                {
                                    if( ch_method == 1)
                                        sgridDS2 = gs1DS.convert(c0, xi, sgridDS);

                                    if( ch_method == 2)
                                        sgridDS2 = gs2DS.convert(c0, xi, sgridDS);
                                }
                                else{
                                    sgridDS2 = sgridDS;
                                }

                                if(d_chop)
                                    num::chop(sgridDS2);


                                //#mb writing the fftpad to grid.bin
                                int fftpad;
                                //cp_FMDS->get_ephemeris_length()
                                fftpad = (log2(nfft) - ceil(log2(data_length)) + 1);

                                cout << "fftpad: " << fftpad << endl ;

                                gridbin.write((char*) &fftpad, sizeof(int));

                                unsigned int datal = data_length;//cp_FM->get_ephemeris_length();
                                /// 'datal' is this same as 'data_length'
                                M.resize(dim*dim);
                                Mn.resize(dim*dim);
                                for(unsigned int i = 0; i<sgridDS2.size(); i++) {
                                    M[i] = sgridDS2[i]/datal;
                                    Mn[i] = sgridDS2[i];
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
                                    cout << densityDS << "\n";
                                    for(unsigned int i = 0; i<sgridDS2.size(); i++)
                                    {
                                        cout << setw(15) << setprecision(11) << left << sgridDS2[i];
                                        if( (i+1)%dim==0 )
                                            cout << endl;
                                        else
                                            cout << "\t";
                                    }

                                    //fs.unget();
                                    //fs.close();
                                }
                                else{
                                    string error = "Can't open stream!\n";
                                    throw runtime_error(error);
                                }
                                //cout << name1 << endl;
                            }
                            cout << endl;
                        } // End of 'for(double c0=CovarianceMin ...'
                    }

                    if (SaveFisherMatrix == "True")
                    {
                        for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                        {
                            vector<double> fmDS = cp_FMDS->postrmf(xi);

                            if( cout.good() == true )
                            {
                                cout << scientific  << left;//<< setw(16) << setprecision(15) << fixed
                                for(unsigned int i = 0; i<fmDS.size(); i++)
                                {
                                    cout << fmDS[i];
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


                    /// Default Michal Bejger's options
                    if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False"){
                        for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                        {
                            vector<double> fmDS = cp_FMDS->postrmf(xi);

                            if( cout.good() == true )
                            {
                                cout << setw(13) << setprecision(12) << fixed << left;

                                cout << "Fisher matrix:" << endl;
                                for(unsigned int i = 0; i<fmDS.size(); i++)
                                {

                                    //#mb Printout of Fisher matrix
                                    cout << scientific << fmDS[i];
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


                        //#mb Fisher matrix written to the grid.bin file
                        vector<double> fish_mat = cp_FMDS->postrmf(InitialTimeMin);
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

                    delete cp_FMDS;
                }
                if(TypeSearch == "Directed" && Spindown >= 2){
                    // for test:
                    //cout << "dim= " << dim << endl;

                    unsigned int data_length;
                    if(DataLength==0){
                        data_length = DataLengthDS; /// Default length of data 'DataLengthDS'.
                    }
                    else{
                        data_length = DataLength;
                    }

                    //#mb binary file grid bin
                    gridout = PathSave + "grid.bin";
                    std::ofstream gridbin(gridout, std::ios::binary);
                    if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False")
                        cout << "grid.bin will be saved in " << gridout << endl;

                    FisherRMDS const *const cp_FMDS = new FisherRMDS(Spindown);

                    GridS1DS gs1DS = GridS1DS(cp_FMDS, Spindown);

                    if (SaveDensityOfCovering == "True")
                    {
                        if( cout.good() == true ){
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
                                case 2:
                                    cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                                    cerr << "'s1' will be applied instead." ;
                                    ch_method = 1;
                                    break;
                                case 3:
                                    cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                                    cerr << "'s1' will be applied instead." ;
                                    ch_method = 1;
                                    break;
                                default:
                                    cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                                    cerr << "'s1' will be applied instead." ;
                                    ch_method = 1;
                                    break;
                            }
                        }

                        if(ChooseMethod == "Automatic")
                            ch_method = 1;

                        if(ch_method == 0)
                        {
                            string error = "Bad option value: " + ChooseMethod + ".\n"
                            "Available options: s1, a (s2 only for spindown = 1).\n";
                            throw domain_error(error);
                        }
                        ///cout << "test 1" << endl;
                        double densityDS;
                        vector<double> sgridDS;
                        if(ch_method == 1)
                        {
                            sgridDS = gs1DS.grid_prim(c0, nfft, data_length);
                            densityDS = gs1DS.DensityS1DS::density(sgridDS);
                        }
                        ///cout << "test 2" << endl;

                        if (SaveDensityOfCovering == "True")
                        {
                            cout.fill('0');
                            cout.width(4);
                            cout << left << c0 << " " << setw(13) << setprecision(12)
                            << densityDS << "\n";
                        }

                        if (SaveGrid == "True")
                        {
                            for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                            {
                                vector<double> sgridDS2;
                                if( conv )
                                {
                                    if( ch_method == 1)
                                       sgridDS2 = gs1DS.convert(c0, xi, sgridDS);
                                }
                                else{
                                    sgridDS2 = sgridDS;
                                }

                                if(d_chop)
                                    num::chop(sgridDS2);

                                //fstream fs;
                                //fs.open( PathSave + "\/" + name1, fstream::out | fstream::trunc);
                                if( cout.good() == true )
                                {   //fixed << left
                                    cout << setw(15) << setprecision(12) << left;
                                    cout << densityDS << "\n";
                                    for(unsigned int i = 0; i<sgridDS2.size(); i++)
                                    {
                                        cout << setw(15) << setprecision(11) << left << sgridDS2[i];
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
                    } // End of 'for(double c0=CovarianceMin ...'

                    /// Default Michal Bejger's options
                    if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False"){
                        for(double c0=CovarianceMin; c0<=CovarianceMax; c0+=CovarianceStep)
                        {
                            int ch_method = 0;
                            if(ChooseMethod == "s1")
                                ch_method = 1;

                            if(ChooseMethod == "s2")
                            {
                                switch(Spindown)
                                {
                                    case 2:
                                        cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                                        cerr << "'s1' will be applied instead." ;
                                        ch_method = 1;
                                        break;
                                    case 3:
                                        cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                                        cerr << "'s1' will be applied instead." ;
                                        ch_method = 1;
                                        break;
                                    default:
                                        cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                                        cerr << "'s1' will be applied instead." ;
                                        ch_method = 1;
                                        break;
                                }
                            }

                            if(ChooseMethod == "Automatic")
                                ch_method = 1;

                            if(ch_method == 0)
                            {
                                string error = "Bad option value: " + ChooseMethod + ".\n"
                                "Available options: s1, a (s2 only for spindown = 1).\n";
                                throw domain_error(error);
                            }
                            ///cout << "test 1" << endl;
                            double densityDS;
                            vector<double> sgridDS;
                            if(ch_method == 1)
                            {
                                sgridDS = gs1DS.grid_prim(c0, nfft, data_length);
                                densityDS = gs1DS.DensityS1DS::density(sgridDS);
                            }
                            //#mb
                            for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                            {
                                vector<double> sgridDS2;
                                if( conv )
                                {
                                    if( ch_method == 1)
                                        sgridDS2 = gs1DS.convert(c0, xi, sgridDS);
                                }
                                else{
                                    sgridDS2 = sgridDS;
                                }

                                if(d_chop)
                                    num::chop(sgridDS2);


                                //#mb writing the fftpad to grid.bin
                                int fftpad;
                                fftpad = (log2(nfft) - ceil(log2(data_length)) + 1);

                                cout << "fftpad: " << fftpad << endl ;

                                gridbin.write((char*) &fftpad, sizeof(int));

                                unsigned int datal = data_length;
                                /// 'datal' is this same as 'data_length'
                                M.resize(dim*dim);
                                Mn.resize(dim*dim);
                                for(unsigned int i = 0; i<sgridDS2.size(); i++) {
                                    M[i] = sgridDS2[i]/datal;
                                    Mn[i] = sgridDS2[i];
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
                                    cout << densityDS << "\n";
                                    for(unsigned int i = 0; i<sgridDS2.size(); i++)
                                    {
                                        cout << setw(15) << setprecision(11) << left << sgridDS2[i];
                                        if( (i+1)%dim==0 )
                                            cout << endl;
                                        else
                                            cout << "\t";
                                    }

                                    //fs.unget();
                                    //fs.close();
                                }
                                else{
                                    string error = "Can't open stream!\n";
                                    throw runtime_error(error);
                                }
                                //cout << name1 << endl;
                            }
                            cout << endl;
                        } // End of 'for(double c0=CovarianceMin ...'
                    } // End of if

                    if (SaveFisherMatrix == "True")
                    {
                        for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                        {
                            vector<double> fmDS = cp_FMDS->postrmf(xi);

                            if( cout.good() == true )
                            {
                                cout << scientific  << left;//<< setw(16) << setprecision(15) << fixed
                                for(unsigned int i = 0; i<fmDS.size(); i++)
                                {
                                    cout << fmDS[i];
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

                    /// Default Michal Bejger's options
                    if(SaveGrid=="False" && SaveFisherMatrix=="False" && SaveDensityOfCovering=="False"){
                        for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                        {
                            vector<double> fmds = cp_FMDS->postrmf(xi);

                            if( cout.good() == true )
                            {
                                cout << setw(13) << setprecision(12) << fixed << left;

                                cout << "Fisher matrix:" << endl;
                                for(unsigned int i = 0; i<fmds.size(); i++)
                                {

                                    //#mb Printout of Fisher matrix
                                    cout << scientific << fmds[i];
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

                        //#mb Fisher matrix written to the grid.bin file
                        vector<double> fish_mat = cp_FMDS->postrmf(InitialTimeMin);
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

                    delete cp_FMDS;

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
            TypeSearch=vov[23].m_s;     /// v0.3.01
            DataS2Use=vov[24].m_s;
            DataS2SaveBest=vov[25].m_s;
            Convert=vov[26].m_s;
            ObserveTimeSymmetric=vov[27].m_s;   /// v0.3.01
            Quiet=vov[28].m_s;
            SayHello=vov[29].m_s;

            sizeCo = vov[2].m_s.size(); /// Remove sizeC0 -> .m_s is a string!
            sizeIt = vov[5].m_s.size()+1;
            sizeNEL = to_string(NfftExpLength).size();
            sizeDLn = to_string(DataLength).size();
            snel = to_string(NfftExpLength).substr(0,sizeNEL);
            sdln = to_string(DataLength).substr(0,sizeDLn);

            if(TypeSearch == "AllSky")
                TypeSearch = "All-Sky"; /// Change in misspell case

            if(TypeSearch == "All-Sky")
                dim = Spindown + 3;

            if( TypeSearch=="Directed")
                dim = Spindown + 1;

            if(ObserveTimeSymmetric=="False"){
                InitialTimeMin-=0.5;
                InitialTimeMax-=0.5;
            }

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

            if(TypeSearch == "All-Sky" && Spindown == 1){
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
                        if( fsd.good() == true ){
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
                            ch_method = 2;

                        if(ChooseMethod == "Automatic")
                            ch_method = method(c0, nfft, data_length, Spindown);

                        if(ch_method == 0)
                        {
                            string error = "Bad option value: " + ChooseMethod + ".\n"
                            "Available options: s1, s2, a.\n";
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
                else{

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
                                    ch_method = 2;

                                if(ChooseMethod == "Automatic")
                                    ch_method = method(c0, nfft, DataLength, Spindown);

                                if(ch_method == 0)
                                {
                                    string error = "Bad option value: " + ChooseMethod + ".\n"
                                    "Available options: s1, s2, a.\n";
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
            if(TypeSearch == "All-Sky" && Spindown == 2){
                /// 06 All-Sky. Spindown = 2. Options read from file
                if(DataLength==0)    // DataLength == 0 means DataLength should be equal to ephemeris length
                {
                    vector<string> paths_to_ephemeris;
                    paths_to_ephemeris.push_back(PathSSB);
                    paths_to_ephemeris.push_back(PathDet);
                    paths_to_ephemeris.push_back(PathDetSSB);

                    ReadEphemeris RE = ReadEphemeris(paths_to_ephemeris);
                    FisherRM const *const cp_FM = new FisherRM(RE.get_ephemeris1(), RE.get_ephemeris2(), Spindown);

                    GridS1 gs1 = GridS1(cp_FM, Spindown);

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
                        if( fsd.good() == true ){
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

                        if(ChooseMethod == "s2"){
                            cerr << "Algorithm 's2' is implemented only for spindown = 1.\n";
                            cerr << "'s1' will be applied instead.";
                            ch_method = 1;
                        }

                        if(ChooseMethod == "Automatic")
                            ch_method = 1;

                        if(ch_method == 0)
                        {
                            string error = "Bad option value: " + ChooseMethod + ".\n"
                            "Available options: s1, a (s2 only for spindown = 1).\n";
                            throw domain_error(error);
                        }

                        double density;
                        vector<double> sgrid;
                        if(ch_method == 1)
                        {
                            sgrid = gs1.grid_prim(c0, nfft, data_length);
                            density = gs1.DensityS1::density(sgrid);
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
                else{

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
                                    cerr << "Algorithm 's2' is implemented only for spindown = 1.\n";
                                    cerr << "'s1' will be applied instead." ;
                                    ch_method = 1;
                                }

                                if(ChooseMethod == "Automatic")
                                    ch_method = 1;

                                if(ch_method == 0)
                                {
                                    string error = "Bad option value: " + ChooseMethod + ".\n"
                                    "Available options: s1, a (s2 only for spindown = 1).\n";
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
            if(TypeSearch == "Directed" && Spindown == 1){
                /// 07 Directed. Spindown = 1. Options read from file
                if(DataLength==0){    // DataLength == 0 means DataLength should be equal to ephemeris length
                    DataLength = DataLengthDS;
                }

                FisherRMDS const *const cp_FMDS = new FisherRMDS(Spindown);

                GridS1DS gs1DS = GridS1DS(cp_FMDS, Spindown);
                GridS2DS gs2DS = GridS2DS(cp_FMDS);

                unsigned int data_length = DataLength;
                size_t dL_Eph_size = to_string(data_length).size();
                string sd_Eph= to_string(data_length).substr(0,dL_Eph_size);
                string named = FilePatternDoC;
                stringmanip::sreplace(named, "%N", snel);
                stringmanip::sreplace(named, "%D", sd_Eph);
                ofstream fsd;

                if (SaveDensityOfCovering == "True")
                {
                    fsd.open( PathSave + "/" + named, fstream::out | fstream::trunc);
                    if( fsd.good() == true ){
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
                        ch_method = 2;

                    if(ChooseMethod == "Automatic")
                        ch_method = methodDS(c0, nfft, data_length, Spindown);

                    if(ch_method == 0)
                    {
                        string error = "Bad option value: " + ChooseMethod + ".\n"
                        "Available options: s1, s2, a.\n";
                        throw domain_error(error);
                    }

                    double densityDS;
                    vector<double> sgridDS;
                    if(ch_method == 1)
                    {
                        sgridDS = gs1DS.grid_prim(c0, nfft, data_length);
                        densityDS = gs1DS.DensityS1DS::density(sgridDS);
                    }

                    if(ch_method == 2)
                    {
                        sgridDS = gs2DS.grid_prim(c0, nfft, data_length, mode);
                        densityDS = gs2DS.density(sgridDS);
                    }

                    if (SaveDensityOfCovering == "True")
                    {
                        fsd.fill('0');
                        fsd.width(4);
                        fsd << left << c0 << " " << setw(13) << setprecision(12)
                        << densityDS << "\n";
                    }

                    if (SaveGrid == "True")
                    {
                        string nameg = FilePatternGrid;
                        string sc0 = to_string(c0).substr(0,sizeCo); ///! sizeCo
                        stringmanip::sreplace(nameg, "%C", sc0);
                        stringmanip::sreplace(nameg, "%N", snel);

                        for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                        {
                            vector<double> sgridDS2;
                            if( conv )
                            {
                                if( ch_method == 1)
                                    sgridDS2 = gs1DS.convert(c0, xi, sgridDS);

                                if( ch_method == 2)
                                    sgridDS2 = gs2DS.convert(c0, xi, sgridDS);
                            }
                            else{
                                sgridDS2 = sgridDS;
                            }

                            if(d_chop)
                                num::chop(sgridDS2);

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
                                fs << densityDS << "\n";
                                for(unsigned int i = 0; i<sgridDS2.size(); i++)
                                {
                                    fs << setw(15) << setprecision(12) << left << sgridDS2[i];
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
                    vector<double> fmDS = cp_FMDS->postrmf(xi);

                    fstream fs;
                    fs.open( PathSave + "/" + name, fstream::out | fstream::trunc);
                    if( fs.good() == true )
                    {
                        fs << setw(16) << setprecision(15) << fixed << left;
                        for(unsigned int i = 0; i<fmDS.size(); i++)
                        {
                            fs << fmDS[i];
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

                delete cp_FMDS;
            }
            if(TypeSearch == "Directed" && Spindown >= 2){
                /// 08 Directed. Spindown >= 2. Options read from file
                if(DataLength==0){    // DataLength == 0 means DataLength should be equal to ephemeris length
                    DataLength = DataLengthDS;
                }
                FisherRMDS const *const cp_FMDS = new FisherRMDS(Spindown);

                GridS1DS gs1DS = GridS1DS(cp_FMDS, Spindown);

                unsigned int data_length = DataLength;
                size_t dL_Eph_size = to_string(data_length).size();
                string sd_Eph= to_string(data_length).substr(0,dL_Eph_size);
                string named = FilePatternDoC;
                stringmanip::sreplace(named, "%N", snel);
                stringmanip::sreplace(named, "%D", sd_Eph);
                ofstream fsd;

                if (SaveDensityOfCovering == "True")
                {
                    fsd.open( PathSave + "/" + named, fstream::out | fstream::trunc);
                    if( fsd.good() == true ){
                        fsd << "Covariance, Density of covering\n";
                    }
                    else{
                        string error = "Can't open: " + named;
                        throw runtime_error(error);
                    }
                }
                /// - 8 -
                int ch_method = 0;
                if(ChooseMethod == "s1")
                    ch_method = 1;

                if(ChooseMethod == "s2")
                {
                    switch(Spindown)
                    {
                        case 2:
                            cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                            cerr << "'s1' will be applied instead." ;
                            ch_method = 1;
                            break;
                        case 3:
                            cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                            cerr << "'s1' will be applied instead." ;
                            ch_method = 1;
                            break;
                        default:
                            cerr << "Algorithm 's2' is implemented only for spindown=1.\n";
                            cerr << "'s1' will be applied instead." ;
                            ch_method = 1;
                            break;
                    }
                }

                if(ChooseMethod == "Automatic")
                    ch_method = 1;

                if(ch_method == 0)
                {
                    string error = "Bad option value: " + ChooseMethod + ".\n"
                    "Available options: s1, a (s2 only for spindown = 1).\n";
                    throw domain_error(error);
                }
                /// = 8 =
                if (SaveGrid == "True" || SaveDensityOfCovering == "True")
                for(double c0=CovarianceMin; c0<=CovarianceMax; c0+=CovarianceStep)
                {
                    //cout << "[c0=" << c0 <<"]\n";

                    double densityDS;
                    vector<double> sgridDS;
                    if(ch_method == 1)
                    {
                        sgridDS = gs1DS.grid_prim(c0, nfft, data_length);
                        densityDS = gs1DS.DensityS1DS::density(sgridDS);
                    }

                    if (SaveDensityOfCovering == "True")
                    {
                        fsd.fill('0');
                        fsd.width(4);
                        fsd << left << c0 << " " << setw(13) << setprecision(12)
                        << densityDS << "\n";
                    }

                    if (SaveGrid == "True")
                    {
                        string nameg = FilePatternGrid;
                        string sc0 = to_string(c0).substr(0,sizeCo); ///! sizeCo
                        stringmanip::sreplace(nameg, "%C", sc0);
                        stringmanip::sreplace(nameg, "%N", snel);

                        for(double xi=InitialTimeMin; xi<=InitialTimeMax; xi+=InitialTimeStep)
                        {
                            vector<double> sgridDS2;
                            if( conv )
                            {
                                if( ch_method == 1)
                                   sgridDS2 = gs1DS.convert(c0, xi, sgridDS);
                            }
                            else{
                                sgridDS2 = sgridDS;
                            }

                            if(d_chop)
                                num::chop(sgridDS2);

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
                                fs << densityDS << "\n";
                                for(unsigned int i = 0; i<sgridDS2.size(); i++)
                                {
                                    fs << setw(15) << setprecision(12) << left << sgridDS2[i];
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
                    vector<double> fmDS = cp_FMDS->postrmf(xi);

                    fstream fs;
                    fs.open( PathSave + "/" + name, fstream::out | fstream::trunc);
                    if( fs.good() == true )
                    {
                        fs << setw(16) << setprecision(15) << fixed << left;
                        for(unsigned int i = 0; i<fmDS.size(); i++)
                        {
                            fs << fmDS[i];
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

                delete cp_FMDS;
            }
        }
    }// End of 'try'
    catch (const runtime_error& e){ cout << e.what() << endl; }
    catch (const domain_error& e) { cout << e.what() << endl; }
    catch (...) { cout << "Error unknown type"; }

    return 0;
}


/// Function to choose in automatic way best type of algorithm for grid generation
/// (for All-Sky only).
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

/// Function to choose in automatic way best type of algorithm for grid generation
/// (for Directed searches only).
int methodDS(double c0, unsigned int nfft, unsigned int data_length, unsigned int spindown)
{
    int s = 0;  // grid: s=1 (grid s1), s=2 (grid s2, 2D case only: spindown = 1).
    switch(spindown){
        case 1:
            {
                double threshold = 3.64744;
                double x = num::delta_omega_zero_prim(c0, nfft, data_length);

                if(x<threshold){
                    s=2;
                }
                else{
                    s=1;
                }
            }
            break;
        case 2:
            s=1;
            break;
        case 3:
            s=1;
            break;
    }

    return s;
}

bool if_exist(map<string, vector<string> >& map_name, string command)
{
    return map_name.find(command) != map_name.end();
}


