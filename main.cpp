///////////////////////////////////////
// Name:        main.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     13/10/2016
///////////////////////////////////////

/// File 'main.cpp' still need work
/// (mainly reorganizing code and more testing)
/// This (alpha) version seems to work good

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

int method(double c0, unsigned int nfft, unsigned int data_length);
bool if_exist(map<string, vector<string>>&, string);


int main(int argc, char* argv[])
{
    //HWND console = GetConsoleWindow();
    //RECT r;
    //GetWindowRect(console, &r); //stores the console's current dimensions

    //MoveWindow(console, r.left, r.top, 1000, 600, TRUE); // 800 width, 100 height

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
        options_available.push_back("-p");  // print g - grid, f - fisher matrix, d - density of covering, a - all
        options_available.push_back("-ch"); // chop (set to zero if number is less than num::epsilon)
        options_available.push_back("-cv"); // convert grid from hyper-sphere to hyper-ellipsoid space
        options_available.push_back("-v");  // program version
        //options_available.push_back("-s");  // stop program before close
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
        options_available_long.push_back("--print");        // print g - grid, f - fisher matrix, d - density of covering, a - all
        options_available_long.push_back("--chop");         // chop (set to zero if number is less than num::epsilon)
        options_available_long.push_back("--convert");      // convert grid from hyper-sphere to hyper-ellipsoid space
        options_available_long.push_back("--version");      // program version
        //options_available_long.push_back("--stop");         // stop program before close
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
        /*
            {
            string current = commands[i];
            string found = string(1, current.at(0));
            if( (found=="-" && !(current.find_first_of(".0123456789") != std::string::npos))
               || (current.find("--")!=std::string::npos && !(current.find_first_of(".0123456789") != std::string::npos)) )
                position.push_back(i);
            }
        */

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
                    && options_available[i]!="-c" && options_available[i]!="-m")
                        options_found[options_available[i]].push_back("");
            }

        }
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
        vector<string> text_of_help;
        text_of_help.push_back("SYNOPSIS\n");
        text_of_help.push_back("\t ./gg [-d|--directory <path>] [-c|--covariance or -m|--match <min> <max> <step>]\n");
        text_of_help.push_back("\t      [-i|--initial <min> <max> <step>] [-a|--algorithm <type>] [-n|--nfft <int>]\n");
        text_of_help.push_back("\t      [-nd|--ndata <int>] [-na|-nalpha <int>] [-nr|--nradius <int>]\n");
        text_of_help.push_back("\t      [-cv|--convert <bool>] [-ch|--chop <bool>] [-p|--print <string>]\n");
        text_of_help.push_back("\n");   //To sign empty line
        text_of_help.push_back("\t ./gg [-h|--help] [-ht|--helptxt] [-v|--version] [-aa|--author]\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("DESCRIPTION\n");
        text_of_help.push_back("\t GridsGenerator (GG) is designated to be used in all-sky narrow-band\n");
        text_of_help.push_back("\t searches of continuous gravitational waves. Program allow to:\n");
        text_of_help.push_back("\t - generate efficient grid(s) for chosen initial time of observation (1).\n");
        text_of_help.push_back("\t - generate reduced Fisher matrix for chosen initial time of observation (1),\n");
        text_of_help.push_back("\t - generate density of covering (2).\n");
        text_of_help.push_back("\t (1) To get result, ephemeris must be provided.\n");
        text_of_help.push_back("\t (2) Result can be get without ephemeris (for more information see flags: -nd\n");
        text_of_help.push_back("\t (--ndata)).\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("FLAGS\n");
        text_of_help.push_back("   Flags with (optional) argument(s):\n");
        text_of_help.push_back("\t -c or --covariance <min> <max> <step>\n");
        text_of_help.push_back("\t\t Covariance. Flag -c is required (even without argument) to get result.\n");
        text_of_help.push_back("\t\t <min> - minimum value of covariance but not less than 0;\n");
        text_of_help.push_back("\t\t default set: 0.75.\n");
        text_of_help.push_back("\t\t <max> - optional maximum value of covariance but less than 1.\n");
        text_of_help.push_back("\t\t <step> - optional step value of covariance;\n");
        text_of_help.push_back("\t\t default set: 0.01.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t # Calculation are preform only in two cases:\n");
        text_of_help.push_back("\t\t # 1. No flag are provided. Sets are read from file 'gg.ini'.\n");
        text_of_help.push_back("\t\t # Result(s) is(are) printed to file(s) - see configuration file: 'gg.ini'.\n");
        text_of_help.push_back("\t\t # 2. Flag -c (--covariance) or -m (--match) is used.\n");
        text_of_help.push_back("\t\t # Result(s) is(are) printed to tty.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -m or --match <min> <max> <step>\n");
        text_of_help.push_back("\t\t Minimal match (MM^2 == covariance).\n");//Minimum versus Minimal ?
        text_of_help.push_back("\t\t <min> - minimum value of minimal match but not less than 0;\n");
        text_of_help.push_back("\t\t default set: MM^2 = 0.75.\n");
        text_of_help.push_back("\t\t <max> - optional maximum value of minimal match but less than 1.\n");
        text_of_help.push_back("\t\t <step> - optional step value of minimal match;\n");
        text_of_help.push_back("\t\t default set: 0.1.\n");
        text_of_help.push_back("\t\t ## If flags -c (--covariance) and -m (--match) provided simultaneously,\n");
        text_of_help.push_back("\t\t ## program will read options from -c flag only.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -d or --directory <path>\n");
        text_of_help.push_back("\t\t Directory containing ephemeris (need to contain binary files:\n");
        text_of_help.push_back("\t\t 'rSSB.bin', 'rDet.bin', 'rDetSSB.bin').\n");
        text_of_help.push_back("\t\t <path> - take path to directory.\n");
        text_of_help.push_back("\t\t E.g. -d 001: directory '001' located inside folder with GridsGenerator.\n");
        text_of_help.push_back("\t\t If directory is not provided, program will try to find ephemeris\n");
        text_of_help.push_back("\t\t in directory with GridsGenerator.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -i or --initial <min> <max> <step>\n");
        text_of_help.push_back("\t\t Initial time of observation.\n");
        text_of_help.push_back("\t\t <min> - minimum value of initial time;\n");
        text_of_help.push_back("\t\t default set: 0.5.\n");
        text_of_help.push_back("\t\t <max> - optional maximum value of initial time.\n");
        text_of_help.push_back("\t\t <step> - optional step value of minimal match;\n");
        text_of_help.push_back("\t\t if not provided step will be set on step = max - min.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -a or --algorithm <type>\n");
        text_of_help.push_back("\t\t Algorithm type to choose.\n");
        text_of_help.push_back("\t\t <type> take option: s1, s2, a. Algorithms:\n");
        text_of_help.push_back("\t\t s1 - based on fully analytic formula,\n");
        text_of_help.push_back("\t\t s2 - partially numeric formula.\n");
        text_of_help.push_back("\t\t Accuracy for algorithm 's2' depended on -na (--nalpha) and -nr (--nradius)\n");
        text_of_help.push_back("\t\t flags.\n");
        text_of_help.push_back("\t\t a - automatically choose algorithm 's1' or 's2'. Use this argument to allow\n");
        text_of_help.push_back("\t\t GridsGenerator to decide which algorithm (for given parameters) should be\n");
        text_of_help.push_back("\t\t used to get grid with better density of covering.\n");
        text_of_help.push_back("\t\t Information about implemented algorithms can be found in article:\n");
        text_of_help.push_back("\t\t http://dx.doi.org/10.1088/0264-9381/32/14/145014\n");
        text_of_help.push_back("\t\t Default set: a.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -n or --nfft <int>\n");
        text_of_help.push_back("\t\t Number of Fourier bins.\n");
        text_of_help.push_back("\t\t <int> take positive integer number without zero (exponent).\n");
        text_of_help.push_back("\t\t E.g. to get grid destined to work with discreet\n");
        text_of_help.push_back("\t\t Fourier transform (DFT, FFT) with length 1024 put same exponent of 2: 10\n");
        text_of_help.push_back("\t\t (1024 == 2^10).\n");
        text_of_help.push_back("\t\t Default set: 20 (1048576 = 2^20).\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -nd or --ndata <int>\n");
        text_of_help.push_back("\t\t Number of data (data length collected by detector, equal to length\n");
        text_of_help.push_back("\t\t of ephemeris).\n");
        text_of_help.push_back("\t\t <int> take positive integer number including zero.\n");
        text_of_help.push_back("\t\t If <int> is set to zero data length will be read from ephemeris*.\n");
        text_of_help.push_back("\t\t If <int> is bigger than zero program will able to obtain density of\n");
        text_of_help.push_back("\t\t covering only**.\n");
        text_of_help.push_back("\t\t Default set: 0.\n");
        text_of_help.push_back("\t\t * With this set (-nd 0 or -ndata 0) density is obtained without using\n");
        text_of_help.push_back("\t\t information stored in file 'dataS2.txt'.\n");
        text_of_help.push_back("\t\t ### With this set density of covering for algorithm 's2' is always\n");
        text_of_help.push_back("\t\t ### obtaining with maximal accuracy, depending only from flags: -na, -nr.\n");
        text_of_help.push_back("\t\t ### Density of covering for algorithm 's1' is always obtained with maximal\n");
        text_of_help.push_back("\t\t ### accuracy regardless to sets of -nd flag.\n");
        text_of_help.push_back("\t\t ** Ephemeris are not used and because of that grid(s) and reduced Fisher\n");
        text_of_help.push_back("\t\t matrix can't be obtained. Density of covering for algorithm 's2' will be\n");
        text_of_help.push_back("\t\t be obtained in approximated way based on information*** collected (in \n");
        text_of_help.push_back("\t\t 'dataS2.txt' file) during previous runs with ephemeris. File 'dataS2.txt'\n");
        text_of_help.push_back("\t\t collect information about coverings only for algorithm 's2'.\n");
        text_of_help.push_back("\t\t Program inform which algorithm has been used only when work without\n");
        text_of_help.push_back("\t\t ephemeris (<int> bigger that zero).\n");
        text_of_help.push_back("\t\t E.g. (without ephemeris) -nd 344656:\n");
        text_of_help.push_back("\t\t Covariance, Density, Algorithm\n");
        text_of_help.push_back("\t\t 0.86      1.82      s2\n");
        text_of_help.push_back("\t\t 0.87      1.8494      s1\n");
        text_of_help.push_back("\t\t 0.88      1.77685      s1\n");
        text_of_help.push_back("\t\t E.g. (data length taken from ephemeris) -nd 0:\n");
        text_of_help.push_back("\t\t Covariance, Density of covering\n");
        text_of_help.push_back("\t\t 0.86 1.82299890516\n");
        text_of_help.push_back("\t\t 0.87 1.84940429814\n");
        text_of_help.push_back("\t\t 0.88 1.77685017541\n");
        text_of_help.push_back("\t\t *** Information stored in file 'dataS2.txt' are independent of ephemeris.\n");
        text_of_help.push_back("\t\t Based on this information density of covering can be obtained fast but\n");
        text_of_help.push_back("\t\t in approximated way (generally speaking more collected data allows\n");
        text_of_help.push_back("\t\t to get more accurate results).\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -na or --nalpha <int>\n");
        text_of_help.push_back("\t\t Number of loops (incrementation deep) in root finding algorithm.\n");
        text_of_help.push_back("\t\t <int> take positive integer number without zero.\n");
        text_of_help.push_back("\t\t This flag affect only on 's2' algorithm.\n");
        text_of_help.push_back("\t\t Default set: 35.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -nr or --nradius <int>\n");
        text_of_help.push_back("\t\t Number of loops in covering radius (deep hole) finding algorithm.\n");
        text_of_help.push_back("\t\t <int> take positive integer number without zero.\n");
        text_of_help.push_back("\t\t This flag affect only on 's2' algorithm.\n");
        text_of_help.push_back("\t\t Default set: 20.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -cv or --convert <bool>\n");
        text_of_help.push_back("\t\t Convert grid from space with hyper-sphere to hyper-ellipsoid space.\n");
        text_of_help.push_back("\t\t <bool> take argument: t (true), f (false).\n");
        text_of_help.push_back("\t\t Default set: t.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -ch or --chop <bool>\n");
        text_of_help.push_back("\t\t Chop result to zero if value of result is smaller than |10^-12|.\n");
        text_of_help.push_back("\t\t <bool> take argument: t (true), f (false).\n");
        text_of_help.push_back("\t\t Default set: f.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -p or --print <string>\n");
        text_of_help.push_back("\t\t Print result(s).\n");
        text_of_help.push_back("\t\t <string> take argument(s): d, f, g.\n");
        text_of_help.push_back("\t\t d - density of covering,\n");
        text_of_help.push_back("\t\t f - Fisher reduced matrix,\n");
        text_of_help.push_back("\t\t g - grid.\n");
        text_of_help.push_back("\t\t Option can be joined (order is not important), e.g. df, dg, fg, dfg.\n");
        text_of_help.push_back("\t\t Default set: g.\n");//????????
        text_of_help.push_back("\t\t #### If argument of flag -nd is set to be bigger than zero, flag -p\n");
        text_of_help.push_back("\t\t #### can accept only 'd' argument.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("   Flags with no arguments:\n");
        text_of_help.push_back("\t -h or --help\n");
        text_of_help.push_back("\t\t Display this help.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -ht or --helptxt\n");
        text_of_help.push_back("\t\t Print this help to text file.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -v or --version\n");
        text_of_help.push_back("\t\t Display information about version of GridsGenerator.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t -aa or --author\n");
        text_of_help.push_back("\t\t About author(s). Display information about developer(s).\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("EXAMPLES\n");
        text_of_help.push_back("\t ./gg\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t Settings are read from configuration file: 'gg.ini'.\n");
        text_of_help.push_back("\t\t Result(s) is(are) printed to file(s) denoted in 'gg.ini'\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t ./gg -c\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t GridsGenerator work with defaults sets (covariance set on 0.75,\n");
        text_of_help.push_back("\t\t initial time on 0.5). Result (grid) is printed to tty.\n");
        text_of_help.push_back("\t\t Ephemeris located at directory with GridsGenerator.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t ./gg -c 0.7 -i -10\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t Generate grid for covariance set to 0.7 and initial time set on -10.\n");
        text_of_help.push_back("\t\t Ephemeris located at directory with GridsGenerator.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t ./gg -c 0.7 -i -10 >> result.txt\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t Result is appended to end of the file 'result.txt'.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t ./gg -c 0.7 0.95 -i -10 10 0.1 -n 21 -nd 344656 -p d\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t Generate density of covering (only) without using ephemeris (344656!=0).\n");
        text_of_help.push_back("\t\t GridsGenerator print density obtained for:\n");
        text_of_help.push_back("\t\t - covariance from 0.7 to 0.95 with default step = 0.01,\n");
        text_of_help.push_back("\t\t - initial time from -10 to 10 with step 0.1,\n");
        text_of_help.push_back("\t\t - Fourier length (nfft) set to 2^21 (2097152),\n");
        text_of_help.push_back("\t\t - data length set on 344656.\n");
        text_of_help.push_back("\t\t Ephemeris are located at directory with GridsGenerator.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t ./gg -d 001 -c 0.7 0.95 -i -10 10 0.1 -n 21 -p dfg\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t GG print densities, fisher matrices, and grids.\n");
        text_of_help.push_back("\t\t Ephemeris are located at directory '001'.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t ./gg -d 001 -c 0.7 0.95 0.001 -i -10 10 0.1 -n 21 -p dfg\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t Covariance step is set to 0.001,\n");
        text_of_help.push_back("\t\t Nfft length is set to 2^21 (2097152).\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t ./gg -d 001 -c 0.7 0.95 0.001 -i -10 10 0.1 -a a -na 40 -nr 30 -cv f\n");
        text_of_help.push_back("\t      -ch t -p dfg\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t Nfft length is set to default value 2^20 (1048576).\n");
        text_of_help.push_back("\t\t GridsGenerator decide which algorithm should be used.\n");
        text_of_help.push_back("\t\t Number of loops in root finding algorithm is set to 40.\n");
        text_of_help.push_back("\t\t Number of loops in covering radius finding algorithm is set to 30.\n");
        text_of_help.push_back("\t\t Grids are not converted to hyper-ellipsoid space.\n");
        text_of_help.push_back("\t\t All values are chopped to zero if are smaller than |10^-12|.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t ./gg -h -v -aa\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t Display help, version of GridsGenerator and show info about author.\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\t\t # Note that GridsGenerator provide result more effective (faster) in way:\n");
        text_of_help.push_back("\t ./gg -c 0.70 0.73 -d 001\n");
        text_of_help.push_back("\t\t # than when is running several times for the same ephemeris:\n");
        text_of_help.push_back("\t ./gg -c 0.70 -d 001\n");
        text_of_help.push_back("\t ./gg -c 0.71 -d 001\n");
        text_of_help.push_back("\t ./gg -c 0.72 -d 001\n");
        text_of_help.push_back("\t ./gg -c 0.73 -d 001\n");
        text_of_help.push_back("\t\t # That clue apply to both: covariance (-c)/minimal match (-m) and\n");
        text_of_help.push_back("\t\t # initial time (-i).\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("CONTACT\n");
        text_of_help.push_back("\t **********************************************************\n");
        text_of_help.push_back("\t * POLGRAW Group                                          *\n");
        text_of_help.push_back("\t * University of Bialystok, Faculty of Physics            *\n");
        text_of_help.push_back("\t * Developer: Andrzej Pisarski <a.pisarski at uwb.edu.pl> *\n");
        text_of_help.push_back("\t * Tutor: Piotr Jaranowski <p.jaranowski at uwb.edu.pl>   *\n");
        text_of_help.push_back("\t **********************************************************\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");
        text_of_help.push_back("\n");


        vector<string> text_of_version;
        text_of_version.push_back("\t *******************************************");
        text_of_version.push_back("\t * Build 0.2.27 (alpha).                   *");
        text_of_version.push_back("\t * Copyright 2015 Andrzej Pisarski.        *");
        text_of_version.push_back("\t * License of GridsGenerator: CC-BY-NC-ND. *");
        text_of_version.push_back("\t * Additional credits to:                  *");
        text_of_version.push_back("\t * Piotr Jaranowski                        *");
        text_of_version.push_back("\t *******************************************");

        double CovarianceMin, CovarianceMax, CovarianceStep=0.01;
        double InitialTimeMin=0.5, InitialTimeMax, InitialTimeStep;
        unsigned int NfftExpLength=20, DataLength=0;
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
        //snel = to_string(NfftExpLength).substr(0,sizeNEL);
        // sdln = to_string(DataLength).substr(0,sizeDLn);

        //bool read_config_file = true;
        //bool make_calculation = true;

        //if( if_exist(options_found, "-c") )
        //    read_config_file = false;

        //if( if_exist(options_found, "-h") || if_exist(options_found, "-v") || if_exist(options_found, "-aa")  )
        //    make_calculation = false;

        if( commands.size() != 0) /// if flag in command line provided
        {
            ///Help
            if( if_exist(options_found, "-h") )
            {
                system("mode 93, 30");
                for(size_t k=0; k<=10; k++)
                    for(size_t i=0; i<24; i++)
                    {
                        if((24*k+i)<text_of_help.size())
                        {
                            cout<< text_of_help[24*k+i];
                            if(i==23)
                            {
                                if(k!=10)
                                    cout << "\n\nPage: " << k+1 << "/11" << ". Press Enter to see next page.\n" <<endl;
                                else
                                    cout << "\n\nPage: " << k+1 << "/11" << ". Press Enter to exit.\n" <<endl;
                                getchar();
                                system("cls");
                            }
                        }

                    }
            }

            ///Help print to text file
            if( if_exist(options_found, "-ht") )
            {
                    fstream fsh;
                    fsh.open( "help.txt", fstream::out | fstream::trunc);
                    if( fsh.good() == true )
                    {
                        fsh << left;
                        for(unsigned int i = 0; i<text_of_help.size(); i++)
                            fsh << text_of_help[i];

                        fsh.close();
                    }
                    else{
                        string error = "Can't open: 'help.txt.";
                        throw runtime_error(error);
                    }
            }


            ///Version
            if( if_exist(options_found, "-v") )
            {
                for(size_t i=0; i<text_of_version.size(); i++)
                    cout << text_of_version[i] << endl;
            }

            ///About author
            if( if_exist(options_found, "-aa") )
            {
                cout << "\t **********************************************************\n"
                "\t * POLGRAW Group                                          *\n"
                "\t * University of Bialystok, Faculty of Physics            *\n"
                "\t * Developer: Andrzej Pisarski <a.pisarski at uwb.edu.pl> *\n"
                "\t * Tutor: Piotr Jaranowski <p.jaranowski at uwb.edu.pl>   *\n"
                "\t **********************************************************\n" << endl;
            }
/*
            bool c_exist = false;
            if( if_exist(options_found, "-c") )
            {
                map<string, vector<string> >::const_iterator it = options_found.find("-c");
                vector<string>::const_iterator line_it = it->second.begin();

                if(*line_it!="")
                    c_exist=true;

            }
*/

            if( if_exist(options_found, "-c") )
            {
                /// Directory; PathSSB, PathDet, PathDetSSB
                //Default path for linux (GG location; POSIX) not implement yet.
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

                    if(*(it->second.begin())!="")
                    {
                        //cout << "CovarianceMin= " << CovarianceMin << endl;
                        //cout << setw(13) << setprecision(12) << "CovarianceMax= " << CovarianceMax << endl;
                        //cout << "CovarianceStep= " << CovarianceStep << endl;
                    }

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

                    if(*(it->second.begin())!="")
                    {
                        //cout << "InitialTimeMin= " << InitialTimeMin << endl;
                        //cout << "InitialTimeMax= " << InitialTimeMax << endl;
                        //cout << "InitialTimeStep= " << InitialTimeStep << endl;
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
                        
                        //#mb was unsigned int - warning: comparison of unsigned expression >= 0 is always true [-Wtype-limits] 
                        int tt = stoi(*line_it);
                        if(tt >= 0)
                            DataLength = tt;
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
                        unsigned int tt = stoi(*line_it);
                        if(tt > 0)
                            NAlpha = tt;
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
                        unsigned int tt = stoi(*line_it);
                        if(tt > 0)
                            NRadius = tt;
                        else{
                            string error = "NAlpha can not be less or equal zero!\n";
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
                    if(ts_end > 3)
                        ts_end = 3;

                    if(*line_it=="")
                        {SaveGrid="True"; no_match = false;}
                    else{
                        string temp = string((*line_it).begin(),(*line_it).begin()+ts_end); //[,)

                        if(temp.find_first_of("d")!=string::npos)
                            SaveDensityOfCovering="True"; no_match = false;

                        if(temp.find_first_of("f")!=string::npos)
                            SaveFisherMatrix="True"; no_match = false;

                        if(temp.find_first_of("g")!=string::npos)
                            SaveGrid="True"; no_match = false;
                    }

                    if(no_match){
                        string error = "Arguments for -p (--print) can be only: d (density of covering),\n"
                        "f (Fisher reduced matrix) and g (grid).\n";
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
string gridout = PathSave + "grid.bin";
cout << "grid.bin will be saved in " << gridout << endl ;

std::ofstream gridbin(gridout, std::ios::binary);

vector<double> M(16), Mn(16); 

                    vector<string> paths_to_ephemeris;
                    paths_to_ephemeris.push_back(PathSSB);
                    paths_to_ephemeris.push_back(PathDet);
                    paths_to_ephemeris.push_back(PathDetSSB);

                    ReadEphemeris RE = ReadEphemeris(paths_to_ephemeris);
                    FisherRM const *const cp_FM = new FisherRM(RE.get_ephemeris1(), RE.get_ephemeris2());

                    GridS1 gs1 = GridS1(cp_FM);
                    GridS2 gs2 = GridS2(cp_FM, NAlpha, NRadius);

                    unsigned int data_length = cp_FM->get_ephemeris_length();
                    //ofstream fsd;

                    if (SaveDensityOfCovering == "True")
                    {
                        if( cout.good() == true )
                        {
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
                            ch_method = method(c0, nfft, RE.get_length());

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
                            cout.fill('0');
                            cout.width(4);
                            cout << "Density of covering: " << left << c0 << " " << setw(13) << setprecision(12)
                            << density << "\n";
                        }

//                        if (SaveGrid == "True")
//                        {
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

for(unsigned int i = 0; i<sgrid2.size(); i++) { 
  M[i] = sgrid2[i]/datal; 
  Mn[i] = sgrid2[i]; 
} 

M[1] = M[1]/datal;
M[5] = M[5]/datal;
M[9] = M[9]/datal;
M[13] = M[13]/datal;

  cout << "Normalized grid matrix:" << endl;
    for(unsigned int i = 0; i<M.size(); i++) {
        cout << scientific << M[i];
             if( (i+1)%4==0 )
                    cout << "\n";
                         else
                                cout << "\t";
                                  }


gridbin.write((char*)&M[0], M.size()*sizeof(double));


                                if( cout.good() == true)
                                {   // fixed << left
//                                    cout << setw(15) << setprecision(12) << left;
//                                    cout << density << "\n";

                                    for(unsigned int i = 0; i<sgrid2.size(); i++)
                                    {
                                     
                                        //#mb Printout of grid matrix 
/*                                        cout << setw(15) << setprecision(12) << left 
                                        << scientific; 
                                        cout << sgrid2[i];
                                        if( (i+1)%4==0 )
                                                 cout << "\n";
                                        else
                                            cout << "\t";
*/                                            
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
//                        } // End of 'if(SaveGrid == "True")'
                    } // End of 'for(double c0=CovarianceMin ...'

                  //  if (SaveFisherMatrix == "True")
                  //  {
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
                                    if( (i+1)%4==0 )
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
                        }                        cout << endl;
                  //  }

//#mb Fisher matrix written to the grid.bin file 
vector<double> fish_mat = cp_FM->postrmf(InitialTimeMin); 
gridbin.write((char*)&fish_mat[0], fish_mat.size()*sizeof(double));

//#mb Mn matrix

  cout << "Grid matrix:" << endl;
  for(unsigned int i = 0; i<Mn.size(); i++) { 
    cout << scientific << Mn[i];
     if( (i+1)%4==0 )
       cout << "\n";
     else 
       cout << "\t";
  }

gridbin.write((char*)&Mn[0], Mn.size()*sizeof(double));
gridbin.close(); 

                    delete cp_FM;

                } // End of 'if(DataLength==0)'
                else
                {

                    if (SaveDensityOfCovering == "True")
                    {
                        //ofstream fs;
                        //fs.open( PathSave + "\\" + name, fstream::out | fstream::trunc);
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
                                    ch_method = method(c0, nfft, DataLength);

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
                                    DensityS1 ds1 = DensityS1();
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
            PathSSB=vov[10].m_s;
            PathDet=vov[11].m_s;
            PathDetSSB=vov[12].m_s;
            PathSave=vov[13].m_s;
            FilePatternGrid=vov[14].m_s;
            FilePatternFM=vov[15].m_s;
            FilePatternDoC=vov[16].m_s;
            SaveGrid=vov[17].m_s;
            SaveFisherMatrix=vov[18].m_s;
            SaveDensityOfCovering=vov[19].m_s;
            DataChop=vov[20].m_s;
            ChooseMethod=vov[21].m_s;
            DataS2Use=vov[22].m_s;
            DataS2SaveBest=vov[23].m_s;
            Convert=vov[24].m_s;
            Quiet=vov[25].m_s;
            SayHello=vov[26].m_s;

            sizeCo = vov[2].m_s.size(); /// Remove sizeC0 -> .m_s is a string !
            sizeIt = vov[5].m_s.size()+1;
            sizeNEL = to_string(NfftExpLength).size();
            sizeDLn = to_string(DataLength).size();
            snel = to_string(NfftExpLength).substr(0,sizeNEL);
            sdln = to_string(DataLength).substr(0,sizeDLn);

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
                for(size_t i=0; i<text_of_version.size(); i++)
                    cout << text_of_version[i] << endl;

                system("pause");
            }

            if(DataLength==0)    // DataLength == 0 means DataLength should be equal to ephemeris length
            {
                vector<string> paths_to_ephemeris;
                paths_to_ephemeris.push_back(PathSSB);
                paths_to_ephemeris.push_back(PathDet);
                paths_to_ephemeris.push_back(PathDetSSB);

                ReadEphemeris RE = ReadEphemeris(paths_to_ephemeris);
                FisherRM const *const cp_FM = new FisherRM(RE.get_ephemeris1(), RE.get_ephemeris2());

                GridS1 gs1 = GridS1(cp_FM);
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
                    fsd.open( PathSave + "\\" + named, fstream::out | fstream::trunc);
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
                        ch_method = 2;

                    if(ChooseMethod == "Automatic")
                        ch_method = method(c0, nfft, RE.get_length());

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
                            fs.open( PathSave + "\\" + name1, fstream::out | fstream::trunc);
                            if( fs.good() == true )
                            {   //fixed << left
                                fs << setw(15) << setprecision(12) << left;
                                fs << density << "\n";
                                for(unsigned int i = 0; i<sgrid2.size(); i++)
                                {
                                    fs << setw(15) << setprecision(12) << left << sgrid2[i];
                                    if( (i+1)%4==0 )
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
                    fs.open( PathSave + "\\" + name, fstream::out | fstream::trunc);
                    if( fs.good() == true )
                    {
                        fs << setw(13) << setprecision(12) << fixed << left;
                        for(unsigned int i = 0; i<fm.size(); i++)
                        {
                            fs << fm[i];
                            if( (i+1)%4==0 )
                            {
                                fs << "\n";
                            }
                            else
                            {
                                fs << "\t";
                            }
                        }

                        fs.close();
                    }
                    else{
                        string error = "Can't open: " + PathSave + "\\" + name;
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
                    fs.open( PathSave + "\\" + name, fstream::out | fstream::trunc);
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
                                ch_method = method(c0, nfft, DataLength);

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
                                DensityS1 ds1 = DensityS1();
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

int method(double c0, unsigned int nfft, unsigned int data_length)
{
    double sq2= sqrt(2);
    double sq2_min = sq2 - num::epsilon();
    double sq2_max = sq2 + num::epsilon();

    double x = num::delta_omega_zero_prim(c0, nfft, data_length);
    double s = 0;

    if(x<=sqrt(3))
    {
        if(x<1.64)
        {
            if( x>sq2_min && x<sq2_max ) //exception grid S1 == grid S2 = grid A4*
                s=1;
            else
                s=2;
        }
        else
            s=1;
    }
    else
    {
        if(x<1.81)
            s=2;
        else
            s=1;
    }

    return s;
}

bool if_exist(map<string, vector<string>>& map_name, string command)
{
    return map_name.find(command) != map_name.end();
}

