///////////////////////////////////////
// Name:        Manual.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     15/11/2015
// Modification:14/06/2020 A.Pisarski
///////////////////////////////////////

#include "Manual.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

Manual::Manual(std::string ver="", unsigned int pages): m_help(init_m_help()), m_version(init_m_version(ver)),m_pagesNr(pages){}

std::vector<std::string> Manual::help()
{
    return m_help;
}

std::vector<std::string> Manual::version()
{
    return m_version;
}

void Manual::save_help()
{
    std::fstream fsh;
    fsh.open( "help.txt", std::fstream::out | std::fstream::trunc);
    if( fsh.good() == true )
    {
        fsh << std::left;
        for(unsigned int i = 0; i<m_help.size(); i++)
            fsh << m_help[i];

        fsh.close();
    }
    else{
        std::string error = "Can't open: help.txt.";
        throw std::runtime_error(error);
    }
}

void Manual::print_help()
{
    //system("'\e[8;30;93t'");
    size_t current_page=1;
    while(current_page<=m_pagesNr)
    {
        for(size_t i=0; i<24; i++)
        {
            if((24*(current_page-1)+i)<m_help.size())
            {
                std::cout<< m_help[24*(current_page-1)+i];
                if(i==23)
                {
                    std::cout << "\n\nPage: " << std::setw(2) << current_page << "/" << m_pagesNr;
                    if(current_page!=m_pagesNr){
                        if(current_page==1)
                            std::cout << ". Press 'Enter' to see next page.\n" <<std::endl;
                        else
                            std::cout << ". Press 'Enter' to see next page; 'u' + 'Enter' to previous page.\n" <<std::endl;
                    }
                    else
                        std::cout << ". Press Enter to exit.\n" <<std::endl;
                }
            }
        }

	char sign = getchar();
        if(sign=='u'){
            if(current_page>1)
                current_page-=2; ///fflush(stdout); system("clear");
        }
	else if(sign=='q'){
	    break;
	}
        else{
            ++current_page;system("clear");
	}

    }
}

void Manual::print_version()
{
    for(size_t i=0; i<m_version.size(); i++)
        std::cout << m_version[i] << std::endl;

}

void Manual::print_author()
{
    std::cout << "\t **********************************************************\n"
    "\t * POLGRAW Group                                          *\n"
    "\t * University of Bialystok, Faculty of Physics            *\n"
    "\t * Developer: Andrzej Pisarski <a.pisarski at uwb.edu.pl> *\n"
    "\t * Tutor: Piotr Jaranowski <p.jaranowski at uwb.edu.pl>   *\n"
    "\t **********************************************************\n" << std::endl;
}

std::vector<std::string> Manual::init_m_help() const
{
    std::vector<std::string> text_of_help;
    text_of_help.reserve(270);

    text_of_help.push_back("SYNOPSIS\n");
    text_of_help.push_back("\t ./gg [-d|--directory <path>] [-c|--covariance or -m|--match <min> <max> <step>]\n");
    text_of_help.push_back("\t      [-i|--initial <min> <max> <step>] [-a|--algorithm <type>] [-n|--nfft <int>]\n");
    text_of_help.push_back("\t      [-nd|--ndata <int>] [-na|-nalpha <int>] [-nr|--nradius <int>]\n");
    text_of_help.push_back("\t      [-cv|--convert <bool>] [-ch|--chop <bool>] [-p|--print <string>]\n");
    text_of_help.push_back("\t      [-t|--type <type>] [-u|--unbiased <bool>]\n");
    text_of_help.push_back("\n");   //To sign empty line
    text_of_help.push_back("\t ./gg [-h|--help] [-ht|--helptxt] [-v|--version] [-aa|--author]\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("DESCRIPTION\n");
    text_of_help.push_back("\t GridsGenerator (GG) is designated to be used in All-Sky narrow-band and Directed\n");
    text_of_help.push_back("\t searches of continuous gravitational waves. Program allow to:\n");
    text_of_help.push_back("\t - generate efficient grid(s) for chosen initial time of observation (1).\n");
    text_of_help.push_back("\t - generate reduced Fisher matrix for chosen initial time of observation (1),\n");
    text_of_help.push_back("\t - generate density of covering (2).\n");
    text_of_help.push_back("\t (1) To get result, ephemeris must be provided (All-Sky); grids always have form of lower\n");
    text_of_help.push_back("\t triangular matrix.\n");
    text_of_help.push_back("\t (2) Result can be get without ephemeris (for more information see flags: -nd\n");
    text_of_help.push_back("\t (--ndata)).\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("FLAGS\n");
    text_of_help.push_back("   Flags with (optional) argument(s):\n");
    text_of_help.push_back("\t -a or --algorithm <type>\n");
    text_of_help.push_back("\t\t Algorithm type to choose.\n");
    text_of_help.push_back("\t\t <type> take option: s1, s2, a.\n");
    text_of_help.push_back("\t\t -- All-Sky searches algorithms: --\n");
    text_of_help.push_back("\t\t s1 - based on analytic formula (for spindowns: 1, 2),\n");
    text_of_help.push_back("\t\t s2 - partially numeric formula (only for spindown 1).\n");
    text_of_help.push_back("\t\t Accuracy for algorithm 's2' depended on -na (--nalpha) and -nr (--nradius)\n");
    text_of_help.push_back("\t\t flags.\n");
    text_of_help.push_back("\t\t -- Directed searches algorithms: --\n");
    text_of_help.push_back("\t\t s1 - based on analytic formula (for spindowns: 1, 2, 3),\n");
    text_of_help.push_back("\t\t s2 - based on fully analytic formula (only for spindown 1).\n");
    text_of_help.push_back("\t\t a - automatically choose algorithm 's1' or 's2'. Use this argument to allow\n");
    text_of_help.push_back("\t\t GridsGenerator to decide which algorithm (for given parameters) should be\n");
    text_of_help.push_back("\t\t used to get grid with better density of covering.\n");
    text_of_help.push_back("\t\t Default set: a.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t   >>> Options for All-Sky searches only <<<\n");
    text_of_help.push_back("\t -d or --directory <path>\n");
    text_of_help.push_back("\t\t Directory containing ephemeris (need to contain binary files:\n");
    text_of_help.push_back("\t\t 'rSSB.bin', 'rDet.bin', 'rDetSSB.bin').\n");
    text_of_help.push_back("\t\t <path> - take path to directory.\n");
    text_of_help.push_back("\t\t E.g. -d 001/: directory '001/' located inside folder with GridsGenerator.\n");
    text_of_help.push_back("\t\t If directory is not provided, program will try to find ephemeris\n");
    text_of_help.push_back("\t\t in directory with GridsGenerator.\n");
    text_of_help.push_back("\t\t To get this same effect you can also apply second form:\n");
    text_of_help.push_back("\t\t -d <path> <SegmentNo> <detectorX>. In this form\n");
    text_of_help.push_back("\t\t you should to provide only one folder -- name of detector.\n");
    text_of_help.push_back("\t\t To get result for more than one detector (with applying\n");
    text_of_help.push_back("\t\t sum of reduced fisher matrices):\n");
    text_of_help.push_back("\t\t -d <path> <SegmentNo> <detectorX> ... <detectorZ> <band>\n");
    text_of_help.push_back("\t\t Provide at least two or more detectors.\n");
    text_of_help.push_back("\t\t E.g. -d /home/ligo/ 001 H1 L1 0371; -d /home/ligo/ 001 H1 L1 V1 0371.\n");
    text_of_help.push_back("\t\t To read data from file different than default 'xdats', provide appropriate\n");
    text_of_help.push_back("\t\t name of the file on the last position:\n");
    text_of_help.push_back("\t\t -d /home/ligo/ 001 H1 L1 0371 xdatsc_001_0067.bin\n");
    text_of_help.push_back("\t\t or more simpler form:\n");
    text_of_help.push_back("\t\t -d /home/ligo/ 001 H1 L1 0371 xdatsc\n");
    text_of_help.push_back("\t\t Data file name need to have prefix 'xdatas'.\n");
    text_of_help.push_back("\t\t Accepted detectors names: H1, L1, V1.\n");
    text_of_help.push_back("\t\t Sum of Fisher matrices is obtain with applying unbiased estimator of variation.\n");
    text_of_help.push_back("\t\t To change for biased version type: -u f, f - false, t - true (default set).\n");
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
    text_of_help.push_back("\t -u or --unbiased <bool>\n");
    text_of_help.push_back("\t\t Wariance estimator can be unbiased (true), or biased (false).\n");
    text_of_help.push_back("\t\t <bool> take argument: t (true), f (false).\n");
    text_of_help.push_back("\t\t Default set: t.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t   >>> Options for both All-Sky and Directed searches <<<\n");
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
    text_of_help.push_back("\t -i or --initial <min> <max> <step>\n");
    text_of_help.push_back("\t\t Initial time (for equations where ti: <-1/2 To, 1/2 To>. To - observation time).\n");
    text_of_help.push_back("\t\t <min> - minimum value of initial time;\n");
    text_of_help.push_back("\t\t default set: 0.5.\n");
    text_of_help.push_back("\t\t <max> - optional maximum value of initial time.\n");
    text_of_help.push_back("\t\t <step> - optional step value of initial time;\n");
    text_of_help.push_back("\t\t if not provided step will be set on step = max - min.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t -I or --Initial <min> <max> <step>\n");
    text_of_help.push_back("\t\t Initial time (for equations where ti: <0, To>. To - observation time).\n");
    text_of_help.push_back("\t\t <min> - minimum value of initial time;\n");
    text_of_help.push_back("\t\t default set: 0 (see default set for flag -i).\n");
    text_of_help.push_back("\t\t <max> - optional maximum value of initial time.\n");
    text_of_help.push_back("\t\t <step> - optional step value of initial time;\n");
    text_of_help.push_back("\t\t if not provided step will be set on step = max - min.\n");
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
    text_of_help.push_back("\t\t (for All-Sky searches; for Directed searches data length will be set to 344656).\n");
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
    text_of_help.push_back("\t\t *** Information stored in file 'dataS2.txt' are independent from \n");
    text_of_help.push_back("\t\t an ephemeris. Based on this information density of covering can \n");
    text_of_help.push_back("\t\t be obtained fast but in approximated way (generally speaking more \n");
    text_of_help.push_back("\t\t collected data allows to get more accurate results).\n");
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
    text_of_help.push_back("\t\t if any option provided, program will work with Michal Bejger's sets\n");
    text_of_help.push_back("\t\t (default set; result saved to 'grid.bin'); a - this same as: dfg.\n");
    text_of_help.push_back("\t\t Option can be joined (order is not important), e.g. df, dg, fg, dfg.\n");
    text_of_help.push_back("\t\t #### If argument of flag -nd is set to be bigger than zero, flag -p\n");
    text_of_help.push_back("\t\t #### can accept only 'd' argument (All-Sky).\n");
    text_of_help.push_back("\t -s or --spindown <int>\n");
    text_of_help.push_back("\t\t Number of spindowns to use.\n");
    text_of_help.push_back("\t\t <int> accept values: 1 or 2.\n");
    text_of_help.push_back("\t\t If number of spindowns is set to be 2, algorithm 's1' will be used.\n");
    text_of_help.push_back("\t\t (regardless to option used with flag '-a').\n");
    text_of_help.push_back("\t\t Default set: 1.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t -t or --type <type>\n");
    text_of_help.push_back("\t\t type of searches to choose.\n");
    text_of_help.push_back("\t\t <type> take option: a, d.\n");
    text_of_help.push_back("\t\t a - all sky, d - directed.\n");
    text_of_help.push_back("\t\t Default set: a.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("   >>> Flags with no arguments: <<<\n");
    text_of_help.push_back("\t -h or --help\n");
    text_of_help.push_back("\t\t Display this help.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t -ht or --helptxt\n");
    text_of_help.push_back("\t\t Print this help to text file.\n");
    text_of_help.push_back("\n");
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
    text_of_help.push_back("\n");
    text_of_help.push_back("\t ./gg -c 0.7 -i -10 >> result.txt\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t\t Result is appended to end of the file 'result.txt'.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t ./gg -c 0.7 0.95 -i -10 10 0.1 -n 21 -nd 344656 -p d\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t\t Generate density of covering (only) without using ephemeris (344656!=0).\n");
    text_of_help.push_back("\t\t GridsGenerator print density (-p d) obtained for:\n");
    text_of_help.push_back("\t\t - covariance from 0.7 to 0.95 with default step = 0.01,\n");
    text_of_help.push_back("\t\t (-i -10 10 0.1) initial time from -10 to 10 with step 0.1,\n");
    text_of_help.push_back("\t\t (-n 21)         Fourier length (--nfft) set to 2^21 (2097152),\n");
    text_of_help.push_back("\t\t (-nd 344656)    data length set on 344656.\n");
    text_of_help.push_back("\t\t Ephemeris are located at directory with GridsGenerator.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t ./gg -d 001/ -c 0.7 0.95 -i -10 10 0.1 -n 21 -p dfg\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t\t GG prints: (d) densities, (f) fisher matrices, (g) grids.\n");
    text_of_help.push_back("\t\t (-d 001/) Ephemeris are located at local directory '001/'.\n");
    text_of_help.push_back("\t ./gg -d 001/ -c 0.7 0.95 0.001 -i -10 10 0.1 -n 21 -p dfg\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t\t Covariance step is set to 0.001,\n");
    text_of_help.push_back("\t\t Nfft length is set to 2^21 (2097152).\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t ./gg -d 001/ -c 0.7 0.95 0.001 -i -10 10 0.1 -a a -na 40 -nr 30 -cv f\n");
    text_of_help.push_back("\t      -ch t -p dfg\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t\t Nfft length is set to default value 2^20 (1048576).\n");
    text_of_help.push_back("\t\t (-a a)   GridsGenerator decide which algorithm should be used. \n");
    text_of_help.push_back("\t\t          If GG decide to use algorithm 's2' next two flags\n");
    text_of_help.push_back("\t\t          will be used: \n");
    text_of_help.push_back("\t\t (-na 40) Number of loops in root finding algorithm is set to 40.\n");
    text_of_help.push_back("\t\t (-nr 30) Number of loops in covering radius finding algorithm is set to 30.\n");
    text_of_help.push_back("\t\t (-cv f)  Grids are not converted to hyper-ellipsoid space.\n");
    text_of_help.push_back("\t\t (-ch t)  All values are chopped to zero if are smaller than |10^-12|.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t ./gg -c 0.2 -a s1 -s 1 -n 19 -cv t -p d \n");
    text_of_help.push_back("\t\t (-c 0.2) Covariance is set to 0.2,\n");
    text_of_help.push_back("\t\t (-a s1)  algorithm 's1' will be used,\n");
    text_of_help.push_back("\t\t (-s 1)   first spindown will be used,\n");
    text_of_help.push_back("\t\t (-n 19)  Fourier length (--nfft) set to 2^19 (524288).\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t ./gg -c 0.2 -a s1 -s 2 -n 19 -cv t -p d \n");
    text_of_help.push_back("\t\t (-s 2) Number of spindown to use (first and second spindowns will be used).\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t ./gg -c 0.2 -a s2 -s 2 -n 19 -cv t -p d \n");
    text_of_help.push_back("\t\t (-s 2)   Number of spindown to use.\n");
    text_of_help.push_back("\t\t (-a s2)  ! Algorithm 's1' will be used ('s2' is implement only for\n");
    text_of_help.push_back("\t\t          first spindown (-s 1) ).\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t ./gg -h -v -aa\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t\t Display help, version of GridsGenerator and show info about author.\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t\t # Note that GridsGenerator provide result more effective (faster) in way:\n");
    text_of_help.push_back("\t ./gg -c 0.70 0.73 -d 001/\n");
    text_of_help.push_back("\t\t # than when is running several times for the same ephemeris:\n");
    text_of_help.push_back("\t ./gg -c 0.70 -d 001/\n");
    text_of_help.push_back("\t ./gg -c 0.71 -d 001/\n");
    text_of_help.push_back("\t ./gg -c 0.72 -d 001/\n");
    text_of_help.push_back("\t ./gg -c 0.73 -d 001/\n");
    text_of_help.push_back("\t\t # That clue apply to both: covariance (-c)/minimal match (-m) and\n");
    text_of_help.push_back("\t\t # initial time (-i).\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t EXAMPLES for All-Sky searches (more than one detector).\n");
    text_of_help.push_back("\t\t To get result for more than one detector (with applying\n");
    text_of_help.push_back("\t\t sum of reduced fisher matrices):\n");
    text_of_help.push_back("\t\t ./gg -d /home/ligo/ 001 H1 L1 0067\n");
    text_of_help.push_back("\t\t ./gg -d /home/ligo/ 001 H1 L1 V1 0067.\n");
    text_of_help.push_back("\t\t To read data from file different than default 'xdats', provide appropriate\n");
    text_of_help.push_back("\t\t name of the file on the last position:\n");
    text_of_help.push_back("\t\t ./gg -d /home/ligo/ 001 H1 L1 0067 xdatsc_001_0067.bin\n");
    text_of_help.push_back("\t\t or more simpler form:\n");
    text_of_help.push_back("\t\t ./gg -d /home/ligo/ 001 H1 L1 0067 xdatsc\n");
    text_of_help.push_back("\t\t Data file name need to have prefix 'xdats'.\n");
    text_of_help.push_back("\t\t Accepted detectors names: H1, L1, V1.\n");
    text_of_help.push_back("\t\t Sum of Fisher matrices is obtain with applying unbiased estimator of variation.\n");
    text_of_help.push_back("\t\t To change for biased version type: -u f, f - false, t - true (default set).\n");
    text_of_help.push_back("\t\t Both notation:\n");
    text_of_help.push_back("\t ./gg -c 0.75 -a s1 -s 1 -n 19 -d ../O3_6d_test/001/V1/\n");
    text_of_help.push_back("\t\t and\n");
    text_of_help.push_back("\t ./gg -c 0.75 -a s1 -s 1 -n 19 -d ../O3_6d_test/ 001 V1\n");
    text_of_help.push_back("\t\t are equivalent (prints this same result).\n");
    text_of_help.push_back("\t ./gg -c 0.75 -a s1 -s 1 -n 19 -d ../O3_6d_test/ 001 H1 L1 0067 xdatsc -u f\n");
    text_of_help.push_back("\t\t Output will be saved in folder: ../O3_6d_test/001/grids/ (\"grids\" folder\n");
    text_of_help.push_back("\t\t have to exist). Name of output file: grid_001_0067_H1L1c.bin\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t EXAMPLES for Directed searches\n");
    text_of_help.push_back("\t ./gg -c 0.75 -i 0.0 -t d -nd 300000 -a s1 -s 1 -cv t\n");
    text_of_help.push_back("\t ./gg -c 0.75 -i 0.5 -t d  -a s1 -s 1 -p d\n");
    text_of_help.push_back("\t ./gg -c 0.75 -i 0.5 -t d  -a s1 -s 2 -cv t\n");
    text_of_help.push_back("\t ./gg -c 0.75 -i 0.5 -t d  -a s1 -s 3 -cv f\n");
    text_of_help.push_back("\t ./gg -c 0.75 0.99 -i 0.5 -t d -a s1 -s 2 -p d\n");
    text_of_help.push_back("\t ./gg -c 0.01 0.999 0.001 -i 0.5 -t d -a s1 -s 3 -n 19 -p d\n");
    text_of_help.push_back("\t ./gg -c 0.01 0.999 0.001 -i 0.5 -t d -a s2 -s 1 -n 19 -p d\n");
    text_of_help.push_back("\t ./gg -c 0.01 0.999 0.001 -i 0.5 -t d -a s2 -s 1 -n 19 -nd 689312 -p d\n");
    text_of_help.push_back("\n");
    text_of_help.push_back("\t **********************************************************\n");
    text_of_help.push_back("\t * POLGRAW Group                                          *\n");
    text_of_help.push_back("\t * University of Bialystok, Faculty of Physics            *\n");
    text_of_help.push_back("\t * Developer: Andrzej Pisarski <a.pisarski at uwb.edu.pl> *\n");
    text_of_help.push_back("\t * Tutor: Piotr Jaranowski <p.jaranowski at uwb.edu.pl>   *\n");
    text_of_help.push_back("\t **********************************************************\n");

    for(int i=0; i<21; i++)
        text_of_help.push_back("\n");


    return text_of_help;
}

std::vector<std::string> Manual::init_m_version(std::string ver="") const
{
    std::vector<std::string> text_of_version;
    text_of_version.reserve(7);

    text_of_version.push_back("\t *******************************************");
    text_of_version.push_back(ver);
    text_of_version.push_back("\t * Copyright 2015 - 2020 Andrzej Pisarski. *");
    text_of_version.push_back("\t * License of GridsGenerator: CC-BY-NC-SA. *");
    text_of_version.push_back("\t * Additional credits to:                  *");
    text_of_version.push_back("\t * Piotr Jaranowski                        *");
    text_of_version.push_back("\t *******************************************");

    return text_of_version;
}
