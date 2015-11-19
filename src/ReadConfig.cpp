///////////////////////////////////////
// Name:        ReadConfig.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
///////////////////////////////////////

#include "ReadConfig.h"
#include <fstream> //std::ifstream
#include <iostream>
#include <algorithm>
#include <stdexcept>

//#include <iostream>

// vector of available options (order is important)
const std::vector<std::string> ReadConfig::m_options{
    "CovarianceMin=", "CovarianceMax=", "CovarianceStep=",
    "InitialTimeMin=", "InitialTimeMax=", "InitialTimeStep=",
    "NfftExpLength=",
    "SetDataLength=",
    "NAlpha=",                              //  8
    "NRadius=",
    "PathSSB=", "PathDet=", "PathDetSSB=",
    "PathSave=",
    "FilePatternGrid=",
    "FilePatternFM=",
    "FilePatternDoC=",                      // 16
    "SaveGrid=",
    "SaveFisherMatrix=",
    "SaveDensityOfCovering=",
    "Chop=",                                // 20
    "ChooseMethod=",
    "DataS2Use=",
    "DataS2SaveBest=",
    "Convert=",                             // 24
    "Quiet=",
    "SayHello="                             // 26
};

ReadConfig::ReadConfig(): m_values(read_config_file())
{
    //ctor
}

ReadConfig::ReadConfig(const std::string& path): m_values(read_config_file(path))
{

}

ReadConfig::ReadConfig(const std::vector<OptionValue> opt_values): m_values(opt_values)
{
    //ctor
}

ReadConfig::~ReadConfig()
{
    //destructor
}

// returns vector of values for options founded in configuration file
// (if some options are missing or badly define options will be set on default values)
std::vector<OptionValue> ReadConfig::read_config_file(const std::string& path) const
{
    std::vector<std::string> veclin;
    std::ifstream in(path.c_str());
    if(!in.is_open())
    {
        //std::cout<<"Can't open "<<path<<std::endl;
        //flag = 0;
        std::string error = "Can't open file: "+path;
        throw std::runtime_error(error);
    }
    else
    {
        std::string line;
        ///std::cout<<"File: '"<<path<<"' opened successful."<<std::endl;
        typedef std::string::const_iterator iter;
        int mark1=0;
        while(getline(in, line))
        {
            iter b = line.begin(), e = line.end();
            if(set_line(b, e))
                mark1=1;

            if(mark1==1)
                veclin.push_back(line);
        }
        in.close();
        //flag = 1;

    }

    std::vector<OptionValue> tt=set_command(veclin);

    return tt;
}

// To read configuration beginning from line containing string '[GridsGenerator.Settings]'
bool ReadConfig::set_line(std::string::const_iterator b, std::string::const_iterator e) const
{
    static const std::string sli = "[GridsGenerator.Settings]";
    typedef std::string::const_iterator iter;

    // i marks where 'sli' was found
    iter i = b;
    i = search(i, e, sli.begin(), sli.end());

    return i != e;
}

// return vector of values associated with option vector: 'm_options',
// returns not only values for options founded in configuration file
// (if some are missing or badly define they will be set on default values)
std::vector<OptionValue> ReadConfig::set_command(const std::vector<std::string>& config_lines) const
{
    typedef std::string::const_iterator iter;
    std::vector<OptionValue> tempV;
    tempV.resize(27);

    unsigned int k;
    for(unsigned int line=0; line<config_lines.size(); line++)
    {
        iter b = config_lines[line].begin(), e = config_lines[line].end(), i;
        for(k=0; k<m_options.size(); k++)
        {
            i = b;
            i = search(i, e, m_options[k].begin(), m_options[k].end());
            if(i != e)
            {
                i+=m_options[k].size();
                std::string temp(i, e);
                    //std::cout << temp << std::endl;
                ///double v = std::atof(temp.c_str());
                if(k<10)
                {
                    if(k<6)
                    {
                        double vd = std::atof(temp.c_str());
                        tempV[k] = OptionValue(k, -1, vd, temp);
                    }
                    else{
                        int vi = std::atoi(temp.c_str());
                        tempV[k] = OptionValue(k, vi, -1.0, temp);
                    }
                }
                else
                    tempV[k] = OptionValue(k, -1, -1.0, temp);

            }
        }
    }

    // Checking if values are well defined
    if(tempV[0].m_d < 0.0 || tempV[0].m_d >= 1.0)
        tempV[0].m_d = 0.75;            // CovarianceMin

    double co_diff = tempV[1].m_d - tempV[0].m_d;
    if(tempV[1].m_d < 0.0 || tempV[1].m_d >= 1.0 ||  co_diff < 0)
        tempV[1].m_d = tempV[0].m_d;    // CovarianceMax

    co_diff = tempV[1].m_d - tempV[0].m_d;
    if(tempV[2].m_d <= 0.0)             // CovarianceStep
    {
        if(co_diff>0)
            tempV[2].m_d =co_diff;
        else
            tempV[2].m_d =1.0;
    }

    double it_diff = tempV[4].m_d - tempV[3].m_d;
    if(it_diff < 0.0)
        tempV[4].m_d = tempV[3].m_d;    // InitialTimeMax

    it_diff = tempV[4].m_d - tempV[3].m_d;
    if(tempV[5].m_d <= 0.0)             // InitialTimeStep
    {
        if(it_diff>0)
            tempV[5].m_d =it_diff;
        else
            tempV[5].m_d =1.0;
    }

    if(tempV[6].m_i < 1)
        tempV[6].m_i = 20;              // NfftExpLength

    if(tempV[7].m_i < 0)
        tempV[7].m_i = 0;               // zero means SetDataLength equal to ephemeris length

    if(tempV[8].m_i <= 0)
        tempV[8].m_i = 18;              // NAlpha

    if(tempV[9].m_i <= 0)
        tempV[9].m_i = 18;              // NRadius

    for(int i=10; i<13; i++)            // Source files
        if(tempV[i].m_s.size() == 0)
        {
            std::string error = "Can't find " + tempV[i].m_s;
            throw std::runtime_error(error);
        }

    if(tempV[13].m_s.size() == 0)       // Output directory
    {
        std::string error = "Can't open " + tempV[13].m_s;
        throw std::runtime_error(error);
    }

    if(tempV[14].m_s.size() == 0)       // File name pattern (grid)
        tempV[14].m_s = "grid_Co%C_IT%I_N%N.txt";

    if(tempV[15].m_s.size() == 0)       // File name pattern (Fisher matrix)
        tempV[15].m_s = "fisher_IT%I.txt";

    if(tempV[16].m_s.size() == 0)
        tempV[16].m_s = "DoC_DL%D_N%N.txt"; // File name pattern (density of covering)

    if(tempV[17].m_s != "True" && tempV[17].m_s != "False")
        tempV[17].m_s = "True";

    for(int i=18; i<21; i++)            // Choose if obtain grid, Fisher matrix or density of covering
        if(tempV[i].m_s != "True" && tempV[i].m_s != "False")   // or chop result
            tempV[i].m_s = "False";

    if(tempV[21].m_s != "S1" && tempV[21].m_s != "S2"
       && tempV[21].m_s != "Automatic") // Choose algorithm to obtain grid or Fisher matrix or density
        tempV[21].m_s = "Automatic";

    for(int i=22; i<25; i++)            // (22) DataS2Use, (23) DataS2SaveBest
        if(tempV[i].m_s != "True" && tempV[i].m_s != "False")
            tempV[i].m_s = "True";

    if(tempV[24].m_s != "True" && tempV[24].m_s != "False") // (24) Convert
        tempV[24].m_s = "True";

    for(int i=25; i<27; i++)            // (25) "Quiet=", (26) "SayHello="
        if(tempV[i].m_s != "True" && tempV[25].m_s != "False") // , (25) Future option, not implement yet
            tempV[i].m_s = "False";
/*
    for(unsigned int i =0; i<6; i++)
        std::cout<< "OptionValue[" << i << "]=" << tempV[i].m_d << "\n";

    for(unsigned int i =6; i<10; i++)
        std::cout<< "OptionValue[" << i << "]=" << tempV[i].m_i << "\n";

    for(unsigned int i =10; i<tempV.size(); i++)
        std::cout<< "OptionValue[" << i << "]=" << tempV[i].m_s << "\n";
*/
    return tempV;
}

std::vector<std::string> ReadConfig::get_options() const
{
    return m_options;
}

std::vector<OptionValue> ReadConfig::get_values() const
{
    return m_values;
}
