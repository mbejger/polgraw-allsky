///////////////////////////////////////
// Name:        ReadEphemeris.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
// Modification:07/06/2020 A.Pisarski
///////////////////////////////////////

#include "ReadEphemeris.h"
#include <stdexcept>
#include <fstream>
#include <cmath>

ReadEphemeris::ReadEphemeris(const std::vector<std::string> paths)//&
{
    if(paths.size()%3!=0)
    {
        //std::cout << "paths.size()=" << paths.size() << std::endl;
        std::string error="Missing at least one of paths files!";
        throw std::runtime_error(error);
    }

    std::string fullpathSSB = paths[0];
    std::string fullpathDet = paths[1];
    std::string fullpathDetSSB = paths[2];

    std::ifstream is1((fullpathSSB).c_str(), std::ios::binary);  //False -> 0
    m_flag.resize(3, 1);
    if(!is1.is_open())
    {
        m_flag[0]=0;    // fail
        std::string error="Can not find: "+fullpathSSB;
        throw std::runtime_error(error);
    }

    std::ifstream is2((fullpathDet).c_str(), std::ios::binary);
    if(!is2.is_open())
    {
        m_flag[1]=0;    // fail
        std::string error="Can not find: "+fullpathDet;
        throw std::runtime_error(error);
    }

    std::ifstream is3((fullpathDetSSB).c_str(), std::ios::binary);
    if(!is3.is_open())
    {
        m_flag[2]=0;    // fail
        std::string error="Can not find: "+fullpathDetSSB;
        throw std::runtime_error(error);
    }

    if(m_flag[0]&&m_flag[1]&&m_flag[2])     /// probably it is safe to remove 'If' statement
    {
        m_DetSSB.resize(3);
        unsigned int i=0;

        // to read last 3 entries
        while(is3.read( reinterpret_cast<char*>( &m_DetSSB[i++] ), sizeof(m_DetSSB[0]) ))
        {
            if(i==3)    i=0;
        }
        //std::cout<<m_DetSSB[0]<<std::endl;
        //std::cout<<m_DetSSB[1]<<std::endl;
        //std::cout<<m_DetSSB[2]<<std::endl;

        m_epsm=m_DetSSB[1];
        m_ce=cos(m_epsm);

        std::vector<double> rSSB(3);
        std::vector<double> rDet(3);

        is1.seekg(0, is1.end);
        unsigned int length = is1.tellg();
        is1.seekg(0, is1.beg);
        length/=(3*sizeof(rSSB[0]));
        m_length=length;
        ///std::cout<<"length="<<length<<std::endl;

        // reserve - "method only allocates memory, but leaves it uninitialized.
        // It only affects capacity(), but size will be unchanged."
        // http://stackoverflow.com/questions/7397768
        m_ephemeris1.reserve(m_length+4);
        m_ephemeris2.reserve(m_length+4);

        double ws=1.0/m_length; //scale factor

        unsigned int j=0;
        while(is1&&is2&&j<m_length) // read to the end of the stream s1 or s2 but not no more than 'n' lines
        {
            for(unsigned int i=0; i<3; i++)
            {
                is1.read( reinterpret_cast<char*>( &rSSB[i] ), sizeof(rSSB[0]) );
                is2.read( reinterpret_cast<char*>( &rDet[i] ), sizeof(rDet[0]) );
            }

            m_ephemeris1.push_back(ws*(rSSB[1]/m_ce+rDet[1]*m_ce));   // rSSB[1]/ce ??
            m_ephemeris2.push_back(ws*(rSSB[0]+rDet[0]));

            j++;
        }
    }
    is1.close();
    is2.close();
    is3.close();
}

const std::vector<double>& ReadEphemeris::get_ephemeris1() const
{
    return m_ephemeris1;
}

const std::vector<double>& ReadEphemeris::get_ephemeris2() const
{
    return m_ephemeris2;
}

unsigned int ReadEphemeris::get_ephemeris_length() const
{
    return m_length;
}
