///////////////////////////////////////
// Name:        ReadData.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:
// Created:     03/06/2020
///////////////////////////////////////

#include "ReadData.h"
#include <stdexcept>
#include <fstream>

ReadData::ReadData(const std::string& path)//&
{

    std::ifstream file(path.c_str(), std::ios::binary | std::ios::in);  //False -> 0
    if(file.fail())
    {
        std::string error="Can not find: " + path;
        throw std::runtime_error(error);
    }

    if(file.good())
    {
        file.seekg(0, file.end);
        unsigned int length = file.tellg();
        file.seekg(0, file.beg);
        length/=sizeof(double);
        m_length=length;
        m_data.resize(m_length);

        for(unsigned int i=0; i<m_length; i++)
        {
            file.read( reinterpret_cast<char*>( &m_data[i] ), sizeof(m_data[0]) );

        }
    }
    file.close();
}

const std::vector<double>& ReadData::get_data() const
{
    return m_data;
}

unsigned int ReadData::get_length() const
{
    return m_length;
}
