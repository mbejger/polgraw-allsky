///////////////////////////////////////
// Name:        ReadData.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:
// Created:     03/06/2020
///////////////////////////////////////

#ifndef READDATA_H
#define READDATA_H

#include <vector>
#include <string>

class ReadData
{
    public:
        ReadData(const std::string& path);

        const std::vector<double>& get_data() const;
        unsigned int get_length() const;

    private:
        std::vector<double> m_data;
        unsigned int m_length;                    // number of elements in ephemeris1, ephemeris2
};


#endif // READDATA_H
