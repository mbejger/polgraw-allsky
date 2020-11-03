///////////////////////////////////////
// Name:        OptionValue.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
// Modification:28/08/2017 A.P.
///////////////////////////////////////

#ifndef OPTIONVALUE_H
#define OPTIONVALUE_H
//#define _GLIBCXX_USE_CXX11_ABI 0

#include <string>

class OptionValue
{
    public:
        OptionValue();
        OptionValue(int, int, double, std::string);

        int m_k;
        int m_i;
        double m_d;
        std::string m_s;
};

#endif // OPTIONVALUE_H


