///////////////////////////////////////
// Name:        stringmanip.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
///////////////////////////////////////

#include "stringmanip.h"

#include <algorithm>
//#include <iostream>

void stringmanip::sreplace(std::string& text2modify, const std::string& text2find, std::string text_new)
{
    if( text2modify.size() >= text2find.size()) // text2find obviously can not be longer that modified text
    {
        size_t tbeg = 0;//, tend = text2modify.size();
        while( tbeg < text2modify.size() )
        {
            tbeg = text2modify.find(text2find, tbeg);

            if( tbeg != std::string::npos)
                text2modify.replace(tbeg, text2find.size(), text_new);
        }
    }
}
