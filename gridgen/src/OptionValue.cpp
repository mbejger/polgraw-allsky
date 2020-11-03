///////////////////////////////////////
// Name:        OptionValue.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
///////////////////////////////////////

#include "OptionValue.h"

OptionValue::OptionValue(): m_k(-1), m_i(0), m_d(-1.0), m_s(""){}

OptionValue::OptionValue(int k, int i=0, double d=0.0, std::string s=""): m_k(k), m_i(i), m_d(d), m_s(s){}
