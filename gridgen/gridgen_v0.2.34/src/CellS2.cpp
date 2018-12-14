///////////////////////////////////////
// Name:        CellS2.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#include "CellS2.h"

CellS2::CellS2(): m_alpha(0.), m_radius(0.), m_hole(4)
{
    //ctor
}

CellS2::CellS2(double alpha, double radius, std::vector<double> hole)
: m_alpha(alpha), m_radius(radius), m_hole(hole)
{
    //ctor
}
/*
CellS2::CellS2(const CellS2& C): m_alpha(C.m_alpha), m_radius(C.m_radius), m_hole(C.get_hole())
{

}

void CellS2::operator=(const CellS2& C)
{
   m_alpha=C.m_alpha;
   m_radius=C.m_radius;
   m_hole=C.get_hole();
}
*/

std::vector<double> CellS2::get_hole() const
{
    return m_hole;
}

bool compare_CellS2_alpha(const CellS2& sa, const CellS2& sb)
{
    return sa.m_alpha < sb.m_alpha;
}


bool compare_CellS2_radius(const CellS2& sa, const CellS2& sb)
{
    return sa.m_radius < sb.m_radius;
}


bool compare_CellS2_alpha_decrease(const CellS2& sa, const CellS2& sb)
{
    return sa.m_alpha > sb.m_alpha;
}


bool compare_CellS2_radius_decrease(const CellS2& sa, const CellS2& sb)
{
    return sa.m_radius > sb.m_radius;
}



