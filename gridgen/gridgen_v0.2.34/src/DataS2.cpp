///////////////////////////////////////
// Name:        DataS2.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     16/06/2015
///////////////////////////////////////

#include "DataS2.h"

DataS2::DataS2(): m_omega0p(0.), m_density(0.), m_alpha(0.), m_radius(0.), m_NAlpha(0), m_NRadius(0){}

DataS2::DataS2(double omega, double density, double alpha, double radius, int n_alpha, int n_radius):
    m_omega0p(omega), m_density(density), m_alpha(alpha), m_radius(radius),
    m_NAlpha(n_alpha), m_NRadius(n_radius)
{
}

int DataS2::get_nalpha() const
{
    return m_NAlpha;
}

int DataS2::get_nradius() const
{
    return m_NRadius;
}

bool compare_DataS2_omega0p(const DataS2& da, const DataS2& db)
{
    return da.m_omega0p < db.m_omega0p;
}


bool compare_DataS2_density(const DataS2& da, const DataS2& db)
{
    return da.m_density < db.m_density;
}

bool compare_DataS2_alpha(const DataS2& da, const DataS2& db)
{
    return da.m_alpha < db.m_alpha;
}

bool compare_DataS2_radius(const DataS2& da, const DataS2& db)
{
    return da.m_radius  < db.m_radius;
}

///decrease
bool compare_DataS2_omega0p_decrease(const DataS2& da, const DataS2& db)
{
    return da.m_omega0p > db.m_omega0p;
}


bool compare_DataS2_density_decrease(const DataS2& da, const DataS2& db)
{
    return da.m_density > db.m_density;
}

bool compare_DataS2_alpha_decrease(const DataS2& da, const DataS2& db)
{
    return da.m_alpha > db.m_alpha;
}

bool compare_DataS2_radius_decrease(const DataS2& da, const DataS2& db)
{
    return da.m_radius  > db.m_radius;
}

