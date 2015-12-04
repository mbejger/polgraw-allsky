///////////////////////////////////////
// Name:        GridS1.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#include "GridS1.h"
#include "num.h"
#include <stdexcept>
#include <cmath>
#include <string>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

GridS1::GridS1(const FisherRM *const fm=nullptr): DensityS1(), m_fm(fm)
{
    if(m_fm==nullptr)
    {
        std::string error="Fisher Matrix is missing.\n";
        throw std::runtime_error(error);
    }

    if(m_fm->dim()!=4)
    {
        std::string error="Dimension of Fisher Reduced Matrix must be = 4.\n";
        throw std::domain_error(error);
    }
}

std::vector<double> GridS1::grid(double c0=0.75, double xi=0., unsigned int nfft=524288) const
{
    std::vector<double> ws1(16, 0.0);
    unsigned int data_length = m_fm->get_ephemeris_length();

    double r_scale=sqrt(1.0-c0);
    std::vector<double> q0 = DensityS1::grid_prim(c0, nfft, data_length);

    //ws1[0]=DensityS1::density(q0);
    std::vector<double> temp;
    temp = num::mult(q0, num::inverse( num::cholesky( m_fm->postrmf(xi), 4), 4), 4, 4, 4);
    for(int i=0; i<16; i++)
        ws1[i]=r_scale*temp[i]; //s1[i+1]=

    return ws1;
}

// convert vector from hyper-sphere space to hyper-ellipsoid space
std::vector<double> GridS1::convert(double c0, double xi, const std::vector<double>& grid_prim)
{
    std::vector<double> ws1(16);
    double r_scale=sqrt(1.0-c0);

    std::vector<double> temp;
    temp = num::mult(grid_prim, num::inverse( num::cholesky( m_fm->postrmf(xi), 4), 4), 4, 4, 4);
    for(int i=0; i<16; i++)
        ws1[i]=r_scale*temp[i]; //ws1[i]=

    return ws1;
}

double GridS1::density(double c0=0.75,  unsigned int nfft=524288) const
{
    unsigned int data_length=m_fm->get_ephemeris_length();
    return DensityS1::density(c0, nfft, data_length);
}

double GridS1::density(double c0,  unsigned int nfft, unsigned int data_length) const
{
    return DensityS1::density(c0, nfft, data_length);
}

