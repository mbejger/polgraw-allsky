///////////////////////////////////////
// Name:        DensityS2.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#include "DensityS2.h"
#include "num.h"
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <iostream>

DensityS2::DensityS2(): m_data(dataS2_load())
{

}

std::vector<DataS2> DensityS2::dataS2_load(std::string path) const
{
    double omega0p, density, alpha, radius;
    int n_alpha, n_radius;

	std::ifstream in(path.c_str());
	if(!in.is_open())
	{
        std::string error = "Can't open file: "+path;
        throw std::runtime_error(error);
	}

    std::vector<DataS2> t_data;
    t_data.reserve(100);
    std::string line;
    in >> std::skipws;
    std::getline(in,line);
    while(in >> omega0p >> density >> alpha  >> radius >> n_alpha >> n_radius){
        t_data.push_back(DataS2(omega0p, density, alpha,  radius,n_alpha, n_radius));
    }
    in.close();

    //for(unsigned int i=0; i<t_data.size(); i++)
    //    std::cout<< i << "\t" << t_data[i].m_omega0p <<"\n";

    return t_data;
}

std::vector<DataS2> DensityS2::get_data() const
{
    return m_data;
}

double DensityS2::density_approx(double deltaOmegaZeroPrim) const
{
    double x = deltaOmegaZeroPrim;
    if(x>=8.-num::epsilon())
    {
        std::string error="Grid S2 do not exist for this data sets.";
        throw std::domain_error(error);
    }
    unsigned int data_size = m_data.size();
    //std::cout << data_size << "\n";
    if(x<m_data[0].m_omega0p || x>m_data[data_size-1].m_omega0p)
    {
        std::string sx = std::to_string(x);
        std::string error="Density can not be obtain: do not enough data.\n";
        error+="DelOmZeroP=";
        error+=sx;
        error+=" is outside of range of available data in dataS2.txt.\n";
        throw std::domain_error(error);
    }
    unsigned int i=0;
    while(i<data_size && m_data[i].m_omega0p<=x)
        i++;

    double x1=m_data[i-1].m_omega0p, x2=m_data[i].m_omega0p;
    double y1=m_data[i-1].m_density, y2=m_data[i].m_density;

    return (x-x1)*(y2-y1)/(x2-x1)+y1;
}

double DensityS2::density_approx(double c0, unsigned int nfft, unsigned int data_length) const
{
    double x = num::delta_omega_zero_prim(c0, nfft, data_length);
    return density_approx(x);
}


