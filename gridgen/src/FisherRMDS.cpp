///////////////////////////////////////
// Name:        FisherRMDS.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     30/07/2019
// Modification:13/10/2019
///////////////////////////////////////

#include "FisherRMDS.h"
#include "num.h"
#include <stdexcept>
#include <cmath>

FisherRMDS::FisherRMDS(): m_spindown(-1)
{
    //m_prefisher.resize(0);
    //ctor
}
/// Ephemeris in directed searches are not used
FisherRMDS::FisherRMDS(int spindown_): m_spindown(spindown_)
{
    //ctor
}

std::vector<double> FisherRMDS::postrmf(double xi) const
{
    return fisherMDS(xi, m_spindown);
}

unsigned int FisherRMDS::dim() const
{
    return m_spindown + 1;
}


unsigned int FisherRMDS::get_spindown() const
{
    return m_spindown;
}


/// Below functions are designated for directed searches and can not be used
/// directly (without some changes) in all-sky searches (dimension of the
/// grid is equal to: spindown + 1)

/// Auxiliary function - simple way to add Fisher matrix elements (in directed searches case).
std::vector<double> FisherRMDS::fisher_elemetsDS(double xi, unsigned int spindown) const
{
    double xi2 = pow(xi, 2), xi4 = pow(xi, 4), xi6 = pow(xi, 6);
    double xi3 = pow(xi, 3), xi5 = pow(xi, 5);

    std::vector<double> temp;
    temp.push_back(1/12.0);               ///[0][0], -std=c++11

    if(spindown >= 1){
        temp.push_back(xi/6.0);             ///[1][0]
        temp.push_back(1/180.0 + xi2/3.0);  ///[1][1]
    }

    if(spindown >= 2){
        temp.push_back(1/80.0 + xi2/4.0);                ///[2][0]
        temp.push_back((xi + 12*xi3)/24.0);              ///[2][1]
        temp.push_back(1/448.0 + xi2/8.0 + 3*xi4/4.0);   ///[2][2]
    }

    if(spindown == num::max_spindownDS){
        temp.push_back(xi*(3 + 20*xi2)/60.0);                   ///[3][0]
        temp.push_back(1/840.0 + 2*(xi2/15.0 + xi4/3.0) );      ///[3][1]
        temp.push_back(xi/80.0 + 3*xi3/10.0 + xi5 );               ///[3][2]
        temp.push_back(1/3600.0 + xi2/20.0 + 3*xi4/5.0 + 4*xi6/3.0 );    ///[3][3]
     ///temp.push_back((1/60.0 + xi2*(3.0 + xi2*(36.0+80*xi2)))/60.0 ); ///[3][3]
    }

    if(spindown > num::max_spindownDS){
        std::string string_err = "This version of Grid Generator do not support spindown bigger than " + std::to_string(num::max_spindownDS) + ".\n";
        throw std::domain_error(string_err);
    }

    return temp;
}

std::vector<double> FisherRMDS::fisherMDS(double xi, int spindown) const
{
    double dim = spindown + 1;
    std::vector<double> ta(dim*dim, 0.0);

    std::vector<double> fe = fisher_elemetsDS(xi, spindown);

    for(int i=0; i<dim; i++){
        int index = ((1+i)*i)/2;///arithmetic sum
        int line = dim * i;
        for(int j=index; j<index+i+1; j++){
            ta[line++]=fe[j];
        }
    }

    for(int i=0; i<dim; i++){
        for(int j=0; j<dim-i-1; j++){
            ta[i*(dim+1)+j+1]=ta[dim*(i+j+1)+i];
        }
    }
/*
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            std::cout<<std::setw(10)<< ta[(dim*i)+j]<<", ";
        }
        std::cout<<std::endl;
    }
*/
    return ta;
}

