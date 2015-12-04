///////////////////////////////////////
// Name:        DensityS1.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#include "DensityS1.h"
#include "Node.h"
#include "num.h"
#include <string>
#include <stdexcept>    // domain_error
#include <cmath>        // sqrt
#include <algorithm>    // sort

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

DensityS1::DensityS1()
{
    //ctor
}

std::vector<double> DensityS1::grid_prim(double c0, unsigned int nfft, unsigned int data_length) const
{
    if(data_length<=0)
    {
        std::string error="Data length can not be <= 0.\n";
        throw std::domain_error(error);
    }

    std::vector<double> q0;
    //unsigned int inn=pow(2, nfft);
    if(data_length <= nfft)
    {
        double dwp=num::delta_omega_zero_prim(c0, nfft, data_length);
        double dlqq=0.0;
        int wi = 0.0, hd = static_cast<int>(ceil(dwp/sqrt(2.0)))+3;
        std::vector<double> a4(16, 0.0);
        std::vector<int> wek;
        std::vector<double> qq(4, 0.0);

        std::vector<Node> vi_temp;

        a4[0]=sqrt(2.0);

        a4[4]=a4[8]=a4[12]=1.0/(2.0 * sqrt(2));
        a4[5]=sqrt(7.5)/2.0;   // sqrt(15.0/2.0)

        a4[9]=a4[13]=-sqrt(5.0/6.0)/2.0;
        a4[10]=sqrt(5.0/3.0);

        a4[14]=-a4[10]/2.0;
        a4[15]=sqrt(5)/2.0;

        for(int i=0; i<=hd; i++)
            for(int j=0; j<=hd; j++)
                for(int k=0; k<=hd; k++)
                    for(int l=0; l<=hd; l++)
                    {
                        if( i==1 || j==1 || k==1 || l==1 )
                        {
                            std::vector<double> vi(16, 0.0);
                            for(int q=0; q<4; q++)
                            {
                                vi[q]+=i*a4[q];
                                vi[4*1+q]+=j*a4[4*1+q];
                                vi[4*2+q]+=k*a4[4*2+q];
                                vi[4*3+q]+=l*a4[4*3+q];
                            }

                            std::vector<double> vitemp(4,0.0);
                            for(int p=0; p<4; p++)
                                for(int q=0; q<4; q++)
                                    vitemp[p]+=vi[4*q+p];

                            double dtemp=num::length(vitemp);
                            if( dtemp >= dwp )//= dwp
                            {
                                std::vector<int> vtemp(4,0);//={i,j,k,l};
                                vtemp[0]=i;
                                vtemp[1]=j;
                                vtemp[2]=k;
                                vtemp[3]=l;

                                Node vtp=Node(dtemp-dwp, vtemp);
                                vi_temp.push_back(vtp);
                            }
                        }
                    }

        std::sort(vi_temp.begin(), vi_temp.end(), compare_Node_distance);

        //vi_temp.size()
        /*
        for(int i=0; i<120; i++)
        {
            std::vector<int> opt_values;
            opt_values=vi_temp[i].get_coef();

            std::cout << vi_temp[i].m_distance << " | ";
            for(int j=0; j<4; j++)
                std::cout << opt_values[j] << " | ";

            std::cout << std::endl;
        }
        std::cout << " ============= " << std::endl;
        */

        wek=vi_temp[0].get_coef();

        std::vector<double> vi(16, 0.0);
        for(int q=0; q<4; q++)
        {
            vi[q]+=wek[0]*a4[q];
            vi[4*1+q]+=wek[1]*a4[4*1+q];
            vi[4*2+q]+=wek[2]*a4[4*2+q];
            vi[4*3+q]+=wek[3]*a4[4*3+q];
        }

        for(int p=0; p<4; p++)
            for(int q=0; q<4; q++)
                qq[p]+=vi[4*q+p];

        dlqq=num::length(qq);

        for(int i=3; i >= 0; i--)
            if( wek[i]==1 ) wi = i;

        for(int j=0; j < 4; j++)
            q0.push_back(qq[j]);

        for(int i=0; i < 4; i++)
            if( i!=wi )
            {
                for(int j=0; j<4; j++)
                    q0.push_back(a4[4*i+j]);
            }


        int j=0;        // rotation number
        for(int i=3; i>0; i--)
            q0=num::rot(q0, 0, j++, i);

        j=0;
        for(int i=3; i>1; i--)
            q0=num::rot(q0, 1, j++, i);

        q0=num::rot(q0, 2, 0, 3);
        q0=num::mult(q0, num::diagonal(dwp/dlqq, 4), 4, 4, 4);
    }
    else{
        std::string error="Length of Fourier transform need to be equal\n"
        "or longer than  data length.\n";
        throw std::domain_error(error);
    }

    return q0;
}

double DensityS1::density(const std::vector<double>& basis) const
{
    return std::abs( pow(M_PI, 2.0)/( 2.0*num::det(basis, 4) ) ); //M_PI -std=c99 acos (-1.0)
}

double DensityS1::density(double c0=0.75,  unsigned int nfft=524288, unsigned int data_length=344656) const
{
    std::vector<double> q0(grid_prim(c0, nfft, data_length));                  // 2^19 = 524288
    return density(q0);
}
