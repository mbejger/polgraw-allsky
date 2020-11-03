///////////////////////////////////////
// Name:        DensityS1DS.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
// Modification:07/09/2017 A.P.
///////////////////////////////////////

#include "DensityS1DS.h" //
#include "Node.h"
#include "num.h"
#include <string>
#include <stdexcept>    // domain_error
#include <cmath>        // sqrt
#include <algorithm>    // sort

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif


#include <iostream>

DensityS1DS::DensityS1DS(unsigned int spindown_): m_spindown(spindown_)
{
    m_dim = sd2dim(spindown_); //ctor
}

std::vector<double> DensityS1DS::grid_prim(double c0, unsigned int nfft, unsigned int data_length) const
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

        std::vector<double> an_star = num::a_star_r1(m_dim);///a4->an_star///(16, 0.0);
        double distance2first = an_star[0];

        ///std::cout << "distance2first=" << distance2first << std::endl;
        ///
        /*std::cout << "an_star_r1:" << std::endl;
        for(unsigned int i=0; i<m_dim; i++){
            for(unsigned int j=0; j<m_dim; j++){
                std::cout << an_star[m_dim*i+j] << " ";
            }
            std::cout << std::endl;
        }
        */
        ///
        /*
        an_star[0]=sqrt(2.0);

        an_star[4]=an_star[8]=an_star[12]=1.0/(2.0 * sqrt(2));
        an_star[5]=sqrt(7.5)/2.0;   // sqrt(15.0/2.0)

        an_star[9]=an_star[13]=-sqrt(5.0/6.0)/2.0;
        an_star[10]=sqrt(5.0/3.0);

        an_star[14]=-an_star[10]/2.0;
        an_star[15]=sqrt(5)/2.0;
        */

        unsigned int wi = 0;
        int hd = static_cast<unsigned int>(ceil(dwp/distance2first))+3;
        std::vector<int> wek;
        std::vector<double> qq(m_dim, 0.0);

        std::vector<Node> vi_temp;
        Node vtp_temp = Node(pow(10,12), std::vector<int> (m_dim, 0));

        switch(m_spindown){
            case 1:
            for(int i=0; i<=hd; i++){
                for(int j=0; j<=hd; j++){
                    if( i==1 || j==1)
                    {
                        std::vector<double> vi(m_dim*m_dim, 0.0);
                        for(unsigned int q=0; q<m_dim; q++)
                        {
                            vi[q]=i*an_star[q];
                            vi[m_dim*1+q]=j*an_star[m_dim*1+q];
                            /// test:
                            //std::cout << "vi[q]= " << vi[q] << ", vi[m_dim*1+q]= " << vi[m_dim*1+q] << std::endl;
                        }

                        std::vector<double> vitemp(m_dim,0.0);
                        for(unsigned int p=0; p<m_dim; p++)
                            for(unsigned int q=0; q<m_dim; q++)
                                vitemp[p]+=vi[m_dim*q+p];

                        double dtemp=num::length(vitemp);
                        if( dtemp >= dwp )//= dwp
                        {
                            std::vector<int> vtemp(m_dim,0);//={i,j,k,l};
                            vtemp[0]=i;
                            vtemp[1]=j;

                            Node vtp=Node(dtemp-dwp, vtemp);
                            ///vi_temp.push_back(vtp);///->(*)
                            if(vtp.m_distance<vtp_temp.m_distance)
                                vtp_temp = vtp;
                        }
                    }
                }
            }
            break;

            case 2:
            for(int i=0; i<=hd; i++){
                for(int j=0; j<=hd; j++){
                    for(int k=0; k<=hd; k++){
                            if( i==1 || j==1 || k==1)
                            {
                                std::vector<double> vi(m_dim*m_dim, 0.0);
                                for(unsigned int q=0; q<m_dim; q++)
                                {
                                    vi[q]=i*an_star[q];
                                    vi[m_dim*1+q]=j*an_star[m_dim*1+q];
                                    vi[m_dim*2+q]=k*an_star[m_dim*2+q];
                                }

                                std::vector<double> vitemp(m_dim,0.0);
                                for(unsigned int p=0; p<m_dim; p++)
                                    for(unsigned int q=0; q<m_dim; q++)
                                        vitemp[p]+=vi[m_dim*q+p];

                                double dtemp=num::length(vitemp);
                                if( dtemp >= dwp )//= dwp
                                {
                                    std::vector<int> vtemp(m_dim,0);//={i,j,k,l};
                                    vtemp[0]=i;
                                    vtemp[1]=j;
                                    vtemp[2]=k;

                                    Node vtp=Node(dtemp-dwp, vtemp);
                                    ///vi_temp.push_back(vtp);///->(*)
                                    if(vtp.m_distance<vtp_temp.m_distance)
                                        vtp_temp = vtp;
                                }
                            }
                    }
                }
            }
            break;

            case 3:
            for(int i=0; i<=hd; i++){
                for(int j=0; j<=hd; j++){
                    for(int k=0; k<=hd; k++){
                        for(int l=0; l<=hd; l++)
                        {
                            if( i==1 || j==1 || k==1 || l==1 )
                            {
                                std::vector<double> vi(m_dim*m_dim, 0.0);
                                for(unsigned int q=0; q<m_dim; q++)
                                {
                                    vi[q]=i*an_star[q];
                                    vi[m_dim*1+q]=j*an_star[m_dim*1+q];
                                    vi[m_dim*2+q]=k*an_star[m_dim*2+q];
                                    vi[m_dim*3+q]=l*an_star[m_dim*3+q];
                                }

                                std::vector<double> vitemp(m_dim,0.0);
                                for(unsigned int p=0; p<m_dim; p++)
                                    for(unsigned int q=0; q<m_dim; q++)
                                        vitemp[p]+=vi[m_dim*q+p];

                                double dtemp=num::length(vitemp);
                                if( dtemp >= dwp )//= dwp
                                {
                                    std::vector<int> vtemp(m_dim,0);//={i,j,k,l};
                                    vtemp[0]=i;
                                    vtemp[1]=j;
                                    vtemp[2]=k;
                                    vtemp[3]=l;

                                    Node vtp=Node(dtemp-dwp, vtemp);
                                    ///vi_temp.push_back(vtp);///->(*)
                                    if(vtp.m_distance<vtp_temp.m_distance)
                                        vtp_temp = vtp;
                                }
                            }
                        }
                    }
                }
            }
            break;

            case 4:
            for(int i=0; i<=hd; i++){
                for(int j=0; j<=hd; j++){
                    for(int k=0; k<=hd; k++){
                        for(int l=0; l<=hd; l++){
                            for(int m=0; m<=hd; m++)
                            {
                                if( i==1 || j==1 || k==1 || l==1 || m==1)
                                {
                                    std::vector<double> vi(m_dim*m_dim, 0.0);
                                    for(unsigned int q=0; q<m_dim; q++)
                                    {
                                        vi[q]=i*an_star[q];
                                        vi[m_dim*1+q]=j*an_star[m_dim*1+q];
                                        vi[m_dim*2+q]=k*an_star[m_dim*2+q];
                                        vi[m_dim*3+q]=l*an_star[m_dim*3+q];
                                        vi[m_dim*4+q]=m*an_star[m_dim*4+q];
                                    }

                                    std::vector<double> vitemp(m_dim,0.0);
                                    for(unsigned int p=0; p<m_dim; p++)
                                        for(unsigned int q=0; q<m_dim; q++)
                                            vitemp[p]+=vi[m_dim*q+p];

                                    double dtemp=num::length(vitemp);
                                    if( dtemp >= dwp )//= dwp
                                    {
                                        std::vector<int> vtemp(m_dim,0);//={i,j,k,l};
                                        vtemp[0]=i;
                                        vtemp[1]=j;
                                        vtemp[2]=k;
                                        vtemp[3]=l;
                                        vtemp[4]=m;

                                        Node vtp=Node(dtemp-dwp, vtemp);
                                        ///vi_temp.push_back(vtp);///->(*)
                                        if(vtp.m_distance<vtp_temp.m_distance)
                                            vtp_temp = vtp;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            break;
        } //End of switch(spindown)

        ///std::sort(vi_temp.begin(), vi_temp.end(), compare_Node_distance);///(*)


        ///wek=vi_temp[0].get_coef();
        wek=vtp_temp.get_coef();
        ///

        /*
        /// test:
        for(unsigned int i=0; i<wek.size(); i++)
            std::cout<< "wek[" << i << "]=" << wek[i] << " ";

        std::cout<< "\n A*: " << std::endl;

        for(unsigned int p=0; p<m_dim; p++){
            for(unsigned int q=0; q<m_dim; q++){
                std::cout << an_star[m_dim*p+q] << ", ";//vi[m_dim*p+q];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        */

        ///
        std::vector<double> vi(m_dim*m_dim, 0.0);
        for(unsigned int p=0; p<m_dim; p++)
            for(unsigned int q=0; q<m_dim; q++)
                vi[m_dim*p+q]+=wek[p]*an_star[m_dim*p+q];//vi[m_dim*p+q];

        for(unsigned int p=0; p<m_dim; p++)//vector component
            for(unsigned int q=0; q<m_dim; q++)//vector index
                qq[p]+=vi[m_dim*q+p];//vi[m_dim*q+p];

        dlqq=num::length(qq);

        for(int i=(m_dim-1); i >= 0; i--){
            if( wek[i]==1 ) wi = i;
        }

        for(unsigned int j=0; j < m_dim; j++)
            q0.push_back(qq[j]);

        for(unsigned int i=0; i < m_dim; i++)
            if( i!=wi )
            {
                for(unsigned int j=0; j<m_dim; j++)
                    q0.push_back(an_star[m_dim*i+j]);
            }

        /*
		/// test:
        std::cout<< "\n q0 before rot(): " << std::endl;

        for(unsigned int p=0; p<m_dim; p++){
            for(unsigned int q=0; q<m_dim; q++){
                std::cout << an_star[m_dim*p+q] << ", ";//vi[m_dim*p+q];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
		//
		*/
        q0=num::rot(q0, m_dim);

        /*
        /// test:
		std::cout<< "\n q0 after rot(): " << std::endl;
        for(unsigned int p=0; p<m_dim; p++){
            for(unsigned int q=0; q<m_dim; q++){
                std::cout << an_star[m_dim*p+q] << ", ";//vi[m_dim*p+q];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
		//


		// test:
		std::vector<double> temp11 = num::diagonal(dwp/dlqq, m_dim);
        std::cout<< "\n diagonal(dwp/dlqq, m_dim): " << std::endl;
        for(unsigned int p=0; p<m_dim; p++){
            for(unsigned int q=0; q<m_dim; q++){
                std::cout << temp11[m_dim*p+q] << ", ";//vi[m_dim*p+q];
            }
            std::cout << std::endl;
        }
        std::cout << "m_dim= " << m_dim << std::endl;
        //
        */

        q0=num::multiply_AB(q0, num::diagonal(dwp/dlqq, m_dim), m_dim, m_dim, m_dim);

        /*
        // test:
        std::cout<< "\n q0 after multiply: " << std::endl;
        for(unsigned int p=0; p<m_dim; p++){
            for(unsigned int q=0; q<m_dim; q++){
                std::cout << q0[m_dim*p+q] << ", ";//vi[m_dim*p+q];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        //
        */
    }

    else{
        std::string error="Length of Fourier transform need to be equal\n"
        "or longer than  data length.\n";
        throw std::domain_error(error);
    }

    return q0;
}

double DensityS1DS::density(const std::vector<double>& basis) const
{
    return num::thickness(basis, 1.0, m_dim);
    ///return std::abs( pow(M_PI, 2.0)/( 2.0*num::det(basis, m_dim ) ) ); //M_PI -std=c99 acos (-1.0)
}

double DensityS1DS::density(double c0,  unsigned int nfft, unsigned int data_length) const
{
    std::vector<double> q0(grid_prim(c0, nfft, data_length));// 2^19 = 524288
    return density(q0);
}

// Convert spindown to dimension (directed searches)
unsigned int DensityS1DS::sd2dim(unsigned int spindown) const
{
    if(spindown > 0 && spindown <= num::max_spindownDS)
    {
        return spindown + 1;
    }
    else{
        //std::cerr << "spindown=" << spindown << std::endl;
        std::string error="Allowed spindown values: <1, " + std::to_string(num::max_spindownDS) + ">. \n";
        throw std::domain_error(error);
    }

}


/// Auxiliary function (not used in the program):
// Convert dimension to spindown  (directed searches)
unsigned int DensityS1DS::dim2sd(unsigned int dimension) const
{
    if(dimension > 0 && dimension <= num::max_spindownDS - 1)
    {
        return dimension - 1;
    }
    else{
        std::string error="Allowed dimensions values of matrix An* : <1,";
        std::string maxsd4DS=std::to_string(num::max_spindownDS-1);
        error+=">.\n";
        throw std::domain_error(error);
    }

}

