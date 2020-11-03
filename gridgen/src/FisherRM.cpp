///////////////////////////////////////
// Name:        FisherRM.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
// Modification:21/06/2020
///////////////////////////////////////

#include "FisherRM.h"
#include "num.h"
#include <cmath>
#include <stdexcept>

//#include <iostream>

FisherRM::FisherRM():m_ephemeris_length(0), m_number_of_detectors(-1), m_spindown(-1)
{
    //m_prefisher.resize(0);
    //ctor
}

FisherRM::FisherRM(const std::vector<double>& ephemeris1, const std::vector<double>& ephemeris2, int spindown_)//, int m_spindown
: m_ephemeris_length(ephemeris1.size()), m_number_of_detectors(1), m_spindown(spindown_)
{
    m_prefisher_net.resize(m_number_of_detectors);
    m_sigma.push_back(-1);

    std::vector<double> x(m_ephemeris_length,0);
    std::vector<double> x2(m_ephemeris_length,0);

    for(unsigned int i=0; i<m_ephemeris_length; i++)
        x[i]=(static_cast<double>(i)/static_cast<double>(m_ephemeris_length-1))-0.5;

    for(unsigned int i=0; i<m_ephemeris_length; i++)
        x2[i]=x[i]*x[i];

    switch (m_spindown)
    {
        case 0:
            m_prefisher_net[0].resize(0);
            break;

        case 1:
            m_prefisher_net[0].resize(9);

            m_prefisher_net[0][0]=num::simpson(ephemeris1);
            m_prefisher_net[0][1]=num::simpson(ephemeris2);
            m_prefisher_net[0][2]=num::simpson(x, ephemeris1);
            m_prefisher_net[0][3]=num::simpson(x, ephemeris2);

            m_prefisher_net[0][4]=num::simpson(x2, ephemeris1);
            m_prefisher_net[0][5]=num::simpson(x2, ephemeris2);
            m_prefisher_net[0][6]=num::simpson(ephemeris1, ephemeris1);
            m_prefisher_net[0][7]=num::simpson(ephemeris1, ephemeris2);

            m_prefisher_net[0][8]=num::simpson(ephemeris2, ephemeris2);
            break;

        case 2:
            m_prefisher_net[0].resize(11);

            m_prefisher_net[0][0]=num::simpson(ephemeris1);
            m_prefisher_net[0][1]=num::simpson(ephemeris2);
            m_prefisher_net[0][2]=num::simpson(x, ephemeris1);
            m_prefisher_net[0][3]=num::simpson(x, ephemeris2);

            m_prefisher_net[0][4]=num::simpson(x2, ephemeris1);
            m_prefisher_net[0][5]=num::simpson(x2, ephemeris2);
            m_prefisher_net[0][6]=num::simpson(ephemeris1, ephemeris1);
            m_prefisher_net[0][7]=num::simpson(ephemeris1, ephemeris2);

            m_prefisher_net[0][8]=num::simpson(ephemeris2, ephemeris2);

            std::vector<double> x3(m_ephemeris_length,0);
            for(unsigned int i=0; i<m_ephemeris_length; i++)
                x3[i]=x2[i]*x[i];

            m_prefisher_net[0][9]=num::simpson(x3, ephemeris1);///boole(x3, ephemeris1);
            m_prefisher_net[0][10]=num::simpson(x3, ephemeris2);///boole(x3, ephemeris2);
            break;
    }
}


FisherRM::FisherRM(const std::vector<std::vector<double> >& ephemeris, const std::vector<double>& sigma, int spindown_)//, int m_spindown
: m_ephemeris_length(ephemeris[0].size()), m_sigma(sigma), m_number_of_detectors(sigma.size()), m_spindown(spindown_)
{
    m_prefisher_net.resize(m_number_of_detectors);

    std::vector<double> x(m_ephemeris_length,0);
    std::vector<double> x2(m_ephemeris_length,0);
    /*
    for(unsigned int i=0; i<m_sigma.size(); i++){
        std::cout << "sigma[" << i << "]=" << m_sigma[i] << ", ";
    }
    std::cout << std::endl;
    */
    for(unsigned int i=0; i<m_ephemeris_length; i++)
        x[i]=(static_cast<double>(i)/static_cast<double>(m_ephemeris_length-1))-0.5;

    for(unsigned int i=0; i<m_ephemeris_length; i++)
        x2[i]=x[i]*x[i];

    for(unsigned int i=0; i<m_number_of_detectors; i++){
        switch (m_spindown)
        {
            case 0:
                m_prefisher_net[i].resize(0);
                break;

            case 1:
                m_prefisher_net[i].resize(9);

                m_prefisher_net[i][0]=num::simpson(ephemeris[2*i]);
                m_prefisher_net[i][1]=num::simpson(ephemeris[2*i + 1]);
                m_prefisher_net[i][2]=num::simpson(x, ephemeris[2*i]);
                m_prefisher_net[i][3]=num::simpson(x, ephemeris[2*i + 1]);

                m_prefisher_net[i][4]=num::simpson(x2, ephemeris[2*i]);
                m_prefisher_net[i][5]=num::simpson(x2, ephemeris[2*i + 1]);
                m_prefisher_net[i][6]=num::simpson(ephemeris[2*i], ephemeris[2*i]);
                m_prefisher_net[i][7]=num::simpson(ephemeris[2*i], ephemeris[2*i + 1]);

                m_prefisher_net[i][8]=num::simpson(ephemeris[2*i + 1], ephemeris[2*i + 1]);
                break;

            case 2:
                m_prefisher_net[i].resize(11);

                m_prefisher_net[i][0]=num::simpson(ephemeris[2*i]);
                m_prefisher_net[i][1]=num::simpson(ephemeris[2*i + 1]);
                m_prefisher_net[i][2]=num::simpson(x, ephemeris[2*i]);
                m_prefisher_net[i][3]=num::simpson(x, ephemeris[2*i + 1]);

                m_prefisher_net[i][4]=num::simpson(x2, ephemeris[2*i]);
                m_prefisher_net[i][5]=num::simpson(x2, ephemeris[2*i + 1]);
                m_prefisher_net[i][6]=num::simpson(ephemeris[2*i], ephemeris[2*i]);
                m_prefisher_net[i][7]=num::simpson(ephemeris[2*i], ephemeris[2*i + 1]);

                m_prefisher_net[i][8]=num::simpson(ephemeris[2*i + 1], ephemeris[2*i + 1]);

                std::vector<double> x3(m_ephemeris_length,0);
                for(unsigned int i=0; i<m_ephemeris_length; i++)
                    x3[i]=x2[i]*x[i];

                m_prefisher_net[i][9]=num::simpson(x3, ephemeris[2*i]);///boole(x3, ephemeris[2*i]);
                m_prefisher_net[i][10]=num::simpson(x3, ephemeris[2*i + 1]);///boole(x3, ephemeris[2*i + 1]);
                break;
        }
    }
}

std::vector<double> FisherRM::postrmf(double xi) const
{
    std::vector<double> postfisher_sum;
    switch(m_spindown)
    {
        case 0:
            postfisher_sum.resize(4);
            break;
        case 1:
            postfisher_sum.resize(16);
            break;
        case 2:
            postfisher_sum.resize(25);
            break;
    }


    if(m_number_of_detectors == 1){
        postfisher_sum = postrmf(xi, 0);
    }
    else{
        std::vector<double> sigma_inverse;
        double sigma_inverse_sum = 0.0;

        for(double m: m_sigma){
            sigma_inverse.push_back(1.0/m);
        }
        for(double s: sigma_inverse){
            sigma_inverse_sum+=s;
        }

        std::vector<double> temp;
        for(unsigned int i=0; i<m_number_of_detectors; i++){
            temp = num::multiply_A(postrmf(xi, i), sigma_inverse[i]/sigma_inverse_sum);

            for(unsigned int j=0; j<postfisher_sum.size(); j++){
                postfisher_sum[j]+=temp[j];
            }
        }

    }

    return postfisher_sum;
}

std::vector<double> FisherRM::postrmf(double xi, int detectorNo) const
{
    ///int prefisher_size = m_prefisher_net[detectorNo].size();
    ///std::cout<< "prefisher_size= " << prefisher_size << "\n";
    std::vector<double> postfisher;
    switch(m_spindown)
    {
        case 0:
            postfisher.resize(4);

            postfisher[0]=1.0/12.0;
            postfisher[1]=postfisher[2]=xi/6.0;
            postfisher[3]=1.0/180.0+(xi*xi)/3.0;
            break;

        case 1:
            postfisher.resize(16);

            postfisher[0]=1.0/12.0;
            postfisher[5]=1.0/180.0+(xi*xi)/3.0;
            postfisher[10]=-(m_prefisher_net[detectorNo][0]*m_prefisher_net[detectorNo][0])+m_prefisher_net[detectorNo][6];
            postfisher[15]=-(m_prefisher_net[detectorNo][1]*m_prefisher_net[detectorNo][1])+m_prefisher_net[detectorNo][8];

            postfisher[1]=postfisher[4]=xi/6.0;
            postfisher[2]=postfisher[8]=m_prefisher_net[detectorNo][2];
            postfisher[3]=postfisher[12]=m_prefisher_net[detectorNo][3];

            postfisher[6]=postfisher[9]=m_prefisher_net[detectorNo][4]+(2.0*xi*m_prefisher_net[detectorNo][2])-m_prefisher_net[detectorNo][0]/12.0;
            postfisher[7]=postfisher[13]=m_prefisher_net[detectorNo][5]+(2.0*xi*m_prefisher_net[detectorNo][3])-m_prefisher_net[detectorNo][1]/12.0;
            postfisher[11]=postfisher[14]=-(m_prefisher_net[detectorNo][0]*m_prefisher_net[detectorNo][1])+m_prefisher_net[detectorNo][7];

            break;

        case 2:
            postfisher.resize(25);

            postfisher[0]=1.0/12.0;
            postfisher[6]=1.0/180.0+(xi*xi)/3.0;
            postfisher[18]=-(m_prefisher_net[detectorNo][0]*m_prefisher_net[detectorNo][0])+m_prefisher_net[detectorNo][6];
            postfisher[24]=-(m_prefisher_net[detectorNo][1]*m_prefisher_net[detectorNo][1])+m_prefisher_net[detectorNo][8];

            postfisher[1]=postfisher[5]=xi/6.0;
            postfisher[3]=postfisher[15]=m_prefisher_net[detectorNo][2];
            postfisher[4]=postfisher[20]=m_prefisher_net[detectorNo][3];

            postfisher[8]=postfisher[16]=m_prefisher_net[detectorNo][4]+(2.0*xi*m_prefisher_net[detectorNo][2])-m_prefisher_net[detectorNo][0]/12.0;
            postfisher[9]=postfisher[21]=m_prefisher_net[detectorNo][5]+(2.0*xi*m_prefisher_net[detectorNo][3])-m_prefisher_net[detectorNo][1]/12.0;
            postfisher[19]=postfisher[23]=-(m_prefisher_net[detectorNo][0]*m_prefisher_net[detectorNo][1])+m_prefisher_net[detectorNo][7];

            postfisher[2]=postfisher[10]=1.0/80.0+(xi*xi)/4.0;
            postfisher[7]=postfisher[11]=xi/24.0+pow(xi,3)/2.0;
            postfisher[12]=1.0/448.0+pow(xi,2)/8.0+3.0*pow(xi,4)/4.0;
            postfisher[13]=postfisher[17]=m_prefisher_net[detectorNo][9]+xi*(3.0*(m_prefisher_net[detectorNo][4]+xi*m_prefisher_net[detectorNo][2])-m_prefisher_net[detectorNo][0]/4.0);
            postfisher[14]=postfisher[22]=m_prefisher_net[detectorNo][10]+xi*(3.0*(m_prefisher_net[detectorNo][5]+xi*m_prefisher_net[detectorNo][3])-m_prefisher_net[detectorNo][1]/4.0);
            break;
    }

    return postfisher;
}


std::vector<double> FisherRM::get_prefisher() const
{
    return m_prefisher_net[0];
}

std::vector<double> FisherRM::get_prefisher(int detectorNo) const
{
    int nof = m_number_of_detectors - 1;
    if(detectorNo > nof || detectorNo <= -1){
        size_t size_nof = std::to_string(nof).size();
        std::string snof = std::to_string(nof).substr(0,size_nof);
        throw std::domain_error("Argument of the 'get_prefisher(int )' function need to be in range: <0, " + snof + ">.\n");
    }

    return m_prefisher_net[detectorNo];
}

unsigned int FisherRM::get_ephemeris_length() const
{
    return m_ephemeris_length;
}

unsigned int FisherRM::dim() const
{
    if(m_ephemeris_length==0)
        return 2;
    else
        return m_spindown+3;
}
