///////////////////////////////////////
// Name:        FisherRM.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
// Modification:04/08/2017
///////////////////////////////////////

#include "FisherRM.h"
#include "num.h"
#include <cmath>


FisherRM::FisherRM():m_ephemeris_length(0), spindown(-1)
{
    //m_prefisher.resize(0);
    //ctor
}

FisherRM::FisherRM(const std::vector<double>& ephemeris1, const std::vector<double>& ephemeris2, int spindown_)//, int spindown
: m_ephemeris_length(ephemeris1.size()), spindown(spindown_)
{
    std::vector<double> x(m_ephemeris_length,0);
    std::vector<double> x2(m_ephemeris_length,0);

    for(unsigned int i=0; i<m_ephemeris_length; i++)
        x[i]=(static_cast<double>(i)/static_cast<double>(m_ephemeris_length-1))-0.5;

    for(unsigned int i=0; i<m_ephemeris_length; i++)
        x2[i]=x[i]*x[i];

    switch (spindown)
    {
        case 0:
            m_prefisher.resize(0);
            break;

        case 1:
            m_prefisher.resize(9);

            m_prefisher[0]=num::simpson(ephemeris1);
            m_prefisher[1]=num::simpson(ephemeris2);
            m_prefisher[2]=num::simpson(x, ephemeris1);
            m_prefisher[3]=num::simpson(x, ephemeris2);

            m_prefisher[4]=num::simpson(x2, ephemeris1);
            m_prefisher[5]=num::simpson(x2, ephemeris2);
            m_prefisher[6]=num::simpson(ephemeris1, ephemeris1);
            m_prefisher[7]=num::simpson(ephemeris1, ephemeris2);

            m_prefisher[8]=num::simpson(ephemeris2, ephemeris2);
            break;

        case 2:
            m_prefisher.resize(11);

            m_prefisher[0]=num::simpson(ephemeris1);
            m_prefisher[1]=num::simpson(ephemeris2);
            m_prefisher[2]=num::simpson(x, ephemeris1);
            m_prefisher[3]=num::simpson(x, ephemeris2);

            m_prefisher[4]=num::simpson(x2, ephemeris1);
            m_prefisher[5]=num::simpson(x2, ephemeris2);
            m_prefisher[6]=num::simpson(ephemeris1, ephemeris1);
            m_prefisher[7]=num::simpson(ephemeris1, ephemeris2);

            m_prefisher[8]=num::simpson(ephemeris2, ephemeris2);

            std::vector<double> x3(m_ephemeris_length,0);
            for(unsigned int i=0; i<m_ephemeris_length; i++)
                x3[i]=x2[i]*x[i];

            m_prefisher[9]=num::simpson(x3, ephemeris1);///boole(x3, ephemeris1);
            m_prefisher[10]=num::simpson(x3, ephemeris2);///boole(x3, ephemeris1);
            break;
    }
}

std::vector<double> FisherRM::get_prefisher() const
{
    return m_prefisher;
}

std::vector<double> FisherRM::postrmf(double xi) const
{
    int prefisher_size = m_prefisher.size();
    //std::cout<< "prefisher_size= " << prefisher_size << "\n";
    std::vector<double> postfisher;
    switch(prefisher_size)
    {
        case 0:
            postfisher.resize(4);

            postfisher[0]=1.0/12.0;
            postfisher[1]=postfisher[2]=xi/6.0;
            postfisher[3]=1.0/180.0+(xi*xi)/3.0;
            break;

        case 9:
            postfisher.resize(16);

            postfisher[0]=1.0/12.0;
            postfisher[5]=1.0/180.0+(xi*xi)/3.0;
            postfisher[10]=-(m_prefisher[0]*m_prefisher[0])+m_prefisher[6];
            postfisher[15]=-(m_prefisher[1]*m_prefisher[1])+m_prefisher[8];

            postfisher[1]=postfisher[4]=xi/6.0;
            postfisher[2]=postfisher[8]=m_prefisher[2];
            postfisher[3]=postfisher[12]=m_prefisher[3];

            postfisher[6]=postfisher[9]=m_prefisher[4]+(2.0*xi*m_prefisher[2])-m_prefisher[0]/12.0;
            postfisher[7]=postfisher[13]=m_prefisher[5]+(2.0*xi*m_prefisher[3])-m_prefisher[1]/12.0;
            postfisher[11]=postfisher[14]=-(m_prefisher[0]*m_prefisher[1])+m_prefisher[7];

            break;

        case 11:
            postfisher.resize(25);

            postfisher[0]=1.0/12.0;
            postfisher[6]=1.0/180.0+(xi*xi)/3.0;
            postfisher[18]=-(m_prefisher[0]*m_prefisher[0])+m_prefisher[6];
            postfisher[24]=-(m_prefisher[1]*m_prefisher[1])+m_prefisher[8];

            postfisher[1]=postfisher[5]=xi/6.0;
            postfisher[3]=postfisher[15]=m_prefisher[2];
            postfisher[4]=postfisher[20]=m_prefisher[3];

            postfisher[8]=postfisher[16]=m_prefisher[4]+(2.0*xi*m_prefisher[2])-m_prefisher[0]/12.0;
            postfisher[9]=postfisher[21]=m_prefisher[5]+(2.0*xi*m_prefisher[3])-m_prefisher[1]/12.0;
            postfisher[19]=postfisher[23]=-(m_prefisher[0]*m_prefisher[1])+m_prefisher[7];

            postfisher[2]=postfisher[10]=1.0/80.0+(xi*xi)/4.0;
            postfisher[7]=postfisher[11]=xi/24.0+pow(xi,3)/2.0;
            postfisher[12]=1.0/448.0+pow(xi,2)/8.0+3.0*pow(xi,4)/4.0;
            postfisher[13]=postfisher[17]=m_prefisher[9]+xi*(3.0*(m_prefisher[4]+xi*m_prefisher[2])-m_prefisher[0]/4.0);
            postfisher[14]=postfisher[22]=m_prefisher[10]+xi*(3.0*(m_prefisher[5]+xi*m_prefisher[3])-m_prefisher[1]/4.0);
            break;
    }

    return postfisher;
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
        return spindown+3;
}
