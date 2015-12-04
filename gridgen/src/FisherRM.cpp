///////////////////////////////////////
// Name:        FisherRM.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
///////////////////////////////////////

#include "FisherRM.h"
#include "num.h"


FisherRM::FisherRM():m_ephemeris_length(0)
{
    //ctor
}

FisherRM::FisherRM(const std::vector<double>& ephemeris1, const std::vector<double>& ephemeris2)
: m_ephemeris_length(ephemeris1.size())
{
    std::vector<double> x(m_ephemeris_length,0);
    std::vector<double> x2(m_ephemeris_length,0);

    for(unsigned int i=0; i<m_ephemeris_length; i++)
        x[i]=(static_cast<double>(i)/static_cast<double>(m_ephemeris_length-1))-0.5;

    for(unsigned int i=0; i<m_ephemeris_length; i++)
        x2[i]=x[i]*x[i];

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

    /// experimental simpson method
/*
    m_prefisher[0]=num::simpson_adv(ephemeris1);
    m_prefisher[1]=num::simpson_adv(ephemeris2);
    m_prefisher[2]=num::simpson_adv(x, ephemeris1);
    m_prefisher[3]=num::simpson_adv(x, ephemeris2);

    m_prefisher[4]=num::simpson_adv(x2, ephemeris1);
    m_prefisher[5]=num::simpson_adv(x2, ephemeris2);
    m_prefisher[6]=num::simpson_adv(ephemeris1, ephemeris1);
    m_prefisher[7]=num::simpson_adv(ephemeris1, ephemeris2);

    m_prefisher[8]=num::simpson_adv(ephemeris2, ephemeris2);
*/

}

std::vector<double> FisherRM::get_prefisher() const
{
    return m_prefisher;
}

std::vector<double> FisherRM::postrmf(double xi) const
{
    std::vector<double> postfisher(16,0.0);
    if(m_ephemeris_length==0)     // without ephemeris 'postrmf' will return 2D fisher reduced matrix
    {
        postfisher.resize(4);

        postfisher[0]=1.0/12.0;
        postfisher[1]=postfisher[2]=xi/6.0;
        postfisher[3]=1.0/180.0+(xi*xi)/3.0;
    }
    else
    {
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
        return 4;
}
