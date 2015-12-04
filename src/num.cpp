///////////////////////////////////////
// Name:        num.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#include "num.h"
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <queue>
#include <cmath>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

double num::simpson(const std::vector<double>& v1)
{
    unsigned int n=v1.size();
    if(n>=5)
    {
        double dx=1.0/static_cast<double>(n-1);
        double sum=0.0, sum1=0.0, sum2=0.0, a, b, c;
        a=dx*(4.0/3.0);
        b=dx*(2.0/3.0);
        c=dx/3.0;
        if(n%2==1)
        {
            for(unsigned int i=1; i<n-1; i+=2)
                sum1+=v1[i];

            for(unsigned int i=2; i<n-2; i+=2)
                sum2+=v1[i];

            sum1*=a;
            sum2*=b;

            sum+=sum1;
            sum+=sum2;

            sum+=c*v1[0];
            sum+=c*v1[n-1];
        }
        else{
            for(unsigned int i=1; i<n-2; i+=2)
                sum1+=v1[i];

            for(unsigned int i=2; i<n-3; i+=2)
                sum2+=v1[i];

            sum1*=a;
            sum2*=b;

            sum+=sum1;
            sum+=sum2;

            sum+=c*v1[0];
            sum+=c*v1[n-2];

            sum+=(dx*0.5)*(v1[n-2]+v1[n-1]);
        }

        return sum;
    }
    else
    {
        throw std::domain_error("To short data input (Simpson's algorithm need at last 5 points).\n");
    }

}

double num::simpson(const std::vector<double>& v1, const std::vector<double>& v2)
{
    unsigned int n=v1.size();
    if(n>=5)
    {
        double dx=1.0/static_cast<double>(n-1);
        double sum=0.0, sum1=0.0, sum2=0.0, a, b, c;
        a=dx*(4.0/3.0);
        b=dx*(2.0/3.0);
        c=dx/3.0;
        if(n%2==1)
        {
            for(unsigned int i=1; i<n-1; i+=2)
                sum1+=v1[i]*v2[i];

            for(unsigned int i=2; i<n-2; i+=2)
                sum2+=v1[i]*v2[i];

            sum1*=a;
            sum2*=b;

            sum+=sum1;
            sum+=sum2;

            sum+=c*v1[0]*v2[0];
            sum+=c*v1[n-1]*v2[n-1];
        }
        else{
            for(unsigned int i=1; i<n-2; i+=2)
                sum1+=v1[i]*v2[i];

            for(unsigned int i=2; i<n-3; i+=2)
                sum2+=v1[i]*v2[i];

            sum1*=a;
            sum2*=b;

            sum+=sum1;
            sum+=sum2;

            sum+=c*v1[0]*v2[0];
            sum+=c*v1[n-2]*v2[n-2];

            sum+=(dx*0.5)*(v1[n-2]*v2[n-2]+v1[n-1]*v2[n-1]);
        }

        return sum;
    }
    else
    {
        throw std::domain_error("To short data input (Simpson's algorithm need at last 5 points).\n");
    }
}



bool num::IsBTZ (double i) {
    return (i>=0);
}


int num::sign(unsigned int i, unsigned int j)
{
    int s=0;
    if((i+j)%2==0)
        s=1;
    else
        s=-1;

    return s;
}

int num::sign(double val)
{
    int temp=0;
    if(val>0)
        temp=1;
    else
        temp=-1;

    return temp;
}

double num::length(const std::vector<double>& v)
{
    double sqsum=0.0;
    for(unsigned int i=0; i < v.size(); i++)
        sqsum+=v[i]*v[i];

    return sqrt(sqsum);
}

std::vector<double>  num::norm(const std::vector<double>& v)
{
    double length=num::length(v);
    std::vector<double> unit(v);
    for(unsigned int i=0; i < v.size(); i++)
        unit[i]/=length;

    return unit;
}

std::vector<double> num::minor2(const std::vector<double>& vc, unsigned int p, unsigned int q, unsigned int dim)
{
    std::vector<double> mr;
    mr.reserve((dim-1)*(dim-1));

    for(unsigned int i=0; i<dim; i++)
    {
        for(unsigned int j=0; j<dim; j++)
        {
            if(i!=p&&j!=q)
            {
                mr.push_back(vc[dim*i+j]);
                //std::cout<<vc[dim*i+j] << " | ";
            }

        }
        //std::cout<<std::endl;
    }

    return mr;
}

double num::det(const std::vector<double> vd, unsigned int dim)
{
    double sum=0.0;
    if(dim==2)
    {
        sum=vd[0]*vd[3]-vd[1]*vd[2];
    }
    else
    {
        for(unsigned int j=0; j<dim; j++)
            sum+=sign(0,j)*vd[j]*det(minor2(vd, 0, j, dim), dim-1);
    }

    return sum;
}

std::vector<double> num::cofactor(const std::vector<double>& cm, unsigned int dim)
{
    std::vector<double> tn;
    tn.reserve((dim-1)*(dim-1));

    for(unsigned int p=0; p<dim; p++)
    {
        for(unsigned int q=0; q<dim; q++)
        {
            tn.push_back( sign(p, q) * det(minor2(cm, p, q, dim), dim-1)) ;
            //std::cout<<tn[dim*p+q] << " | ";
        }
        //std::cout<<std::endl;
    }

    return tn;
}

std::vector<double> num::transpose(const std::vector<double>& tm, unsigned int k, unsigned int m)//&
{
    std::vector<double> tn;
    tn.reserve(k*m);

    for(unsigned int j=0; j<m; j++)
    {
        for(unsigned int i=0; i<k; i++)
        {
            tn.push_back( tm[m*i+j] ) ;
        }
    }

/*
    for(int i=0; i<k; i++)
    {
        for(int j=0; j<m; j++)
        {
            std::cout<<tn[m*i+j] << " | ";
        }
        std::cout<< "t" << std::endl;
    }
*/
    return tn;
}

std::vector<double> num::inverse(const std::vector<double>& im, unsigned int dim)//&
{
    std::vector<double> tn(transpose(cofactor(im, dim), dim, dim));
    double ddet=0.0;
    for(unsigned int i=0; i<dim; i++)
        ddet+=im[i]*tn[dim*i];

    //std::cout<<"ddet="<< ddet << " \n ";

    for(unsigned int i=0; i<dim; i++)
    {
        for(unsigned int j=0; j<dim; j++)
        {
            tn[dim*i+j]/=ddet;
            //std::cout<<tn[dim*i+j] << " | ";
        }
        //std::cout<<std::endl;
    }

    return tn;
}

std::vector<double> num::cholesky(const std::vector<double>& v, unsigned int dim)//&
{
    std::vector<double> chl(dim*dim, 0.0);

    for(unsigned int j=0; j<dim; j++)
    {
        for(unsigned int i=j; i<dim; i++)
        {
            double s=0.0;
            unsigned int ni=dim*i, nj=dim*j;

			for(unsigned int k=0; k<j; k++){
				s+=chl[ni+k]*chl[nj+k];
			}

			if(i==j)
                chl[ni+j]=sqrt(v[ni+j]-s);
            else
                chl[ni+j]=(v[ni+j]-s)/chl[nj+j];

		}
	}

	return chl;
}

std::vector<double> num::mult(const std::vector<double>& a, const std::vector<double>& b,
                              unsigned int col1, unsigned int row1, unsigned int row2)//&&
{
    std::vector<double> t;//t(col1*row2,0.0);
    t.reserve(col1*row2);

    for(unsigned int i=0; i<col1; i++)
    {
        for(unsigned int j=0; j<row2; j++)
        {
            double sum=0.0;
            for(unsigned int p=0; p<row1; p++)
            {
                sum+=a[row1*i+p]*b[row2*p+j];
            }
            t.push_back(sum);
            //t[row2*i+j]=sum;
        }
    }

/*
    for(int i=0; i<col1; i++)
    {
        for(int j=0; j<row2; j++)
        {
            std::cout<<t[row2*i+j] << " | ";
        }
        std::cout<< "mult" << std::endl;
    }
*/

    return t;
}

std::vector<double> num::diagonal(double f, unsigned int dim)
{
    std::vector<double> t(dim*dim, 0.0);

    t[0]=f;
    for(unsigned int i=1; i<dim; i++)
        t[dim*i+i]=1.0;

    return t;
}

double num::epsilon()
{
    return 1.0/pow(10.0,12.0);
}

double num::chop(double to_chop, double precision)
{
    if(to_chop<precision && to_chop>-precision)
        to_chop=0.;

    return to_chop;
}
/*
std::vector<double> num::chop(std::vector<double> to_chop, bool mode, double precision)
{
    if( mode )
        for(unsigned int i=0; i<to_chop.size(); i++)
            if(to_chop[i]<precision && to_chop[i]>-precision)
                to_chop[i]=0.;

    return to_chop;
}
*/
void num::chop(std::vector<double>& to_chop, double precision)
{
    for(unsigned int i=0; i<to_chop.size(); i++)
        if(to_chop[i]<precision && to_chop[i]>-precision)
            to_chop[i]=0.;

}


/// Matrix rotation functions:

std::vector<double> num::scb(const std::vector<double>& cathetus)
{
    std::vector<double> tsc(2,0.0);   // to store values of sin and cos functions
    double c=cathetus[0], d=cathetus[1];
    if(std::abs(d)<num::epsilon())
    {
        tsc[0]=0.;
        tsc[1]=1.;
    }
    else{
        double cd2=c*c+d*d;
        double cd=sqrt(cd2);
        if(cd<num::epsilon())
        {
            tsc[0]=0.;
            tsc[1]=1.;
        }
        else{
            tsc[0]=std::abs(d/cd);
            tsc[1]=std::abs(c/cd);
        }
    }

    return tsc;
}

std::vector<double> num::rotright(const std::vector<double>& w, int axis)
{
    std::vector<double> b(scb(w));
    std::vector<double> tr(16, 0.0);
    if(axis<4 && axis>-1)
    {
        if(axis==0)
        {
            tr[0]=tr[5]=1.0;
            tr[10]=b[1]; tr[11]=b[0];
            tr[14]=-b[0]; tr[15]=b[1];
        }

        if(axis==1)
        {
            tr[0]=tr[15]=1.0;
            tr[5]=b[1]; tr[6]=b[0];
            tr[9]=-b[0]; tr[10]=b[1];
        }

        if(axis==2)
        {
            tr[10]=tr[15]=1.0;
            tr[0]=b[1]; tr[1]=b[0];
            tr[4]=-b[0]; tr[5]=b[1];
        }
    }

    return tr;
}


std::vector<double> num::rotleft(const std::vector<double>& w, int axis)
{
    std::vector<double> b(scb(w));
    std::vector<double> tr(16, 0.0);
    if(axis<4 && axis>-1)
    {
        if(axis==0)
        {
            tr[0]=tr[5]=1.0;
            tr[10]=b[1]; tr[11]=-b[0];
            tr[14]=b[0]; tr[15]=b[1];
        }

        if(axis==1)
        {
            tr[0]=tr[15]=1.0;
            tr[5]=b[1]; tr[6]=-b[0];
            tr[9]=b[0]; tr[10]=b[1];
        }

        if(axis==2)
        {
            tr[10]=tr[15]=1.0;
            tr[0]=b[1]; tr[1]=-b[0];
            tr[4]=b[0]; tr[5]=b[1];
        }
    }

    return tr;
}

std::vector<double> num::rot(const std::vector<double>& v, int u, int axis, int q)
{
    std::vector<double> t(2,0.0);
    std::vector<double> r;
    t[0]=v[4*u+q-1];
    t[1]=v[4*u+q];

    r=num::mult( v, num::transpose( rotright(t, axis), 4, 4 ), 4, 4, 4 );
    if(std::abs(r[4*u+q])>num::epsilon())
        r=num::mult( v, num::transpose( rotleft(t, axis), 4, 4 ), 4, 4, 4 );

    return r;
}


/// Fourier frequency ('_prim' - at space with hyper-spheres)

double num::delta_omega_zero(unsigned int nfft, unsigned int data_length)
{
    return 2.*M_PI*static_cast<double>(data_length)/static_cast<double>(nfft);
}

double num::delta_omega_zero_prim(double c0, unsigned int nfft, unsigned int data_length)
{
    double omega0 = num::delta_omega_zero(nfft, data_length);
    double omega0_prim = omega0/(2.*sqrt(3.0*(1.0-c0)));

    return omega0_prim;
}


/// Experimental functions:

double num::simpson_adv(const std::vector<double>& v1)
{
    unsigned int n=v1.size();
    if(n>=5)
    {
        double dx=1.0/static_cast<double>(n-1);
        double sum=0.0, h=dx/3.0;
        if(n%2==1)
        {
            std::vector<double> vt;
            vt.resize(n);

            for(unsigned int i=1; i<n-1; i+=2)
                vt[i]=4.0*v1[i];

            for(unsigned int i=2; i<n-2; i+=2)
                vt[i]=2.0*v1[i];

            vt[0]=v1[0];
            vt[n-1]=v1[n-1];

            sum+=num::sum(vt);

            sum*=h;
        }
        else{

            std::vector<double> vt;
            vt.resize(n-1);

            for(unsigned int i=1; i<n-2; i+=2)
                vt[i]=4.0*v1[i];

            for(unsigned int i=2; i<n-3; i+=2)
                vt[i]=2.0*v1[i];

            vt[0]=v1[0];
            vt[n-2]=v1[n-2];

            sum+=num::sum(vt);

            sum+=(0.5)*(v1[n-2]+v1[n-1]);
            sum*=h;
        }

        return sum;
    }
    else
    {
        throw std::domain_error("To short data input (Simpson's algorithm need at last 5 points).\n");
    }
}

double num::simpson_adv(const std::vector<double>& v1, const std::vector<double>& v2)
{

    unsigned int n=v1.size();
    if(n>=5)
    {
        double dx=1.0/static_cast<double>(n-1);
        double sum=0.0, h=dx/3.0;
        if(n%2==1)
        {
            std::vector<double> vt;
            vt.resize(n);

            for(unsigned int i=1; i<n-1; i+=2)
                vt[i]=4.0*v1[i]*v2[i];

            for(unsigned int i=2; i<n-2; i+=2)
                vt[i]=2.0*v1[i]*v2[i];

            vt[0]=v1[0]*v2[0];
            vt[n-1]=v1[n-1]*v2[n-1];

            sum+=num::sum(vt);

            sum*=h;
        }
        else{

            std::vector<double> vt;
            vt.resize(n-1);

            for(unsigned int i=1; i<n-2; i+=2)
                vt[i]=4.0*v1[i]*v2[i];

            for(unsigned int i=2; i<n-3; i+=2)
                vt[i]=2.0*v1[i]*v2[i];

            vt[0]=v1[0]*v2[0];
            vt[n-2]=v1[n-2]*v2[n-2];

            sum+=num::sum(vt);

            sum+=(0.5)*(v1[n-2]*v2[n-2]+v1[n-1]*v2[n-1]);
            sum*=h;
        }

        return sum;
    }
    else
    {
        throw std::domain_error("To short data input (Simpson's algorithm need at last 5 points).\n");
    }
}

double num::sum(std::vector<double>& v)
{
    double sum=0.0;
    sort(v.begin(), v.end());
    std::vector<double>::const_iterator it;
    it=find_if(v.begin(), v.end(), IsBTZ);
    std::priority_queue<double, std::vector<double>, std::greater<double> > pqp;
    std::priority_queue<double, std::vector<double>, std::less<double> > pqm;
    for(std::vector<double>::const_iterator i=it; i!=v.end(); i++)
    {
        pqp.push(*i);
        //std::cout<<pqp.top() <<"====\n";
    }
    for(std::vector<double>::const_iterator i=v.begin(); i!=it; i++)
    {
        pqm.push(*i);
        //std::cout<<pqm.top() <<"\n";
    }

    while(pqp.size() > 1)
    {
        double x = pqp.top(); pqp.pop();
        double y = pqp.top(); pqp.pop();
        //std::cout<<"["<<x<<"]["<<y<<"]plus\n";
        pqp.push(x+y);
    }

    while(pqm.size() > 1)
    {
        double x = pqm.top(); pqm.pop();
        double y = pqm.top(); pqm.pop();
        //std::cout<<"["<<x<<"]["<<y<<"]mines\n";
        pqm.push(x+y);
    }

    if(pqp.size() == 1)
        sum+=pqp.top();

    if(pqm.size() == 1)
        sum+=pqm.top();

    return sum;
}

