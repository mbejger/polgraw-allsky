///////////////////////////////////////
// Name:        num.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
// Modification: 26/07/2017
///////////////////////////////////////

#include "num.h"
#include <stdexcept>
#include <algorithm>
#include <queue>
#include <cmath>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

///All below integration methods are based on the closed Newton–Cotes’ formulas
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

//Boole's rule:
double num::boole(const std::vector<double>& v1)
{
    unsigned int n=v1.size();
    if(n>=7)
    {
        double dx=1.0/static_cast<double>(n-1);
        double sum=0.0, sum1=0.0, sum2=0.0, sum3=0.0, a, b, c;
        a=dx*(64.0/45.0);
        b=dx*(24.0/45.0);
        c=dx*(14.0/45.0);
        if(n%2==1)
        {
            for(unsigned int i=1; i<n-3; i+=4)
                sum1+=v1[i];

            for(unsigned int i=2; i<n-2; i+=4)
                sum2+=v1[i];

            for(unsigned int i=3; i<n-1; i+=4)
                sum3+=v1[i];

            sum1*=a;
            sum2*=b;
            sum3*=a;

            sum+=sum1;
            sum+=sum2;
            sum+=sum3;

            sum+=c*v1[0];
            sum+=c*v1[n-1];
        }
        else{
            for(unsigned int i=1; i<n-4; i+=4)
                sum1+=v1[i];

            for(unsigned int i=2; i<n-3; i+=4)
                sum2+=v1[i];

            for(unsigned int i=3; i<n-2; i+=4)
                sum3+=v1[i];

            sum1*=a;
            sum2*=b;
            sum3*=a;

            sum+=sum1;
            sum+=sum2;
            sum+=sum3;

            sum+=c*v1[0];
            sum+=c*v1[n-2];

            sum+=(dx*0.5)*(v1[n-2]+v1[n-1]);
        }

        return sum;
    }
    else
    {
        throw std::domain_error("To short data input (Boole's integration rule need at last 7 points).\n");
    }
}

//Boole's rule
double num::boole(const std::vector<double>& v1, const std::vector<double>& v2)
{
    unsigned int n=v1.size();
    if(n>=7)
    {
        double dx=1.0/static_cast<double>(n-1);
        double sum=0.0, sum1=0.0, sum2=0.0, sum3=0.0, a, b, c;
        a=dx*(64.0/45.0);
        b=dx*(24.0/45.0);
        c=dx*(14.0/45.0);
        if(n%2==1)
        {
            for(unsigned int i=1; i<n-3; i+=4)
                sum1+=v1[i]*v2[i];

            for(unsigned int i=2; i<n-2; i+=4)
                sum2+=v1[i]*v2[i];

            for(unsigned int i=3; i<n-1; i+=4)
                sum3+=v1[i]*v2[i];

            sum1*=a;
            sum2*=b;
            sum3*=a;

            sum+=sum1;
            sum+=sum2;
            sum+=sum3;

            sum+=c*v1[0]*v2[0];
            sum+=c*v1[n-1]*v2[n-1];
        }
        else{
            for(unsigned int i=1; i<n-4; i+=4)
                sum1+=v1[i]*v2[i];

            for(unsigned int i=2; i<n-3; i+=4)
                sum2+=v1[i]*v2[i];

            for(unsigned int i=3; i<n-2; i+=4)
                sum3+=v1[i]*v2[i];

            sum1*=a;
            sum2*=b;
            sum3*=a;

            sum+=sum1;
            sum+=sum2;
            sum+=sum3;

            sum+=c*v1[0]*v2[0];
            sum+=c*v1[n-2]*v2[n-2];

            sum+=(dx*0.5)*(v1[n-2]*v2[n-2]+v1[n-1]*v2[n-1]);
        }

        return sum;
    }
    else
    {
        throw std::domain_error("To short data input (Boole's integration rule need at last 7 points).\n");
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

std::vector<double> num::minor1(const std::vector<double>& vc, unsigned int p, unsigned int q, unsigned int dim)
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
            sum+=sign(0,j)*vd[j]*det( minor1(vd, 0, j, dim), dim-1);
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
            tn.push_back( sign(p, q) * det(minor1(cm, p, q, dim), dim-1) ) ;
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

///col1 == length of the column of 1-st matrix 'A',
///row1 == length of the row of 1-st matrix 'A',
///row2 == length of the row of 2-nd matrix 'B'.
std::vector<double> num::multiply_AB(const std::vector<double>& a, const std::vector<double>& b,
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

std::vector<double> num::multiply_A(const std::vector<double>& a, double scalar)
{
    std::vector<double> t(a);
    t.reserve(a.size());

    std::transform(a.begin(), a.end(), t.begin(), std::bind1st(std::multiplies<double>(),scalar));

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

std::vector<double> num::rotation(const std::vector<double>& w, int axis, int left_or_right, unsigned int dim)
{
    std::vector<double> b(scb(w));
    if(left_or_right == -1)//-1 == left, 1 == right
        b[0]*=-1;

    std::vector<double> tr(diagonal(1.0,dim));

    unsigned int axis2 = axis + 1;/// axis == axis1
    if(axis2<dim && axis>-1)
    {

        tr[axis*(dim+1)]=tr[axis2*(dim+1)]=b[1];
        tr[axis*(dim+1)+1]=b[0];
        tr[axis2*(dim+1)-1]=-b[0];

        /*
        std::cout << "tr:\n" << std::endl;
        for(unsigned int i=0; i<dim; i++){
            for(unsigned int j=0; j<dim; j++){
                std::cout<< std::setw(13)<< std::setprecision(12)<< tr[dim*i+j] << " ";
            }
            std::cout << std::endl;
        }
        system("pause");
        */
    }
    else{
        throw std::domain_error("first rotation axis need to be from range <0, dimension - 2>.\n");
    }

    return tr;
}

std::vector<double> num::single_rot(const std::vector<double>& v, int u, int axis, int q, unsigned int dim)
{
    std::vector<double> t(2,0.0);
    std::vector<double> r;
    t[0]=v[dim*u+q-1];
    t[1]=v[dim*u+q];

    r=num::multiply_AB( v, num::transpose( rotation(t, axis, 1, dim), dim, dim ), dim, dim, dim );
    if(std::abs(r[dim*u+q])>num::epsilon())
        r=num::multiply_AB( v, num::transpose( rotation(t, axis, -1, dim), dim, dim ), dim, dim, dim );

    /*
    std::cout << "single_rot:" << t[0] << " " << t[1] << std::endl;
    for(unsigned int i=0; i<dim; i++){
        for(unsigned int j=0; j<dim; j++){
            std::cout<< std::setw(13)<< std::setprecision(12)<< r[dim*i+j] << " ";
        }
        std::cout << std::endl;
    }
    system("pause");
    */

    return r;
}

std::vector<double> num::rot(const std::vector<double>& q0, unsigned int dim)
{
    std::vector<double> r0(q0);
    for(unsigned int u=0; u<dim-1; u++){
        int axis=dim-2;        // rotation number
        for(unsigned int q=dim-1; q>u; q--)
        {
            /*
            std::cout << "q0: rot nr: -->" << u << "\n";
            for(unsigned int k=0; k<dim; k++){
                for(unsigned int p=0; p<dim; p++){
                    std::cout<< std::setw(13)<< std::setprecision(12)<< r0[k*dim+p] << " ";
                }
                std::cout << "\n";
            }
            */

            r0=num::single_rot(r0, u, axis--, q, dim);
        }


       /*
        int j=0;        // rotation number
        for(int i=3; i>0; i--)
            q0=num::rot(q0, 0, j++, i);
        j=0;
        for(int i=3; i>1; i--){
            q0=num::rot(q0, 1, j++, i);
        }
        q0=num::rot(q0, 2, 0, 3);
        */
    }

    //r0=num::single_rot(r0, dim-2, 2, dim-1, dim);

    return r0;
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

/// Sphere coverings:

std::vector<double> num::gram(unsigned int dim)
{
    std::vector<double> mr;
    mr.reserve((dim-1)*(dim-1));

    for(unsigned int i=0; i<dim; i++)
    {
        for(unsigned int j=0; j<dim; j++)
        {
            if(i==j)
            {
                mr.push_back(dim);
                //std::cout<<dim << " | ";
            }
            else{
                mr.push_back(-1);
                //std::cout<<-1 << " | ";
            }

        }
        //std::cout<<std::endl;
    }

    return mr;
}

std::vector<double> num::a_star(unsigned int dim)
{
    std::vector<double> aNstar(cholesky(gram(dim), dim));

	return aNstar;
}

std::vector<double> num::a_star_r1(unsigned int dim)
{
    double rescale = 1.0/covering_radius(dim);
    std::vector<double> aNstarR1(multiply_A(a_star(dim), rescale));

    ///below step is not quite necessary but in some cases
    ///produce better (lower) thickness:
    for(unsigned int i=0; i<dim; i++)
        aNstarR1[dim*i]=fabs(aNstarR1[dim*i]);//fabs()

    /// possibly all elements can not be positive
    ///(algorithm produce then bigger thickness)!
	return aNstarR1;
}

double num::sphere_volume(unsigned int dim, double radius){
    double fr = static_cast<double> (dim/2.0);
    double dpow = pow(M_PI, fr)*pow(radius, dim);
    double dgamma = std::tgamma(fr+1.0);

    return dpow/dgamma;
}

double num::covering_radius(unsigned int dim){
    double d1=sqrt(dim*(dim+2.0));
    double d2=2.0*sqrt(3.0);

    return d1/d2;
}

double num::thickness(const std::vector<double>& astar, double radius, unsigned int dim){
    double sv = sphere_volume(dim, radius);
    double ddet = num::det(astar, dim);

    return sv/fabs(ddet);
}
