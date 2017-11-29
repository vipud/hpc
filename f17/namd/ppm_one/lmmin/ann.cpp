#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>
#include <time.h>
#include <sstream>
using namespace std;

#include <config.h>
#include "ann.h"


CAnn::CAnn()
{
	p=NULL;
	srand (time(NULL));
	#pragma acc enter data create(this)
};


CAnn::~CAnn()
{
	#pragma acc exit data delete(p_save_flat)
	#pragma acc exit data delete(this)
	delete(p_save_flat);
	p_save_size = 0;
};


void CAnn::xapplyminmax(double *xx)
{
	double tmin,tmax;
	int i,j,step;
	for(i = 0; i < n_dim; i++)
	{
		tmin=v_min[i];
		tmax=v_max[i];
		for(j = 0; j < n_dat; j++)
		{
			step = j*n_dim+i;
			xx[step]=(xx[step]-tmin)/(tmax-tmin)*2-1;
		}
	}

	return;
};


void CAnn::xapplyminmax()
{
	double tmin,tmax;
	int i,j,step;
	for(i=0;i<n_dim;i++)
	{
		tmin=v_min[i];
		tmax=v_max[i];
		for(j=0;j<n_dat;j++)
		{
			step=j*n_dim+i;
			x[step]=(x[step]-tmin)/(tmax-tmin)*2-1;
		}

	}

	return;
};


#pragma acc routine vector
double CAnn::myfunc_neuron(int ndim, int nneuron, double *xx, const double *p, int offset)
{
	int i;
	double r2;
	r2 = 0;
	const int p2=offset+ndim*nneuron;
	const int p3=p2+nneuron;
	const int p4=p3+nneuron;

	#pragma acc loop reduction(+:r2) vector
	for(int j = 0; j < nneuron; j++)
	{
		double r = 0.0;
		const int p1 = offset+j*ndim;
		#pragma acc loop reduction(+:r) seq
		for(i = 0; i < ndim; i++)
		{
			r += xx[i] * p[p1+i];
		}

		r += p[p2+j];

		if(r<-30 )
			r=0.0;
		else if(r>30)
			r=1.0;
		else
			r=(1-exp(-2*r))/(1+exp(-2*r));

		r = r*p[p3+j];
		r2+=r;
	}

	r2 += p[p4];

	return r2;
};


void CAnn::loadp(double *pdata)
{
	int i,j;
	ifstream fin;
	string part;
	int n_set;
	vector<double> t;



	n_dim=(int)(*pdata);pdata++;
	n_neuron=(int)(*pdata);pdata++;
	n_par=(int)(*pdata);pdata++;

	n_par=n_dim*n_neuron+n_neuron+n_neuron+1;

	p=new double[n_par];
	v_min=new double[n_dim];
	v_max=new double[n_dim];

	for(i=0;i<n_dim;i++)
	{
		v_min[i]=*pdata;pdata++;
		v_max[i]=*pdata;pdata++;
	}

	y_min=*pdata;pdata++;
	y_max=*pdata;pdata++;


	n_set=(int)(*pdata);pdata++;
	p_save_flat = new double[n_set*n_par];
	p_save_size = n_set;
	for(j=0;j<n_set;j++)
	{
		for(i=0;i<n_par;i++)
		{
			p_save_flat[j*n_par + i] = *pdata;
			pdata++;
		}
	}
	#pragma acc enter data copyin(p_save_flat[0:n_set*n_par])
};


double CAnn::predict_one_first( double *xx, int vec_size, double *next, int next_size, CAnn *next_CAnn )
{
	double tt,out;

	if(vec_size!=n_dim)
		return 0;

	n_dat=1;
	xapplyminmax(xx);
	out=0;

	int ndim = n_dim;
	int nneuron = n_neuron;
	int npar = n_par;
	double ymax = y_max;
	double ymin = y_min;
	double *psaveflat = p_save_flat;
	int psavesize = p_save_size;

	#pragma acc enter data copyin(xx[0:n_dim])

	#pragma acc parallel loop gang reduction(+:out) private(tt) \
		present(xx[0:ndim],psaveflat[0:npar]) async
	for(int j=0;j<psavesize;j++)
	{
		tt=CAnn::myfunc_neuron(ndim, nneuron, xx, psaveflat, npar*j);
		tt=(tt+1)/2*(ymax-ymin)+ymin;
		out+=tt;
	}

	next_CAnn->n_dat = 1;
	next_CAnn->xapplyminmax(next);
	#pragma acc enter data copyin(next[0:next_size]) async

	#pragma acc wait

	#pragma acc exit data delete(xx)

	out/=p_save_size;
	return out;
}


double CAnn::predict_one_next( double *xx, int vec_size, double *next, int next_size, CAnn *next_cann )
{
	double tt,out;

	if(vec_size!=n_dim)
		return 0;

	out=0;

	int ndim = n_dim;
	int nneuron = n_neuron;
	int npar = n_par;
	double ymax = y_max;
	double ymin = y_min;
	double *psaveflat = p_save_flat;
	int psavesize = p_save_size;

	#pragma acc parallel loop gang reduction(+:out) private(tt) \
		present(xx[0:ndim],psaveflat[0:npar]) async
	for(int j=0;j<psavesize;j++)
	{
		tt=CAnn::myfunc_neuron(ndim, nneuron, xx, psaveflat, npar*j);
		tt=(tt+1)/2*(ymax-ymin)+ymin;
		out+=tt;
	}

	next_cann->n_dat = 1;
	next_cann->xapplyminmax(next);

	#pragma acc enter data copyin(next[0:next_size]) async

	#pragma acc wait

	#pragma acc exit data delete(xx)

	out/=p_save_size;
	return out;
}


double CAnn::predict_one_last( double *xx, int vec_size)
{
	double tt,out;

	if(vec_size!=n_dim)
		return 0;

	out=0;

	int ndim = n_dim;
	int nneuron = n_neuron;
	int npar = n_par;
	double ymax = y_max;
	double ymin = y_min;
	double *psaveflat = p_save_flat;
	int psavesize = p_save_size;

	#pragma acc parallel loop gang reduction(+:out) private(tt) \
		present(xx[0:ndim],psaveflat[0:npar])
	for(int j=0;j<psavesize;j++)
	{
		tt=CAnn::myfunc_neuron(ndim, nneuron, xx, psaveflat, npar*j);
		tt=(tt+1)/2*(ymax-ymin)+ymin;
		out+=tt;
	}

	#pragma acc exit data delete(xx)

	out/=p_save_size;
	return out;
}


double CAnn::predict_one( double *xx, int vec_size )
{
	double tt,out;

	if(vec_size!=n_dim)
		return 0;

	n_dat=1;
	xapplyminmax(xx);
	out=0;

	int ndim = n_dim;
	int nneuron = n_neuron;
	int npar = n_par;
	double ymax = y_max;
	double ymin = y_min;
	double *psaveflat = p_save_flat;
	int psavesize = p_save_size;

	#pragma acc enter data copyin(xx[0:n_dim])

	#pragma acc parallel loop gang reduction(+:out) private(tt) \
		present(xx[0:ndim],psaveflat[0:npar])
	for(int j=0;j<psavesize;j++)
	{
		tt=CAnn::myfunc_neuron(ndim, nneuron, xx, psaveflat, npar*j);
		tt=(tt+1)/2*(ymax-ymin)+ymin;
		out+=tt;
	}

	#pragma acc exit data delete(xx)

	out/=p_save_size;
	return out;
};

