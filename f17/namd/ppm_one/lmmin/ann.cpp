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
};

CAnn::~CAnn()
{
	#pragma acc exit data delete(p_save_flat, v_min, v_max)
	#pragma acc exit data delete(this)
	delete(p_save_flat);
	p_save_size = 0;
};


void CAnn::randperm(int n,int *perm)
{
	int i, j, t;

	for(i=0; i<n; i++)
		perm[i] = i;
	for(i=0; i<n; i++) {
		j = rand()%(n-i)+i;
		t = perm[j];
		perm[j] = perm[i];
		perm[i] = t;
	}
}

int CAnn::train(int ntrain, string name1,string name2,int n, double percent)
{
	int i;

	if(loadx(name1)!=0)
	{
		cout<<"Load X fail!"<<endl;
		return 1;
	}
	else if(loady(name2)!=0)
	{
		cout<<"Load Y fail!"<<endl;
		return 2;
	}
	mapminmax();


	for(i=0;i<ntrain;i++)
	{
		initp(n);
		train_it(percent);
		push_p();
		if(i==0)  savehead("temp_save.dat",ntrain);
		savep("temp_save.dat",i);
	}
	return 0;
}

int CAnn::train_md(int ntrain, string name1,string name2,int n, double percent)
{
	int i;

	if(loadx_md(name1)!=0)
	{
		cout<<"Load X fail!"<<endl;
		return 1;
	}
	else if(loady(name2)!=0)
	{
		cout<<"Load Y fail!"<<endl;
		return 2;
	}
	mapminmax_md();


	for(i=0;i<ntrain;i++)
	{
		initp(n);
		train_it_md(percent);
		push_p();
		if(i==0)  savehead("temp_save.dat",ntrain);
		savep("temp_save.dat",i);
	}
	return 0;
}



void CAnn::close()
{
	n_dim=0;
	n_dat=0;
	n_neuron=0;
	n_par=0;
	n_conf=0;

	delete[] p;
	delete[] x;
	delete[] y;
	delete[] v_min;
	delete[] v_max;

	return;
}


void CAnn::push_p(void)
{
	int i;
	vector<double> t;

	t.clear();

	double *p_save_tmp = new double[(p_save_size+1) * n_par];
	memcpy(p_save_tmp, p_save_flat, p_save_size*n_par*sizeof(double));
	delete(p_save_flat);
	p_save_flat = p_save_tmp;

	for(i=0;i<n_par;i++)
	{
		p_save_flat[p_save_size*n_par + i] = p[i];
	}
	p_save_size += 1;

}




void CAnn::initp(int n)
{
	int i;

	n_neuron=n;
	n_par=n_dim*n_neuron+n_neuron+n_neuron+1;

	if(p!=NULL)
		delete [] p;
	p=new double[n_par];


	for(i=0;i<n_par;i++)
	{
		p[i]=((double)rand())/RAND_MAX*2-1;
	}

	for(i=0;i<n_dim*n_neuron;i++)
		p[i]/=sqrt((double) n_dim);
	for(i=n_dim*n_neuron+n_neuron;i<n_dim*n_neuron+n_neuron*2;i++)
		p[i]/=sqrt((double) n_neuron);

	return;
};

bool CAnn::loadx(vector<vector<double> > xx)
{
	int i,j;
	n_dat=xx.size();
	n_dim=xx.at(0).size();

	x=new double[n_dat*n_dim];

	for(i=0;i<n_dat;i++)
	{
		for(j=0;j<n_dim;j++)
		{
			x[i*n_dim+j]=xx.at(i).at(j);
		}
	}

	return 1;
}


bool  CAnn::loadx_md(string name)
{
	int i;
	ifstream fin;
	string line,part;
	istringstream iss;
	vector<double> x_line;
	vector<double> xx;
	bool breturn=0;
	bool bfirst=1;

	fin.open(name.c_str());
	n_dat=0;
	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		if(bfirst==1)
		{
			iss>>n_conf;
			bfirst=0;
			continue;
		}
		x_line.clear();
		while(iss>>part)
		{
			x_line.push_back(atof(part.c_str()));
		}
		if(x_line.size()==0)
			continue;
		if(n_dat==0)
		{
			n_dim=x_line.size();
		}
		if(x_line.size()!=n_dim)
		{
			cerr<<"input data error"<<endl;
			breturn=1;
			break;
		}
		n_dat++;
		xx.insert(xx.end(),x_line.begin(),x_line.end());
	}

	n_dat/=n_conf;

	x=new double[n_dat*n_dim*n_conf];


	for(i=0;i<n_dat*n_dim*n_conf;i++)
		x[i]=xx.at(i);

	return breturn;
};


bool  CAnn::loadx(string name)
{
	int i;
	ifstream fin;
	string line,part;
	istringstream iss;
	vector<double> x_line;
	vector<double> xx;
	bool breturn=0;

	fin.open(name.c_str());
	n_dat=0;
	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		x_line.clear();
		while(iss>>part)
		{
			x_line.push_back(atof(part.c_str()));
		}
		if(n_dat==0)
		{
			n_dim=x_line.size();
		}
		if(x_line.size()!=n_dim)
		{
			cerr<<"input data error"<<endl;
			breturn=1;
			break;
		}
		n_dat++;
		xx.insert(xx.end(),x_line.begin(),x_line.end());
	}


	x=new double[n_dat*n_dim];


	for(i=0;i<n_dat*n_dim;i++)
		x[i]=xx.at(i);


	return breturn;
};


bool CAnn::loady(string name)
{
	int i,n_dat2;
	ifstream fin;
	string line,part;
	istringstream iss;
	vector<double> x_line;
	vector<double> xx;
	bool breturn=0;

	fin.open(name.c_str());
	n_dat2=0;
	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		x_line.clear();
		while(iss>>part)
		{
			x_line.push_back(atof(part.c_str()));
		}
		if(x_line.size()!=1)
		{
			cerr<<"input data error in Y"<<endl;
			breturn=1;
			break;
		}

		n_dat2++;
		xx.insert(xx.end(),x_line.begin(),x_line.end());
	}

	if(n_dat2!=n_dat)
	{
		breturn=1;
	}
	else
	{
		y=new double[n_dat];
		for(i=0;i<n_dat;i++)
			y[i]=xx.at(i);
	}


	return breturn;
};


void CAnn::mapminmax()
{
	double tmin,tmax,v;
	int i,j,step;


	v_min=new double[n_dim];
	v_max=new double[n_dim];

	for(i=0;i<n_dim;i++)
	{
		tmax=-100000.0;
		tmin=100000.0;
		for(j=0;j<n_dat;j++)
		{
			step=j*n_dim+i;
			v=x[step];
			if(v>tmax)
				tmax=v;
			if(v<tmin)
				tmin=v;
		}

		v_min[i]=tmin;
		v_max[i]=tmax;

		for(j=0;j<n_dat;j++)
		{
			step=j*n_dim+i;
			v=x[step];
			x[step]=(v-tmin)/(tmax-tmin)*2-1;
		}

	}

	y_max=-100000.0;
	y_min=100000.0;
	for(j=0;j<n_dat;j++)
	{
		if(y[j]>y_max)
			y_max=y[j];
		if(y[i]<y_min)
			y_min=y[j];
	}
	for(j=0;j<n_dat;j++)
	{
		y[j]=(y[j]-y_min)/(y_max-y_min)*2-1;
	}
	return;
};


void CAnn::mapminmax_md()
{
	double tmin,tmax,v;
	int i,j,step;


	v_min=new double[n_dim];
	v_max=new double[n_dim];

	for(i=0;i<n_dim;i++)
	{
		tmax=-100000.0;
		tmin=100000.0;
		for(j=0;j<n_dat*n_conf;j++)
		{
			step=j*n_dim+i;
			v=x[step];
			if(v>tmax)
				tmax=v;
			if(v<tmin)
				tmin=v;
		}

		v_min[i]=tmin;
		v_max[i]=tmax;

		for(j=0;j<n_dat*n_conf;j++)
		{
			step=j*n_dim+i;
			v=x[step];
			x[step]=(v-tmin)/(tmax-tmin)*2-1;
		}

	}

	y_max=-100000.0;
	y_min=100000.0;
	for(j=0;j<n_dat;j++)
	{
		if(y[j]>y_max)
			y_max=y[j];
		if(y[i]<y_min)
			y_min=y[j];
	}
	for(j=0;j<n_dat;j++)
	{
		y[j]=(y[j]-y_min)/(y_max-y_min)*2-1;
	}
	return;
};


void CAnn::xapplyminmax_md()
{
	double tmin,tmax;
	int i,j,step;


	for(i=0;i<n_dim;i++)
	{
		tmin=v_min[i];
		tmax=v_max[i];

		for(j=0;j<n_dat*n_conf;j++)
		{
			step=j*n_dim+i;
			x[step]=(x[step]-tmin)/(tmax-tmin)*2-1;
		}

	}

	return;
};


// New function for OpenACC
// Vector Routine
void CAnn::xapplyminmax_acc(double *xx)
{
	double tmin,tmax;
	int i,j,step;
	#pragma acc loop vector independent private(tmin,tmax,step)
	for(i = 0; i < n_dim; i++)
	{
		tmin=v_min[i];
		tmax=v_max[i];
		#pragma acc loop seq
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

void CAnn::train_it_md(double percent)
{
	lm_status_struct status;
	lm_control_struct control;
	int i,j;
	int n_dat_val;
	double *x2,*y2;
	int *n;

	n=new int[n_dat];
	x2=new double[n_dat*n_dim*n_conf];
	y2=new double[n_dat];


	randperm(n_dat,n);
	for(i=0;i<n_dat;i++)
	{
		y2[i]=y[n[i]];
		for(j=0;j<n_dim*n_conf;j++)
		{
			x2[i*n_dim*n_conf+j]=x[n[i]*n_dim*n_conf+j];
		}
	}


	n_dat_val=(int)(n_dat*percent);

	control = lm_control_double;
    control.verbosity = 3;

	data_struct_md data;
	data_struct_md data2;

	data.input=x2;
	data.n_dim=n_dim;
	data.n_conf=n_conf;
	data.n_neuron=n_neuron;
	data.y=y2;
	data.f=&(CAnn::myfunc_md);

	data2.input=x2+n_dim*n_conf*(n_dat-n_dat_val);
	data2.n_dim=n_dim;
	data2.n_conf=n_conf;
	data2.n_neuron=n_neuron;
	data2.y=y2+n_dat-n_dat_val;
	data2.f=&(CAnn::myfunc_md);


	lmmin(n_par,p,n_dat-n_dat_val,(const void*) &data,n_dat_val,(const void*) &data2,evaluation_md,&control,&status);

	delete [] x2;
	delete [] y2;

	return;
};


void CAnn::train_it(double percent)
{
	lm_status_struct status;
	lm_control_struct control;
	int i,j;
	int n_dat_val;
	double *x2,*y2;
	int *n;

	n=new int[n_dat];
	x2=new double[n_dat*n_dim];
	y2=new double[n_dat];


	randperm(n_dat,n);
	for(i=0;i<n_dat;i++)
	{
		y2[i]=y[n[i]];
		for(j=0;j<n_dim;j++)
		{
			x2[i*n_dim+j]=x[n[i]*n_dim+j];
		}
	}


	n_dat_val=(int)(n_dat*percent);

	control = lm_control_double;
    control.verbosity = 3;

	data_struct_neuron data;
	data_struct_neuron data2;

	data.input=x2;
	data.n_dim=n_dim;
	data.n_neuron=n_neuron;
	data.y=y2;
	data.f=&(CAnn::myfunc_neuron);

	data2.input=x2+n_dim*(n_dat-n_dat_val);
	data2.n_dim=n_dim;
	data2.n_neuron=n_neuron;
	data2.y=y2+n_dat-n_dat_val;
	data2.f=&(CAnn::myfunc_neuron);


	lmmin(n_par,p,n_dat-n_dat_val,(const void*) &data,n_dat_val,(const void*) &data2,evaluation_neuron,&control,&status);

	delete [] x2;
	delete [] y2;

	return;
};


void CAnn::evaluation_neuron(const double *par, int n_dat, const void *pdata, double *fvect, int *user)
{
	data_struct_neuron *d;

	d=(data_struct_neuron *)pdata;

#pragma omp parallel for
	for(int i=0;i<n_dat;i++)
	{
		int begin=i*(d->n_dim);
		fvect[i]=(d->y[i])-(d->f)(d->n_dim,d->n_neuron,&(d->input[begin]),(const double *)par,0);
	}
	return;
};


double CAnn::myfunc_neuron(int ndim, int nneuron, double *xx, const double *p, int offset)
{
	int i;
	double r2;
	r2 = 0;
	const int p2=offset+ndim*nneuron;
	const int p3=p2+nneuron;
	const int p4=p3+nneuron;

	for(int j = 0; j < nneuron; j++)
	{
		double r = 0.0;
		const int p1 = offset+j*ndim;
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

// New function for OpenACC
// Seq Routine
double CAnn::myfunc_neuron_acc(int ndim, int nneuron, double *xx, const double *p, int offset)
{
	int i;
	double r2;
	r2 = 0;
	const int p2=offset+ndim*nneuron;
	const int p3=p2+nneuron;
	const int p4=p3+nneuron;

	#pragma acc loop seq reduction(+:r2)
	for(int j = 0; j < nneuron; j++)
	{
		double r = 0.0;
		const int p1 = offset+j*ndim;
		#pragma acc loop seq reduction(+:r)
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

void CAnn::evaluation_md(const double *par, int n_dat, const void *pdata, double *fvect, int *user)
{
	data_struct_md *d;

	d=(data_struct_md *)pdata;

#pragma omp parallel for
	for(int i=0;i<n_dat;i++)
	{
		int begin=i*(d->n_dim)*(d->n_conf);
		fvect[i]=(d->y[i])-(d->f)(d->n_dim,d->n_conf,d->n_neuron,&(d->input[begin]),(const double *)par);
	}
	return;
};

double CAnn::myfunc_md(int ndim, int nconf, int nneuron, double *x, const double *p)
{
	int i,j,n;
	double r2,r;
	double sum;
	const double *p2,*p3,*p4;
	const double *p1;

	p2=p+ndim*nneuron;
	p3=p2+nneuron;
	p4=p3+nneuron;

	sum=0.0;
	for(n=0;n<nconf;n++)
	{

		r2=0;
		for(j=0;j<nneuron;j++)
		{
			p1=p+j*ndim;
			r=0.0;
			for(i=0;i<ndim;i++)
			{
				r+=x[i+n*ndim]*p1[i];
			}
			r+=p2[j];
			if(r<-30 )
				r=0.0;
			else if(r>30)
				r=1.0;
			else
				r=(1-exp(-2*r))/(1+exp(-2*r));
			r=r*p3[j];
			r2+=r;
		}
		r2+=p4[0];

		sum+=r2;
	}

	sum/=nconf;
	return sum;
};

void CAnn::savehead(string filename, int n)
{
	int i,j;
	ofstream fout;

	fout.open(filename.c_str());

	fout<<"n_dim "<<n_dim<<endl;
	fout<<"n_neuron "<<n_neuron<<endl;
	fout<<"n_par "<<n_par<<endl;

	for(i=0;i<n_dim;i++)
		fout<<v_min[i]<<" "<<v_max[i]<<endl;

	fout<<y_min<<" "<<y_max<<endl;

	fout<<n<<endl;

	fout.close();
};


void CAnn::savep(string filename,int n)
{
	int i,j;
	ofstream fout;

	fout.open(filename.c_str(),std::ofstream::app);
	for(j=0;j<n_par;j++)
		fout << p_save_flat[n*n_par + j] << " ";
	fout<<endl;

	fout.close();
}


void CAnn::save(string filename)
{
	int i,j;
	ofstream fout;

	fout.open(filename.c_str());

	fout<<"n_dim "<<n_dim<<endl;
	fout<<"n_neuron "<<n_neuron<<endl;
	fout<<"n_par "<<n_par<<endl;

	for(i=0;i<n_dim;i++)
		fout<<v_min[i]<<" "<<v_max[i]<<endl;

	fout<<y_min<<" "<<y_max<<endl;

	fout << p_save_size << endl;

	for(i=0;i<p_save_size;i++)
	{
		for(j=0;j<n_par;j++)
			fout << p_save_flat[i*n_par + j] << " ";
		fout<<endl;
	}

	fout.close();
};


// Updated function for OpenACC
// Contains data directives
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
	p_save_flat = new double[n_set*n_par]; // p_save_flat replaces a 2D vector
	p_save_size = n_set;
	for(j=0;j<n_set;j++)
	{
		for(i=0;i<n_par;i++)
		{
			p_save_flat[j*n_par + i] = *pdata;
			pdata++;
		}
	}

	#pragma acc enter data copyin(this)
	#pragma acc enter data copyin(p_save_flat[0:n_set*n_par])
	#pragma acc enter data copyin(v_min[0:n_dim],v_max[0:n_dim])
};



void CAnn::load(string filename)
{
	int i,j;
	ifstream fin;
	string part;
	int n_set;
	vector<double> t;

	fin.open(filename.c_str());

	fin>>part>>n_dim;
	fin>>part>>n_neuron;
	fin>>part>>n_par;

	n_par=n_dim*n_neuron+n_neuron+n_neuron+1;

	p=new double[n_par];
	v_min=new double[n_dim];
	v_max=new double[n_dim];

	for(i=0;i<n_dim;i++)
		fin>>v_min[i]>>v_max[i];

	fin>>y_min>>y_max;

	fin>>n_set;

	p_save_flat = new double[n_set*n_par];
	p_save_size = n_set;
	for(j=0;j<n_set;j++)
	{
		t.clear();
		for(i=0;i<n_par;i++)
		{
			fin>>p[i];
			p_save_flat[j*n_par + i] = p[i];
		}
	}

	fin.close();
};

double CAnn::assess(string ann_name,string x_name,string y_name)
{
	int i;
	vector<double> out;
	double rms,tt;
	vector<vector<double> > xx;


	out=predict(0,ann_name,x_name,xx);


	if(loady(y_name)!=0)
	{
		cout<<"Load Y error"<<endl;
		return 1000.0;
	}


	rms=0.0;
	for(i=0;i<n_dat;i++)
	{
		tt=out.at(i)-y[i];
		rms+=tt*tt;
	}
	rms=sqrt(rms/n_dat);

	return rms;
};

double CAnn::assess_md(string ann_name,string x_name,string y_name)
{
	int i;
	vector<double> out;
	double rms,tt;
	vector<vector<double> > xx;


	out=predict_md(0,ann_name,x_name,xx);


	if(loady(y_name)!=0)
	{
		cout<<"Load Y error"<<endl;
		return 1000.0;
	}


	rms=0.0;
	for(i=0;i<n_dat;i++)
	{
		tt=out.at(i)-y[i];
		rms+=tt*tt;
	}
	rms=sqrt(rms/n_dat);

	return rms;
};


double CAnn::predict_one( double *xx, int vec_size )
{
	double tt,out;

	if(vec_size!=n_dim)
		return 0;

	n_dat=1;
	xapplyminmax_acc(xx);
	out=0;

	for(int j=0;j<p_save_size;j++)
	{
		tt=CAnn::myfunc_neuron(n_dim, n_neuron, xx, p_save_flat, n_par*j);
		tt=(tt+1)/2*(y_max-y_min)+y_min;
		out+=tt;
	}

	out/=p_save_size;
	return out;
};


// New function for OpenACC
// Vector Routine
double CAnn::predict_one_acc( double *xx, int vec_size )
{
	double tt,out;

	if(vec_size!=n_dim)
		return 0;

	n_dat=1;
	xapplyminmax_acc(xx);
	out=0;

	#pragma acc loop vector independent reduction(+:out) private(tt)
	for(int j=0;j<p_save_size;j++)
	{
		tt=CAnn::myfunc_neuron_acc(n_dim, n_neuron, xx, p_save_flat, n_par*j);
		tt=(tt+1)/2*(y_max-y_min)+y_min;
		out+=tt;
	}

	out/=p_save_size;
	return out;
};


double CAnn::predict_one_md(int n, vector<double> xx )
{
	int j;
	double tt,out;

	n_conf=n;

	if(xx.size()!=n_dim*n_conf)
		return 0;

	n_dat=1;
	x=new double[n_conf*n_dim];
	for(j=0;j<n_dim*n_conf;j++)
	{
		x[j]=xx.at(j);
	}
	xapplyminmax_md();

	out=0;
	for(int j=0;j<(int)p_save_size;j++)
	{
		for(int i=0;i<n_par;i++)
			p[i]=p_save_flat[j*n_par + i];

		tt=CAnn::myfunc_md(n_dim,n_conf,n_neuron,x,p);
		tt=(tt+1)/2*(y_max-y_min)+y_min;
		out+=tt;
	}
	out/=p_save_size;
	return out;
};


vector<double> CAnn::predict(int flag, string ann_name,string x_name, vector< vector<double> > xx)
{
	int t;
	vector<double> out;
	double tt;
	double *pre;

	out.clear();


	load(ann_name);
	t=n_dim;

	if(flag==0)
	{
		if(loadx(x_name)!=0)
		{
			cout<<"load X fail"<<endl;
			return out;
		}
	}
	else //flag==1
	{
		loadx(xx);
	}

	if(t!=n_dim)
	{
		cout<<"Inconsistent dimension between trained ANN and input X data"<<endl;
		return out;
	}

	xapplyminmax();
	pre=new double[n_dat];
	for(int i=0;i<n_dat;i++)
		pre[i]=0.0;
	for(int j=0;j<(int)p_save_size;j++)
	{
		for(int i=0;i<n_par;i++)
			p[i] = p_save_flat[j*n_par + i];

#pragma omp parallel for
		for(int i=0;i<n_dat;i++)
		{
			int begin=i*n_dim;
			tt=CAnn::myfunc_neuron(n_dim,n_neuron,x+begin,p,0);
			tt=(tt+1)/2*(y_max-y_min)+y_min;
			pre[i]+=tt;
		}
	}

	for(int i=0;i<n_dat;i++)
	{
		pre[i]/=p_save_size;
		out.push_back(pre[i]);
	}

	return out;
}


vector<double> CAnn::predict_md(int flag, string ann_name,string x_name, vector< vector<double> > xx)
{
	int t;
	vector<double> out;
	double tt;
	double *pre;

	out.clear();


	load(ann_name);
	t=n_dim;

	if(flag==0)
	{
		if(loadx_md(x_name)!=0)
		{
			cout<<"load X fail"<<endl;
			return out;
		}
	}
	else //flag==1
	{
		loadx(xx);
	}

	if(t!=n_dim)
	{
		cout<<"Inconsistent dimension between trained ANN and input X data"<<endl;
		return out;
	}

	xapplyminmax_md();
	pre=new double[n_dat];
	for(int i=0;i<n_dat;i++)
		pre[i]=0.0;
	for(int j=0;j<(int)p_save_size;j++)
	{
		for(int i=0;i<n_par;i++)
			p[i] = p_save_flat[j*n_par + i];

#pragma omp parallel for
		for(int i=0;i<n_dat;i++)
		{
			int begin=i*n_dim*n_conf;
			tt=CAnn::myfunc_md(n_dim,n_conf,n_neuron,x+begin,p);
			tt=(tt+1)/2*(y_max-y_min)+y_min;
			pre[i]+=tt;
		}
	}

	for(int i=0;i<n_dat;i++)
	{
		pre[i]/=p_save_size;
		out.push_back(pre[i]);
	}

	return out;
}




double CAnn::myfunc_mix(int n_neuron, double *x, const double *p)
{
	int i;
	const double *p1,*p2,*p3,*p4;
	double *x1,*x2,*x3,*x4;
	double r1,r2,r3;

	r1=0;
	p1=p;
	x1=x;
	for(i=0;i<500;i++)
		r1+=x1[i]*p1[i];

	r2=0;
	p2=p1+500;
	x2=x+500;
	for(i=0;i<260;i++)
		r2+=x2[i]*p2[i];

	r3=0;
	p3=p1+500+260;
	x3=x+500+260;
	for(i=0;i<260;i++)
		r3+=x3[i]*p3[i];

	p4=p+500+260+260;
	x4=x+500+260+260;

	x4[0]=r1/11.464-1.7415;     /* rescale and shift r1 to be within [-1,1], as required by next step*/
	x4[1]=r2/1.4956-11.8745;
	x4[2]=r3/1.5156-12.031;

	return myfunc_neuron(10,n_neuron,x4,p4,0);
};

void CAnn::evaluation_mix(const double *par, int m_dat, const void *pdata, double *fvect, int *user)
{
	data_struct_mix *d;

	d=(data_struct_mix *)pdata;

#pragma omp parallel for
	for(int i=0;i<m_dat;i++)
	{
		int begin=i*(d->n_dim);
		fvect[i]=(d->y[i])-(d->f)(d->n_neuron,&(d->input[begin]),(const double *)par);
	}
	return;
};
