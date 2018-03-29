#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <sstream>
using namespace std;

#include "ann.h"
#include "mainbody.h"

#include <openacc.h>

static int blosum[400]=
{4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,
-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,
-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,
-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,
0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,
-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,
-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,
0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,
-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,
-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,
-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,
-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,
-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,
-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,
-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,
1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,
0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,
-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,
-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,
0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4};

// Seperate from dihe_process
// replace out with static sized array of 32
// These are the 3 functions needed to completely seperate GPU kernel from dihe_process

//int MAX = 0;

#pragma acc routine seq
int * pos(int in, int *num_arr, int num_size, int &index_size)
{
	//cout << "Pos" << endl;
	int i;
	int base,stop;
	//vector<int> out;
	//int *out_pos; // seems to be the max value, unconfirmed
	int out_pos[7];


	if(in==1)
		base=1;
	else
		base=num_arr[in-2]+1;
	
	stop=num_arr[in-1];

	if(in<=0)
		base=num_arr[num_size-1]+100;
	if(in>num_size)
		base=num_arr[num_size-1]+100;

	index_size = stop-base+1;
	//out_pos = new int[index_size];

	for(i=base;i<=stop;i++)
	{
		out_pos[i-base]=i;//.push_back(i);
		//cout<<"in is "<<in<<" "<<i<<endl;
	}

	return out_pos;
}

#pragma acc routine seq
int process_static_new(int id,int cut,int order,int order2, double *out_arr, int *num_arr, int num_size, int ndihe, int nframe, double *dihe, int out_index)
{
	//cout << "process_static_new" << endl;
	int p,i,j,k,t,base,stop;
	int index[7];
	int index_size;
	double phi;
	double cosphi[10];
	double sinphi[10];


//	index=pos(id, num_arr, num_size, index_size);

	if(id==1)
		base=1;
	else
		base=num_arr[id-2]+1;
	
	stop=num_arr[id-1];

	if(id<=0)
		base=num_arr[num_size-1]+100;
	if(id>num_size)
		base=num_arr[num_size-1]+100;

	index_size = stop-base+1;

	for(i=base,j=0;i<=stop;i++,j++)
	{
		index[j]=i;//.push_back(i);
		if(j>=7)
			cout << "ERROR" << endl;
		//cout<<"in is "<<in<<" "<<i<<endl;
	}


	if(index_size>=cut)
		t=cut;
	else
		t=index_size;

	for(i=0;i<2;i++)
	{
		p=index[i];
		for(k=1;k<=order;k++)
			cosphi[k]=sinphi[k]=0.0;
		for(j=0;j<nframe;j++)
		{
			base=j*ndihe;
			phi=dihe[base+p-1];
			for(k=1;k<=order;k++)
			{
				cosphi[k]+=cos(phi*k);
				sinphi[k]+=sin(phi*k);
			}
		}
		for(k=1;k<=order;k++)
		{
			out_arr[out_index++]=cosphi[k]/nframe;//.push_back(cosphi[k]/nframe);
			out_arr[out_index++]=sinphi[k]/nframe;//.push_back(sinphi[k]/nframe);
		}
	}

	for(i=2;i<t;i++)
        {
                p=index[i];
                for(k=1;k<=order2;k++)
                        cosphi[k]=sinphi[k]=0.0;
                for(j=0;j<nframe;j++)
                {
                        base=j*ndihe;
                        phi=dihe[base+p-1];
                        for(k=1;k<=order2;k++)
                        {
                                cosphi[k]+=cos(phi*k);
                                sinphi[k]+=sin(phi*k);
                        }
                }
                for(k=1;k<=order2;k++)
                {
                        out_arr[out_index++]=cosphi[k]/nframe;//.push_back(cosphi[k]/nframe);
			out_arr[out_index++]=sinphi[k]/nframe;//.push_back(sinphi[k]/nframe);
                }
        }

	for(i=t;i<cut;i++)
	{
		for(j=0;j<order2*2;j++)
			out_arr[out_index++]=0.0;//.push_back(0.0);
	}
	return out_index;
}

#pragma acc routine seq
void ca_ann(int id, double *out_arr, int ndihe, int nframe, double *dihe, int *num_arr, int num_size)
{
	//cout << "ca_ann" << endl;
	//out.clear();
	int out_index = 0;
	out_index = process_static_new(id-1,4,1,1, out_arr, num_arr, num_size, ndihe, nframe, dihe, out_index);
	out_index = process_static_new(id  ,4,2,2, out_arr, num_arr, num_size, ndihe, nframe, dihe, out_index);
	out_index = process_static_new(id+1,4,1,1, out_arr, num_arr, num_size, ndihe, nframe, dihe, out_index);
  	return;
}


/*char code(CAminoacid *v, int in)
{
	char c;
	in=in-1;
	if(in<0 || in >(int)v.size()-1)
		c='X';
	else
		c=v[in]->OneLetterName;

	return c;
}*/




//class Dihe_process, used by ppm only


void CDihe_process::init(vector<int> innum,vector<double> *indihe)
{
	double st = omp_get_wtime();
	num=innum;
	dihe=indihe;
	if(num.size()>2)
	{
		ndihe=num.at(num.size()-1);
		nframe=dihe->size()/ndihe;
	}
	else
	{
		ndihe=0;
		nframe=0;
	}
	cout << "dihe_process:mainbody::init: " << omp_get_wtime()-st << " seconds" << endl;
	return;
}

void CDihe_process::init(vector<int> innum,vector<double> *indihe, vector<dihe_group> *indihe_index)
{
	num=innum;
	dihe=indihe;
	dihe_index=indihe_index;

	if(num.size()>2)
	{
		ndihe=num.at(num.size()-1);
		nframe=dihe->size()/ndihe;
	}
	else
	{
		ndihe=0;
		nframe=0;
	}

	return;
}

void CDihe_process::init(vector<int> innum,vector<double> *indihe, vector<double> *inangle)
{
	num=innum;
	dihe=indihe;
	angle=inangle;
	if(num.size()>2)
	{
		ndihe=num.at(num.size()-1);
		nframe=dihe->size()/ndihe;
	}
	else
	{
		ndihe=0;
		nframe=0;
	}

	return;
}


void CDihe_process::hb_expand(int type)
{
	int i;
	int begin,stop;
	vector<double> temp;

	begin=table[type]*8;
	stop=table[type]*8+8;

	if(stop==0)
	{
		out.clear();
		out.resize(18*8,0.0);
	}
	else
	{
		temp.resize(18*8,0.0);
		for(i=begin;i<stop;i++)
			temp.at(i)=out.at(i-begin);
		out=temp;
	}

	return;
}

void CDihe_process::hb_expand2(int type)
{
    int i;
    int begin,stop;
    vector< vector<double> > temp;
    vector<double> t;
    
    t.resize(nframe,0.0);
    
    begin=table[type]*8;
    stop=table[type]*8+8;
    
    if(stop==0)
    {
        out2.clear();
        out2.resize(18*8,t);
    }
    else
    {
        temp.resize(18*8,t);
        for(i=begin;i<stop;i++)
            temp.at(i)=out2.at(i-begin);
        out2=temp;
    }
    
    return;
}


vector<int> CDihe_process::pos(int in)
{
	int i;
	int base,stop;
	vector<int> out;


	if(in==1)
		base=1;
	else
		base=num.at(in-2)+1;
	
	stop=num.at(in-1);

	if(in<=0)
		base=num.at(num.size()-1)+100;
	if(in>(int)num.size())
		base=num.at(num.size()-1)+100;

	
	for(i=base;i<=stop;i++)
	{
		out.push_back(i);
		//cout<<"in is "<<in<<" "<<i<<endl;
	}

	return out;
}


vector<int> CDihe_process::pos_angle(int in)
{
	int i;
	int base,stop;
	vector<int> out;


	if(in==1)
		base=1;
	else
		base=num.at(in-2)*2+1;
	
	stop=num.at(in-1)*2;

	if(in<=0)
		base=num.at(num.size()-1)+100;
	if(in>(int)num.size())
		base=num.at(num.size()-1)+100;

	
	for(i=base;i<=stop;i++)
	{
		out.push_back(i);
		//cout<<"in is "<<in<<" "<<i<<endl;
	}

	return out;
}

void CDihe_process::allproton(int id)
{
	out.clear();
	process(id,2,2,2);
	return;
}

void CDihe_process::allproton2(int id)
{
    out2.clear();
    process2(id,2,2,2);
    return;
}

bool CDihe_process::test_good(int id,int cut)
{
	bool r;
	vector<int> index;
	int t,i;

	index=pos(id);

	if((int)index.size()>=cut)
		t=cut;
	else
		t=index.size();

	r=1;
	for(i=0;i<t;i++)
	{
		if(dihe_index->at(index.at(i)-1).bgood==0)
			r=0;
	}

	return r;
}


bool CDihe_process::test(int id, int t1, int t2)
{
	bool r;
	if(test_good(id-1,t2) && test_good(id,t1) && test_good(id+1,t2) )
		r=1;
	else
		r=0;
	return r;
}


bool CDihe_process::test_proton(int id, int type)
{
	return ((table[type]==-1) || test_good(id,2));
}


bool CDihe_process::ca_static_new(int id)
{
	out.clear();
	bool r;

	if(test_good(id-1,3) && test_good(id,4) && test_good(id+1,3) )
	{
		process_static_new(id-1,3,2,2);
		process_static_new(id  ,4,3,3);
		process_static_new(id+1,3,2,2);
		r=1;
	}
	else
		r=0;

	return r;
	
}

void CDihe_process::md_new(int id)
{
	out.clear();

	process_static_new(id-1,4,2,2);
	process_static_new(id  ,4,3,3);
	process_static_new(id+1,4,2,2);
  	return;
}

void CDihe_process::md_new_detail(int id, int n)
{
	out.clear();
    
	process_md_sep(n,id-1,4,2,2);
	process_md_sep(n,id  ,4,3,3);
	process_md_sep(n,id+1,4,2,2);
  	return;
}


void CDihe_process::ca(int id)
{
	out.clear();

	process(id-1,3,3,3);
	process(id  ,4,3,3);
	process(id+1,3,3,3);
  	return;
}



void CDihe_process::ca_ann(int id)
{
	out.clear();

	process_static_new(id-1,4,1,1);
	process_static_new(id  ,4,2,2);
	process_static_new(id+1,4,1,1);
  	return;
}

void CDihe_process::md_ann(int id,int n)
{
	out.clear();
	process_md_sep(n,id-1,4,1,1);
	process_md_sep(n,id  ,4,1,1);
	process_md_sep(n,id+1,4,1,1);
	return;
}


void CDihe_process::process_md_sep(int nn,int id,int cut, int order1, int order2)
{
	int i,ii,j,p;
	int base;
	int t;
	double phi;
	int order;
	vector<int> index;
	vector<double> temp;

	index=pos(id);
	if((int)index.size()>=cut)
		t=cut;
	else
		t=index.size();


	for(i=0;i<t;i++)
	{
		p=index.at(i);
		order=order1;
		if(i>1) order=order2;

		for(ii=1;ii<=order;ii++)
		{
			j=nn;
			base=j*ndihe;
			phi=dihe->at(base+p-1);
			out.push_back(cos(phi*ii));
			
			j=nn;
			base=j*ndihe;
			phi=dihe->at(base+p-1);
			out.push_back(sin(phi*ii));
		}
	}

	for(i=t;i<cut;i++)
	{
		for(ii=1;ii<=order2;ii++)
		{
			out.push_back(0.0);
			out.push_back(0.0);
		}
	}

	return;
}

void CDihe_process::process(int id,int cut,int order,int order2)
{
	int p,i,j,k,t,base;
	vector<int> index;
	double phi;
	double cosphi[10];
	double sinphi[10];

	

	index=pos(id);
	if((int)index.size()>=cut)
		t=cut;
	else
		t=index.size();

	for(i=0;i<2;i++)
	{
		p=index.at(i);
		for(k=1;k<=order;k++)
			cosphi[k]=sinphi[k]=0.0;
		for(j=0;j<nframe;j++)
		{
			base=j*ndihe;
			phi=dihe->at(base+p-1);
			for(k=1;k<=order;k++)
			{
				cosphi[k]+=cos(phi*k);
				sinphi[k]+=sin(phi*k);
			}
		}
		for(k=1;k<=order;k++)
		{
			out.push_back(cosphi[k]/nframe);
		}
		for(k=1;k<=order;k++)
		{
			out.push_back(sinphi[k]/nframe);
		}
	}

	for(i=2;i<t;i++)
        {
                p=index.at(i);
                for(k=1;k<=order2;k++)
                        cosphi[k]=sinphi[k]=0.0;
                for(j=0;j<nframe;j++)
                {
                        base=j*ndihe;
                        phi=dihe->at(base+p-1);
                        for(k=1;k<=order2;k++)
                        {
                                cosphi[k]+=cos(phi*k);
                                sinphi[k]+=sin(phi*k);
                        }
                }
                for(k=1;k<=order2;k++)
                {
                        out.push_back(cosphi[k]/nframe);
                }
                for(k=1;k<=order2;k++)
                {
                        out.push_back(sinphi[k]/nframe);
                }
        }

	for(i=t;i<cut;i++)
	{
		for(j=0;j<order2*2;j++)
			out.push_back(0.0);
	}
}


void CDihe_process::process_static_new(int id,int cut,int order,int order2)
{
	int p,i,j,k,t,base;
	vector<int> index;
	double phi;
	double cosphi[10];
	double sinphi[10];


	index=pos(id);
	if((int)index.size()>=cut)
		t=cut;
	else
		t=index.size();

	for(i=0;i<2;i++)
	{
		p=index.at(i);
		for(k=1;k<=order;k++)
			cosphi[k]=sinphi[k]=0.0;
		for(j=0;j<nframe;j++)
		{
			base=j*ndihe;
			phi=dihe->at(base+p-1);
			for(k=1;k<=order;k++)
			{
				cosphi[k]+=cos(phi*k);
				sinphi[k]+=sin(phi*k);
			}
		}
		for(k=1;k<=order;k++)
		{
			out.push_back(cosphi[k]/nframe);
			out.push_back(sinphi[k]/nframe);
		}
	}

	for(i=2;i<t;i++)
        {
                p=index.at(i);
                for(k=1;k<=order2;k++)
                        cosphi[k]=sinphi[k]=0.0;
                for(j=0;j<nframe;j++)
                {
                        base=j*ndihe;
                        phi=dihe->at(base+p-1);
                        for(k=1;k<=order2;k++)
                        {
                                cosphi[k]+=cos(phi*k);
                                sinphi[k]+=sin(phi*k);
                        }
                }
                for(k=1;k<=order2;k++)
                {
                        out.push_back(cosphi[k]/nframe);
						out.push_back(sinphi[k]/nframe);
                }
        }

	for(i=t;i<cut;i++)
	{
		for(j=0;j<order2*2;j++)
			out.push_back(0.0);
	}
}

void CDihe_process::proton(int id)
{
	out2.clear();
	process_fit(id,7);
  	return;
}


void CDihe_process::ca2(int id)
{
	out2.clear();

	process2(id-1,3,3,3);
	process2(id  ,4,3,3);
	process2(id+1,3,3,3);
  	return;
}

void CDihe_process::for_fit(int id)
{
	out2.clear();

	process_fit(id-1,4);
	process_fit(id,4);
	process_fit(id+1,4);
	return;
}


void CDihe_process::fit_angle(int id)
{
	int i,j,p,base;
	vector<int> index;
	vector<double> temp;
	double phi;

	index=pos_angle(id);

	out2.clear();

	int touse[5]={0,1,3,4,5};

	
	for(i=0;i<5;i++)
	{
		if(touse[i]<(int)index.size())
		{
			p=index.at(touse[i]);
			temp.clear();
			for(j=0;j<nframe;j++)
			{
				base=j*ndihe*2;
				phi=angle->at(base+p-1);
				temp.push_back(phi);
			}
			out2.push_back(temp);
		}
		else
		{
			temp.clear();
			for(j=0;j<nframe;j++)
			{
				temp.push_back(0.0);
			}
			out2.push_back(temp);
		}

	}

}


void CDihe_process::process_fit(int id,int cut)
{
	int i,j,p;
	int base;
	int t;
	double phi;
	vector<int> index;
	vector<double> temp;

	index=pos(id);
	if((int)index.size()>=cut)
		t=cut;
	else
		t=index.size();


	for(i=0;i<t;i++)
	{
		p=index.at(i);
		temp.clear();
		for(j=0;j<nframe;j++)
		{
			base=j*ndihe;
			phi=dihe->at(base+p-1);
			temp.push_back(phi);
		}
		out2.push_back(temp);
	}

	for(i=t;i<cut;i++)
	{
		temp.clear();
		for(j=0;j<nframe;j++)
		{
			temp.push_back(0.0);
		}
		out2.push_back(temp);
	}

	return;
}





void CDihe_process::process2(int id,int cut,int order,int order2)
{
	int p,i,j,k,t,base;
	vector<int> index;
	double phi;
	vector< vector<double> > t1,t2;

	

	index=pos(id);
	if((int)index.size()>=cut)
		t=cut;
	else
		t=index.size();

	for(i=0;i<2;i++)
	{
		p=index.at(i);
		t1.clear();
		t2.clear();
		t1.resize(order);
		t2.resize(order);
		for(j=0;j<nframe;j++)
		{
			base=j*ndihe;
			phi=dihe->at(base+p-1);
			for(k=1;k<=order;k++)
			{
				t1.at(k-1).push_back(cos(phi*k));
				t2.at(k-1).push_back(sin(phi*k));
			}
		}
		for(k=1;k<=order;k++)
		{
			out2.push_back(t1.at(k-1));
		}
		for(k=1;k<=order;k++)
		{
			out2.push_back(t2.at(k-1));
		}
	}

	for(i=2;i<t;i++)
    {
		p=index.at(i);
		t1.clear();
		t2.clear();
		t1.resize(order2);
		t2.resize(order2);
		for(j=0;j<nframe;j++)
        {
			base=j*ndihe;
            phi=dihe->at(base+p-1);
            for(k=1;k<=order2;k++)
            {
				t1.at(k-1).push_back(cos(phi*k));
				t2.at(k-1).push_back(sin(phi*k));
			}
		}
        for(k=1;k<=order2;k++)
		{
			out2.push_back(t1.at(k-1));
		}
        for(k=1;k<=order2;k++)
        {
			out2.push_back(t2.at(k-1));
		}
    }

	for(i=t;i<cut;i++)
	{
		t1.clear();
		t2.clear();
		t1.resize(order2);
		t2.resize(order2);
		for(j=0;j<nframe;j++)
		{
			for(k=1;k<=order2;k++)
            {
				t1.at(k-1).push_back(0.0);
				t2.at(k-1).push_back(0.0);
			}
		}
		for(k=1;k<=order2;k++)
		{
			out2.push_back(t1.at(k-1));
		}
        for(k=1;k<=order2;k++)
        {
			out2.push_back(t2.at(k-1));
		}
	}
}

vector<double> CDihe_process::output()
{
	return out;
}

vector< vector<double> > CDihe_process::output2()
{
	return out2;
}

CDihe_process::CDihe_process(void)
{
	int i;
	table=new int[98+1];
	for(i=0;i<98+1;i++)
		table[i]=-1;

	for(i=0;i<18;i++)
		table[hbs[i]]=i;

};


CDihe_process::~CDihe_process(void)
{
	delete [] table ;
};



CMainbody::CMainbody()
{
	int i;

	bnew=0;

	sep_table=new int[98+1];
	for(i=0;i<98+1;i++)
		sep_table[i]=-1;

	for(i=0;i<19;i++)
		sep_table[sep[i]]=i;
};

CMainbody::~CMainbody()
{
	if(bnew)
	{
		delete pdb;
		delete traj;
	}
	delete [] sep_table;

#pragma acc exit data delete(anistropy_new)
#pragma acc exit data delete(allprotons3_new)
#pragma acc exit data delete(ring_index_new)
#pragma acc exit data delete(hbond_arr)
};



int CMainbody::loadpdb(CPdb *p_pdb, CTraj * p_traj)
{
	pdb=p_pdb;
	traj=p_traj;

	natom=pdb->getnatom();
	nres=pdb->getnres();
	nconf=traj->getnframe();

	return nconf;

}



int CMainbody::loadpdb(string name)
{
	pdb=new CPdb;
	traj=new CTraj;
	bnew=1;

	natom=pdb->loadpdb(name);
	nres=pdb->getnres();
	traj->setnres(nres);
	traj->setnatom(natom);
	nconf=traj->loadcoor(name);
	return nconf;
}


int CMainbody::loadpdb(string name,string name2)
{
	pdb=new CPdb;
	traj=new CTraj;
	bnew=1;

	pdb_name=name;

	natom=pdb->loadpdb(name);
	nres=pdb->getnres();
	traj->setnres(nres);
	traj->setnatom(natom);

		
    nconf=traj->loadcoor(name);
 

	return nconf;
}



// Updated function for OpenACC
// Includes data directives and GPU function call
void CMainbody::load(string bmrbname)
{

	bmrb.process(bmrbname.c_str()); // Read from file
	pdb->attach_bmrb(bmrb);


	pdb->getdihe(&dihe_index,&dihe_num);
	pdb->getring(&ring_index);
	// getring ///
	ring_index_new = ring_index.data();
	ring_index_size = ring_index.size();
#pragma acc enter data copyin(ring_index_new[0:ring_index_size])
	//////////////
	pdb->ani(&anistropy);
	anistropy_new = anistropy.data();
	anistropy_size = anistropy.size();
#pragma acc enter data copyin(anistropy_new[0:anistropy_size])
	//pdb->proton(&protons);
	pdb->proton_acc(&protons);
	pdb->allproton(&allprotons);
	pdb->process_ambig(2);
	pdb->allproton3(&allprotons3);
	heavy=pdb->getheavy();
	pdb->getbb(&bb);
	pdb->bbnh(&bbnh);
	pdb->bbhbond(&hbond);
	pdb->schbond(&hbond);  //This is new ! 
	hbond_arr = hbond.data();
	hbond_size = hbond.size();
#pragma acc enter data copyin(hbond_arr[0:hbond_size])

	

	ndihe=dihe_index.size();


	traj->getdihe(&dihe_index,&dihe);
	dihe_process.init(dihe_num,&dihe,&dihe_index);
	//process bb to remove all entry that has missing part !!
	//bbnh willn't take effect if bb is not there for particular residue
	clear(bb);
	allprotons = clear_filter(allprotons);
	allprotons3 = clear_filter(allprotons3);
	allprotons3_new = allprotons3.data();
	allprotons3_size = allprotons3.size();
#pragma acc enter data copyin(allprotons3_new[0:allprotons3_size])

	//seperate ring group to two, one for internal, one for surface, according to contact sum !
	int i;
	vector<int> ring_atom;
	//vector<float> result;
	ring_atom.clear();
	//result.clear();
	for(i=0;i<(int)ring_index.size();i++)
	{
		ring_atom.push_back(ring_index.at(i).x2);
	}
	vector<float> result(ring_atom.size());
	traj->get_contact(1.00,0.0,ring_atom,heavy,&result);
	ring_index_internal.clear();
	ring_index_external.clear();

	for(i=0;i<(int)ring_index.size();i++)
	{
		if(result.at(i)>2.5)
			ring_index_internal.push_back(ring_index.at(i));
		else
			ring_index_external.push_back(ring_index.at(i));
	}

	return;
}


// New function for OpenACC
// Sequential, but optimized
// Functionally different, but produces the same results
vector<proton> CMainbody::clear_filter(vector<proton> protons){
	double st = omp_get_wtime();
	vector<proton> results;
	results.reserve(protons.size());
    	for (auto& p : protons)
		if((p.id>=1 || p.id<=pdb->getnres()) || (dihe_process.test_proton(p.id,p.type)==0))
            		results.push_back(p);
	return results;
}


void CMainbody::clear(vector<struct proton> &protons)
{
	double st = omp_get_wtime();
	int i;
	int id,type;

	for(i=protons.size()-1;i>=0;i--)
	{
		id=protons.at(i).id;
		type=protons.at(i).type;

		if(id<1 || id>pdb->getnres())
		{
			protons.erase(protons.begin()+i);
			continue;
		}
		if(dihe_process.test_proton(id,type)==0) //==0 means missing dihedral angles in this calculation !!
		{
			protons.erase(protons.begin()+i);
			continue;
		}
	}
	cout << "mainbody::clear(proton): " << omp_get_wtime()-st << " seconds" << endl;
}


void CMainbody::clear(vector<struct bb_group> &bb)
{
	double st = omp_get_wtime();
	int i;
	int id;
	char code,code_pre,code_fol;

	
	for(i=bb.size()-1;i>=0;i--)
	{
		//cout<<i<<endl;
		id=bb.at(i).id;

		//first and last residue are excluded
		if(id<=1)
		{
			bb.erase(bb.begin()+i);
			continue;
		}
		if(id>=pdb->getnres())
		{
			bb.erase(bb.begin()+i);
			continue;
		}

		code=pdb->code(id);
		code_pre=pdb->code(id-1);
		code_fol=pdb->code(id+1);

		//previous or following residue actually belong to another chain. 
		if(pdb->chain(id)!=pdb->chain(id-1) || pdb->chain(id)!=pdb->chain(id+1))
		{
			bb.erase(bb.begin()+i);
			continue;
		}

		//missing or unknow residue should NOT be predicted. This is also true if either previous or following residue is missing (or unknown)
		if(code_pre=='X' || code_pre=='B' || code_fol=='X' || code_fol=='B'|| code=='X' || code=='B')
		{
			bb.erase(bb.begin()+i);
			continue;
		}

		if(dihe_process.test(id,4,4)==0) //==0 means missing dihedral angles in this calculation !!
		{
			bb.erase(bb.begin()+i);
			continue;
		}

		if(bb.at(i).capos<0 ||  bb.at(i).copos<0 || bb.at(i).npos<0 )
		{
			bb.erase(bb.begin()+i);
			continue;
		}

		if(bb.at(i).cbpos<0 && bb.at(i).code!='G')
		{
			bb.erase(bb.begin()+i);
			continue;
		}

		if(bb.at(i).hpos<0 && bb.at(i).code!='P')
		{
			bb.erase(bb.begin()+i);
			continue;
		}

	}
	cout << "mainbody::clear(bb): " << omp_get_wtime()-st << " seconds" << endl;
	return;
}



int CMainbody::set_range(int begin,int stop)
{
	nconf=traj->set_range(begin,stop);
	return nconf;
}





void CMainbody::predict_bb()
{
	int i,j,j2,jj;
	int id;
	char code;
	vector<double> out;
	vector<double> in,in2;
	vector<struct double_five> ring_effect;
	vector<struct ehbond> hbond_effect;
	vector<struct double_four> ani_effect;
	vector<struct index_two> index;
	double pre[6];
	double temp;
	vector<double> eca,ecb,eco,eh,en;
	

	traj->gethbond(&hbond,&hbond_effect);
	traj->getani(&anistropy,&bbnh,&ani_effect);
	traj->getring(&ring_index,&bbnh,&ring_effect);

	
	/*carbon=bb_ca; exp.loadexp_bb("exp_ca.dat",&bb,carbon);
	carbon=bb_cb; exp.loadexp_bb("exp_cb.dat",&bb,carbon);
	carbon=bb_co; exp.loadexp_bb("exp_co.dat",&bb,carbon);
	exp.loadexp_bbn("exp_n.dat",&bbnh);
	exp.loadexp_bbnh("exp_hn.dat",&bbnh);*/


	index.resize(pdb->getnres());
	for(i=0;i<(int)index.size();i++)
		index.at(i).x1=index.at(i).x2=-1;
	for(i=0;i<(int)bb.size();i++)
		index.at(bb.at(i).id-1).x1=i+1;
	for(i=0;i<(int)bbnh.size();i++)
		index.at(bbnh.at(i).id-1).x2=i+1;
	

	for(i=0+1;i<(int)index.size()-1;i++)
	{
		//cout<<i<<endl;
		if(index.at(i).x1<0)
			continue;
		id=i+1;
		code=pdb->code(id);
		for(jj=0;jj<5;jj++)
			pre[jj]=0.0;
		//ca=cb=co=h=n=0.0;
 
		dihe_process.ca(id); 
		out=dihe_process.output();
		for(j=0;j<30;j++)
		{
			for(jj=0;jj<5;jj++)
				pre[jj]+=c_c[jj][j]*out.at(j);
		}

		for(j=42;j<60;j++)
		{
			for(jj=0;jj<5;jj++)
				pre[jj]+=c_c[jj][j-12]*out.at(j);
		}
			

		in.clear();
		for(j=0;j<12;j++)
		{
			in.push_back(out.at(j+30));
		}
		in2=Sequence::expand(code,&in);
		for(j=0;j<(int)in2.size();j++)
		{
			for(jj=0;jj<5;jj++)
				pre[jj]+=c_c[jj][j+48]*in2.at(j);
		}
		
		
		for(jj=0;jj<5;jj++)
		{
			for(j=-2;j<=0;j++)
			{
				j2=(j+2)*6;
				pre[jj]+=c_c[jj][288+j2]*hbond_effect.at(id+j).c_length;
				pre[jj]+=c_c[jj][289+j2]*hbond_effect.at(id+j).c_phi;
				pre[jj]+=c_c[jj][290+j2]*hbond_effect.at(id+j).c_psi;
				pre[jj]+=c_c[jj][291+j2]*hbond_effect.at(id+j).n_length;
				pre[jj]+=c_c[jj][292+j2]*hbond_effect.at(id+j).n_phi;
				pre[jj]+=c_c[jj][293+j2]*hbond_effect.at(id+j).n_psi;
			}
		}

			//sequence information
			code=pdb->code(id-1);
			Sequence::code2array(code,buffer);
			for(j=0;j<20;j++)		
			{
				for(jj=0;jj<5;jj++)
					pre[jj]+=c_c[jj][j+306]*buffer[j];
			}

			code=pdb->code(id);
			Sequence::code2array(code,buffer);
			for(j=0;j<20;j++)		
			{
				for(jj=0;jj<5;jj++)
					pre[jj]+=c_c[jj][j+326]*buffer[j];
			}

			code=pdb->code(id+1);
			Sequence::code2array(code,buffer);
			for(j=0;j<20;j++)		
			{
				for(jj=0;jj<5;jj++)
					pre[jj]+=c_c[jj][j+346]*buffer[j];
			}

			if(index.at(i).x2>0)
			{
				temp=0.0;
				for(j=0;j<5;j++)
					temp+=ring_effect.at(index.at(i).x2-1).x[j];
				temp*=c_h_add[4];
				for(j=0;j<4;j++)
					temp+=ani_effect.at(index.at(i).x2-1).x[j]*c_h_add[j];
				pre[3]+=temp;
			}
			pre[5]=999.0;
			pdb->attach_bbprediction(id,pre);
	}

	cal_error();
}






// Updated function for OpenACC
void CMainbody::predict_bb_static_ann()
{
	double st;
	st = omp_get_wtime();

	int i,j;
	int id;
	char code,code_pre,code_fol;
	vector<double> out;
	vector<double> in,in2;
	vector<struct ehbond> hbond_effect(176);
	vector<struct index_two> index;
	vector<int> c1,c2;
	vector<float> result;
	double pre[6];
	vector<double> eca,ecb,eco,eh,en;


	class CAnn ann_ca,ann_cb,ann_co,ann_n,ann_h,ann_ha;


	//ann_ca.load("ann_ca.dat");
	ann_ca.loadp(p_ann_ca);
	ann_cb.loadp(p_ann_cb);
	ann_co.loadp(p_ann_co);
	ann_n.loadp(p_ann_n);
	ann_h.loadp(p_ann_h);
	ann_ha.loadp(p_ann_ha);


	traj->gethbond_acc(hbond_arr, hbond_size, &hbond_effect);
	ehbond *hbond_effect_arr = hbond_effect.data();
	int hbond_effect_size = hbond_effect.size();

	nh_group *bbnh_new = bbnh.data();
	int bbnh_size = bbnh.size();
	#pragma acc enter data copyin(bbnh_new[0:bbnh_size])

	vector<struct double_four> ani_effect(bbnh_size);
	traj->getani_acc(anistropy_new,anistropy_size,bbnh_new,bbnh_size,&ani_effect);
	double_four *ani_effect_arr = ani_effect.data();
	int ani_effect_size = ani_effect.size();

	vector<struct double_five> ring_effect(bbnh_size);
	traj->getring_acc(ring_index_new, ring_index_size, bbnh_new, bbnh_size, &ring_effect);
	double_five *ring_effect_arr = ring_effect.data();
	int ring_effect_size = ring_effect.size();

	#pragma acc exit data delete(bbnh_new) 

	
	//gather all ha protons to calculate ring and ani.
	vector<struct proton> ha_protons;
	struct proton ha;
	for(i=0;i<(int)bb.size();i++)
	{
		
		if(bb.at(i).code=='G')
		{
			ha.nh=2;
			ha.hpos[0]=bb.at(i).hapos;
			ha.hpos[1]=bb.at(i).hapos2;
			ha.exp=(bb.at(i).exp_ha+bb.at(i).exp_ha2)/2.0;
			ha.exp1=bb.at(i).exp_ha;
			ha.exp2=bb.at(i).exp_ha2;
			ha.type=90;
			ha.name="HB2";
			ha.name2="HB3";
			ha.cname="CA";
			ha.cname2="CA";
		}
		else
		{
			ha.nh=1;
			ha.hpos[0]=bb.at(i).hapos;
			ha.exp=bb.at(i).exp_ha;
			ha.type=91;
			ha.name="HA";
			ha.cname="CA";
		}
		ha.cpos=bb.at(i).capos;
		ha.exp_c=bb.at(i).exp_ca;
		ha.id=bb.at(i).id;
		ha.code=bb.at(i).code;

		ha_protons.push_back(ha);
	}

	proton * ha_protons_new = ha_protons.data();
	int ha_protons_size = ha_protons.size();
#pragma acc enter data copyin(ha_protons_new[0:ha_protons_size])

	vector<double_four> ani_effect_ha(ha_protons_size);
	traj->getani_acc(anistropy_new,anistropy_size,ha_protons_new,ha_protons_size,&ani_effect_ha);
	double_four *ani_effect_ha_arr = ani_effect_ha.data();
	int ani_effect_ha_size = ani_effect_ha.size();

	vector<struct double_five> ring_effect_ha(ha_protons_size);
	traj->getring_acc(ring_index_new, ring_index_size, ha_protons_new, ha_protons_size, &ring_effect_ha);
	double_five *ring_effect_ha_arr = ring_effect_ha.data();
	int ring_effect_ha_size = ring_effect_ha.size();

#pragma acc exit data delete(ha_protons_new)


	index.resize(pdb->getnres());
	for(i=0;i<(int)index.size();i++)
		index.at(i).x1=index.at(i).x2=-1;
	for(i=0;i<(int)bb.size();i++)
		index.at(bb.at(i).id-1).x1=i+1;
	for(i=0;i<(int)bbnh.size();i++)
		index.at(bbnh.at(i).id-1).x2=i+1;

	
	c2=pdb->getselect(":1-%@allheavy");	
	int* c2_arr = c2.data();
	int c2_size = c2.size();
	int results_size = (index.size()-2)*3;
	float *results = new float[results_size];
	#pragma acc enter data create(results[0:results_size])
	traj->get_all_contacts(&bb, &index, index.size(),c2_arr,c2_size,results,results_size);


	index_two *index_arr = index.data();
	int index_size = index.size();

	char *v_oln = pdb->getvoneletter();
	int *v_pos = pdb->code_pos;
	int pos;
	int bi;

	int ndihe = dihe_process.ndihe;
	int nframe = dihe_process.nframe;
	double *dihe = dihe_process.dihe->data();
	int dihe_size = dihe_process.dihe->size();
	int *num_arr = dihe_process.num.data();
	int num_size = dihe_process.num.size();
	double *predictions = new double[(index_size-2)*6];
	CAminoacid **v = pdb->v_arr;
	int v_size = pdb->v_size;
	double pre_ca, pre_cb, pre_co, pre_n, pre_h, pre_ha;

	
#pragma acc parallel default(present) copyin(blosum[0:400],v_oln[0:v_size],index_arr[0:index_size],dihe[0:dihe_size],\
num_arr[0:num_size],v_pos[0:v_size]) copyout(predictions[0:(index_size-2)*6])
{
	#pragma acc loop independent gang private(code,code_pre,code_fol,pos,id,pre_ca,pre_cb,pre_co,pre_n,pre_h,pre_ha)
	for(i=0+1;i<index_size-1;i++)
	{
		double oneline[101];
		double oneline_cb[101];
		double oneline_co[101];
		double oneline_n[101];
		double oneline_h[110];
		double oneline_ha[110];
		double out_arr[32];


		if(index_arr[i].x1<0){
			continue;
		}
		id=i+1;
		if((id-1)<0 || (id-1) >v_size-1)
			code='X';
		else
			code=v_oln[id-1];


		if((id-2)<0 || (id-2) >v_size-1)
			code_pre='X';
		else
			code_pre=v_oln[id-2];


		if((id)<0 || (id) >v_size-1)
			code_fol='X';
		else
			code_fol=v_oln[id];


		if(code_pre=='X')
			pos=7;
		else
			pos=v_pos[id-2];


		#pragma acc loop vector independent
		for(j=0;j<20;j++)
		{
			oneline[j] = blosum[pos*20+j];
			oneline_cb[j] = blosum[pos*20+j];
			oneline_co[j] = blosum[pos*20+j];
			oneline_n[j] = blosum[pos*20+j];
			oneline_h[j] = blosum[pos*20+j];
			oneline_ha[j] = blosum[pos*20+j];
		}


		if(code=='X')
			pos=7;
		else
			pos=v_pos[id-1];

		#pragma acc loop vector independent
		for(j=0;j<20;j++)
		{
			oneline[j+20]=blosum[pos*20+j];
			oneline_cb[j+20] = blosum[pos*20+j];
			oneline_co[j+20] = blosum[pos*20+j];
			oneline_n[j+20] =blosum[pos*20+j];
			oneline_h[j+20] = blosum[pos*20+j];
			oneline_ha[j+20] = blosum[pos*20+j];
		}


		if(code_fol=='X')
			pos=7;
		else
			pos=v_pos[id];

		#pragma acc loop vector independent
		for(j=0;j<20;j++)
		{
			oneline[j+40]=blosum[pos*20+j];
			oneline_cb[j+40] = blosum[pos*20+j];
			oneline_co[j+40] = blosum[pos*20+j];
			oneline_n[j+40] = blosum[pos*20+j];
			oneline_h[j+40] = blosum[pos*20+j];
			oneline_ha[j+40] = blosum[pos*20+j];
		}

		ca_ann(id, out_arr, ndihe, nframe, dihe, num_arr, num_size);

		#pragma acc loop vector independent
		for(j=2;j<26;j++)
		{
			oneline[j+60-2+1]=out_arr[j];
			oneline_cb[j+60-2+1] = out_arr[j];
			oneline_co[j+60-2+1] = out_arr[j];
			oneline_n[j+60-2+1] = out_arr[j];
			oneline_h[j+60-2+1] = out_arr[j];
			oneline_ha[j+60-2+1] = out_arr[j];
		}

		#pragma acc loop seq
		for(j=28;j<32;j++)
		{
			oneline[j+84-28+1]=out_arr[j];
			oneline_cb[j+84-28+1] = out_arr[j];
			oneline_co[j+84-28+1] = out_arr[j];
			oneline_n[j+84-28+1] = out_arr[j];
			oneline_h[j+84-28+1] = out_arr[j];
			oneline_ha[j+84-28+1] = out_arr[j];
		}

		//hbond effect, 12 terms
		oneline[88+1]=hbond_effect_arr[id-2].c_length;
		oneline_cb[88+1] = hbond_effect_arr[id-2].c_length;
		oneline_co[88+1] =hbond_effect_arr[id-2].c_length;
		oneline_n[88+1] = hbond_effect_arr[id-2].c_length;
		oneline_h[88+1] = hbond_effect_arr[id-2].c_length;
		oneline_ha[88+1] = hbond_effect_arr[id-2].c_length;

		oneline[89+1]=hbond_effect_arr[id-2].c_phi;
		oneline_cb[89+1] = hbond_effect_arr[id-2].c_phi;
		oneline_co[89+1] =hbond_effect_arr[id-2].c_phi;
		oneline_n[89+1] = hbond_effect_arr[id-2].c_phi;
		oneline_h[89+1] = hbond_effect_arr[id-2].c_phi;
		oneline_ha[89+1] = hbond_effect_arr[id-2].c_phi;

		oneline[90+1]=hbond_effect_arr[id-2].c_psi;
		oneline_cb[90+1] = hbond_effect_arr[id-2].c_psi;
		oneline_co[90+1] =hbond_effect_arr[id-2].c_psi;
		oneline_n[90+1] = hbond_effect_arr[id-2].c_psi;
		oneline_h[90+1] = hbond_effect_arr[id-2].c_psi;
		oneline_ha[90+1] =hbond_effect_arr[id-2].c_psi;

		

		oneline[91+1]=hbond_effect_arr[id-1].c_length;
		oneline[92+1]=hbond_effect_arr[id-1].n_length;
		oneline[93+1]=hbond_effect_arr[id-1].c_phi;
		oneline[94+1]=hbond_effect_arr[id-1].c_psi;
		oneline[95+1]=hbond_effect_arr[id-1].n_phi;
		oneline[96+1]=hbond_effect_arr[id-1].n_psi;

		oneline_cb[91+1]=hbond_effect_arr[id-1].c_length;
		oneline_cb[92+1]=hbond_effect_arr[id-1].n_length;
		oneline_cb[93+1]=hbond_effect_arr[id-1].c_phi;
		oneline_cb[94+1]=hbond_effect_arr[id-1].c_psi;
		oneline_cb[95+1]=hbond_effect_arr[id-1].n_phi;
		oneline_cb[96+1]=hbond_effect_arr[id-1].n_psi;

		oneline_co[91+1]=hbond_effect_arr[id-1].c_length;
		oneline_co[92+1]=hbond_effect_arr[id-1].n_length;
		oneline_co[93+1]=hbond_effect_arr[id-1].c_phi;
		oneline_co[94+1]=hbond_effect_arr[id-1].c_psi;
		oneline_co[95+1]=hbond_effect_arr[id-1].n_phi;
		oneline_co[96+1]=hbond_effect_arr[id-1].n_psi;

		oneline_n[91+1]=hbond_effect_arr[id-1].c_length;
		oneline_n[92+1]=hbond_effect_arr[id-1].n_length;
		oneline_n[93+1]=hbond_effect_arr[id-1].c_phi;
		oneline_n[94+1]=hbond_effect_arr[id-1].c_psi;
		oneline_n[95+1]=hbond_effect_arr[id-1].n_phi;
		oneline_n[96+1]=hbond_effect_arr[id-1].n_psi;

		oneline_h[91+1]=hbond_effect_arr[id-1].c_length;
		oneline_h[92+1]=hbond_effect_arr[id-1].n_length;
		oneline_h[93+1]=hbond_effect_arr[id-1].c_phi;
		oneline_h[94+1]=hbond_effect_arr[id-1].c_psi;
		oneline_h[95+1]=hbond_effect_arr[id-1].n_phi;
		oneline_h[96+1]=hbond_effect_arr[id-1].n_psi;

		oneline_ha[91+1]=hbond_effect_arr[id-1].c_length;
		oneline_ha[92+1]=hbond_effect_arr[id-1].n_length;
		oneline_ha[93+1]=hbond_effect_arr[id-1].c_phi;
		oneline_ha[94+1]=hbond_effect_arr[id-1].c_psi;
		oneline_ha[95+1]=hbond_effect_arr[id-1].n_phi;
		oneline_ha[96+1]=hbond_effect_arr[id-1].n_psi;


		oneline[97+1]=hbond_effect_arr[id].n_length;
		oneline[98+1]=hbond_effect_arr[id].n_phi;
		oneline[99+1]=hbond_effect_arr[id].n_psi;

		oneline_cb[97+1]=hbond_effect_arr[id].n_length;
		oneline_cb[98+1]=hbond_effect_arr[id].n_phi;
		oneline_cb[99+1]=hbond_effect_arr[id].n_psi;

		oneline_co[97+1]=hbond_effect_arr[id].n_length;
		oneline_co[98+1]=hbond_effect_arr[id].n_phi;
		oneline_co[99+1]=hbond_effect_arr[id].n_psi;

		oneline_n[97+1]=hbond_effect_arr[id].n_length;
		oneline_n[98+1]=hbond_effect_arr[id].n_phi;
		oneline_n[99+1]=hbond_effect_arr[id].n_psi;

		oneline_h[97+1]=hbond_effect_arr[id].n_length;
		oneline_h[98+1]=hbond_effect_arr[id].n_phi;
		oneline_h[99+1]=hbond_effect_arr[id].n_psi;

		oneline_ha[97+1]=hbond_effect_arr[id].n_length;
		oneline_ha[98+1]=hbond_effect_arr[id].n_phi;
		oneline_ha[99+1]=hbond_effect_arr[id].n_psi;


		oneline[60]=results[((i-1)*3)+0];
		oneline_n[60]=results[((i-1)*3)+0];
		oneline_h[60]=results[((i-1)*3)+0];
		oneline_ha[60]=results[((i-1)*3)+0];
		oneline_cb[60]=results[((i-1)*3)+1];
		oneline_co[60]=results[((i-1)*3)+2];

		//hn
		if(index_arr[i].x2>0)
		{
			oneline_h[101]=ring_effect_arr[index_arr[i].x2-1].x[0];
			oneline_h[102]=ring_effect_arr[index_arr[i].x2-1].x[1];
			oneline_h[103]=ring_effect_arr[index_arr[i].x2-1].x[2];
			oneline_h[104]=ring_effect_arr[index_arr[i].x2-1].x[3];
			oneline_h[105]=ring_effect_arr[index_arr[i].x2-1].x[4];

			oneline_h[106]=ani_effect_arr[index_arr[i].x2-1].x[0];
			oneline_h[107]=ani_effect_arr[index_arr[i].x2-1].x[1];
			oneline_h[108]=ani_effect_arr[index_arr[i].x2-1].x[2];
			oneline_h[109]=ani_effect_arr[index_arr[i].x2-1].x[3];

			pre_h=ann_h.predict_one_acc(oneline_h,110);
		} else {
			pre_h=-999.0;
		}

		//ha
		oneline_ha[101]=ring_effect_ha_arr[index_arr[i].x1-1].x[0];
		oneline_ha[102]=ring_effect_ha_arr[index_arr[i].x1-1].x[1];
		oneline_ha[103]=ring_effect_ha_arr[index_arr[i].x1-1].x[2];
		oneline_ha[104]=ring_effect_ha_arr[index_arr[i].x1-1].x[3];
		oneline_ha[105]=ring_effect_ha_arr[index_arr[i].x1-1].x[4];

		oneline_ha[106]=ani_effect_ha_arr[index_arr[i].x1-1].x[0];
		oneline_ha[107]=ani_effect_ha_arr[index_arr[i].x1-1].x[1];
		oneline_ha[108]=ani_effect_ha_arr[index_arr[i].x1-1].x[2];
		oneline_ha[109]=ani_effect_ha_arr[index_arr[i].x1-1].x[3];

		pre_ha=ann_ha.predict_one_acc(oneline_ha,110);
		pre_ca=ann_ca.predict_one_acc(oneline,101);
		pre_cb=ann_cb.predict_one_acc(oneline_cb,101);
		pre_co=ann_co.predict_one_acc(oneline_co,101);
		pre_n=ann_n.predict_one_acc(oneline_n,101);

		predictions[((i-1)*6)+0]=pre_ca;
		predictions[((i-1)*6)+1]=pre_cb;
		predictions[((i-1)*6)+2]=pre_co;
		predictions[((i-1)*6)+4]=pre_n;
		predictions[((i-1)*6)+3]=pre_h;
		predictions[((i-1)*6)+5]=pre_ha;
	}
} // end parallel

	for(i=1; i<index_size-1; i++){
		if(index_arr[i].x1<0){
			continue;
		}
		id = i+1;
		pdb->attach_bbprediction(id,predictions[((i-1)*6)+0],predictions[((i-1)*6)+1],predictions[((i-1)*6)+2],
			predictions[((i-1)*6)+4],predictions[((i-1)*6)+3],predictions[((i-1)*6)+5]);
	}

#pragma acc exit data delete(hbond_effect_arr,results,ani_effect_arr,ani_effect_ha_arr,ring_effect_arr,ring_effect_ha_arr)

	cal_error();
};



void CMainbody::predict_bb_static_new()
{
	int i,j,jj;
	int id;
	int jump;
	char code,code_pre,code_fol;
	vector<double> out;
	vector<double> in,in2;
	vector<struct double_five> ring_effect,ring_effect_ha;
	vector<struct ehbond> hbond_effect;
	vector<struct double_four> ani_effect,ani_effect_ha;
	vector<struct index_two> index;
	vector<int> c1,c2;
	vector<float> result;
	double pre[6];
	double temp,temp1,temp2,tt1,tt2;
	vector<double> eca,ecb,eco,eh,en;
	double rc;
	double vhill;
	


	traj->gethbond(&hbond,&hbond_effect);
	traj->getani(&anistropy,&bbnh,&ani_effect);
	traj->getring(&ring_index,&bbnh,&ring_effect);



	

	vector<struct proton> ha_protons;
	struct proton ha;
	for(i=0;i<(int)bb.size();i++)
	{
		
		if(bb.at(i).code=='G')
		{
			ha.nh=2;
			ha.hpos[0]=bb.at(i).hapos;
			ha.hpos[1]=bb.at(i).hapos2;
			ha.exp=(bb.at(i).exp_ha+bb.at(i).exp_ha2)/2.0;
			ha.exp1=bb.at(i).exp_ha;
			ha.exp2=bb.at(i).exp_ha2;
			ha.type=90;
			ha.name="HB2";
			ha.name2="HB3";
			ha.cname="CA";
			ha.cname2="CA";
		}
		else
		{
			ha.nh=1;
			ha.hpos[0]=bb.at(i).hapos;
			ha.exp=bb.at(i).exp_ha;
			ha.type=91;
			ha.name="HA";
			ha.cname="CA";
		}
		ha.cpos=bb.at(i).capos;
		ha.exp_c=bb.at(i).exp_ca;
		ha.id=bb.at(i).id;
		ha.code=bb.at(i).code;

		ha_protons.push_back(ha);
	}



	index.resize(pdb->getnres());	
	for(i=0;i<(int)index.size();i++)
		index.at(i).x1=index.at(i).x2=-1;
	for(i=0;i<(int)bb.size();i++)
		index.at(bb.at(i).id-1).x1=i+1;
	for(i=0;i<(int)bbnh.size();i++)
		index.at(bbnh.at(i).id-1).x2=i+1;
	

	for(i=0+1;i<(int)index.size()-1;i++)
	{
		//cout<<i<<endl;
		if(index.at(i).x1<0)
			continue;
		id=i+1;
		code=pdb->code(id);
		code_pre=pdb->code(id-1);
		code_fol=pdb->code(id+1);

		for(jj=0;jj<6;jj++)
			pre[jj]=0.0;
		jump=0;

		//ca=cb=co=h=n=0.0;

		//missing or unknow residue shouldn't be predicted. This is also true if either previous or following residue is missing (or unknown)
		if(code_pre=='X' || code_pre=='B' || code_fol=='X' || code_fol=='B'|| code=='X' || code=='B')
		{
			for(jj=0;jj<6;jj++)
				pre[jj]=-999.0;
			continue;
		}

		//sequence information
		Sequence::code2array(code_pre,buffer);
		for(j=0;j<20;j++)		
		{
			for(jj=0;jj<6;jj++)
				pre[jj]+=static_c[jj][j+jump]*buffer[j];
		}
		jump+=20;
		
		code=pdb->code(id);
		Sequence::code2array(code,buffer);
		for(j=0;j<20;j++)		
		{
			for(jj=0;jj<6;jj++)
				pre[jj]+=static_c[jj][j+jump]*buffer[j];
		}
		jump+=20;

		Sequence::code2array(code_fol,buffer);
		for(j=0;j<20;j++)		
		{
			for(jj=0;jj<6;jj++)
				pre[jj]+=static_c[jj][j+jump]*buffer[j];
		}
		jump+=20;



		//dihedral angle contribution!
		if(dihe_process.ca_static_new(id)==0) //==0 means missing dihedral angles in this calculation !!
		{
			for(jj=0;jj<6;jj++)
				pre[jj]=-999.0;
			continue;
		}

		out=dihe_process.output();

		//pre dihe  12*20 terms
		in.clear();
		for(j=0;j<12;j++)
		{
			in.push_back(out.at(j));
		}
		in2=Sequence::expand(code_pre,&in);
		for(j=0;j<(int)in2.size();j++)
		{
			for(jj=0;jj<6;jj++)
				pre[jj]+=static_c[jj][j+jump]*in2.at(j);
		}
		jump+=240;


		//self dihe 24*20 terms
		in.clear();
		for(j=0;j<24;j++)
		{
			in.push_back(out.at(j+12));
		}
		in2=Sequence::expand(code,&in);
		for(j=0;j<(int)in2.size();j++)
		{
			for(jj=0;jj<6;jj++)
				pre[jj]+=static_c[jj][j+jump]*in2.at(j);
		}
		jump+=480;
		
		//follow dihe  12*20 terms
		in.clear();
		for(j=0;j<12;j++)
		{
			in.push_back(out.at(j+36));
		}
		in2=Sequence::expand(code_fol,&in);
		for(j=0;j<(int)in2.size();j++)
		{
			for(jj=0;jj<6;jj++)
				pre[jj]+=static_c[jj][j+jump]*in2.at(j);
		}
		jump+=240;
		
	
		//hbond effect, length only
		for(jj=0;jj<6;jj++)
		{	
			tt1=temp1=hbond_effect.at(id-1).c_length;
			tt2=temp2=hbond_effect.at(id-1).n_length;
			pre[jj]+=static_c[jj][jump+0]*temp1;
			pre[jj]+=static_c[jj][jump+1]*temp2;

			temp1*=temp1;temp2*=temp2;
			pre[jj]+=static_c[jj][jump+2]*temp1;
			pre[jj]+=static_c[jj][jump+3]*temp2;

			pre[jj]+=static_c[jj][jump+4]*temp1*tt1;
			pre[jj]+=static_c[jj][jump+5]*temp2*tt2;

			temp1*=temp1;temp2*=temp2;
			pre[jj]+=static_c[jj][jump+6]*temp1;
			pre[jj]+=static_c[jj][jump+7]*temp2;
		}
		jump+=8;

						
		if(index.at(i).x2>0)
		{
			temp=0.0;
			for(j=0;j<5;j++)
				temp+=ring_effect.at(index.at(i).x2-1).x[j]*static_h[0][j];
			for(j=0;j<4;j++)
				temp+=ani_effect.at(index.at(i).x2-1).x[j]*static_h[0][j+5];
			pre[3]+=temp;
		}

		temp=0.0;
		for(j=0;j<5;j++)
			temp+=ring_effect_ha.at(index.at(i).x1-1).x[j]*static_h[1][j];
		for(j=0;j<4;j++)
			temp+=ani_effect_ha.at(index.at(i).x1-1).x[j]*static_h[1][j+5];
		pre[5]+=temp;


		//get contact sum adjustment
		c1.clear();
		c1.push_back(bb.at(index.at(i).x1-1).capos);
		c1.push_back(bb.at(index.at(i).x1-1).cbpos);
		c1.push_back(bb.at(index.at(i).x1-1).copos);
		c2=pdb->getselect(":1-%@allheavy");
		int* c2_arr = c2.data();
		int c2_size = c2.size();

		result.clear();
		traj->get_contact(c1,c2_arr, c2_size,&result);
		result.push_back(result.at(0));
		result.push_back(result.at(0));
		result.push_back(result.at(0));

		for(jj=0;jj<6;jj++)
		{
			vhill=hill(result.at(jj),hill_para[jj][0],hill_para[jj][1]);
			rc=(1-vhill)/vhill;
			pre[jj]+=rc*static_c[jj][Sequence::code2pos(code)+jump];
			pre[jj]*=vhill;
		}
		
		pdb->attach_bbprediction(id,pre);
	}
	
	cal_error();

}


void CMainbody::cal_error()
{
	int j;
	vector< vector<double> > cas,cbs,cos,hs,ns,has;
	double e,t,w;


	cas.resize(2);cbs.resize(2);cos.resize(2);hs.resize(2);ns.resize(2);has.resize(2);
	bb.clear();
	pdb->getbb(&bb);
	clear(bb);
	for(j=0;j<(int)bb.size();j++)
	{
		cas.at(0).push_back(bb.at(j).exp_ca);
		cbs.at(0).push_back(bb.at(j).exp_cb);
		cos.at(0).push_back(bb.at(j).exp_co);
		hs.at(0).push_back(bb.at(j).exp_h);
		ns.at(0).push_back(bb.at(j).exp_n);
		has.at(0).push_back(bb.at(j).exp_ha);
	
		cas.at(1).push_back(bb.at(j).pre_ca);
		cbs.at(1).push_back(bb.at(j).pre_cb);
		cos.at(1).push_back(bb.at(j).pre_c);
		hs.at(1).push_back(bb.at(j).pre_h);
		ns.at(1).push_back(bb.at(j).pre_n);
		has.at(1).push_back(bb.at(j).pre_ha);

	}
	compare("CA",cas);compare("CB",cbs);compare("C'",cos);compare("HN",hs);compare("N",ns);compare("HA",has);


	FILE *fp=fopen("cs_rmsd.dat","wt");

	for(j=0;j<(int)bb.size();j++)
	{
		e=0;
		w=0;
		if(fabs(bb.at(j).exp_ca)<490.0 && fabs(bb.at(j).pre_ca)<490.0)
		{
			t=bb.at(j).exp_ca-bb.at(j).pre_ca;
			e+=t*t;
			w+=1.0;
		}
		if(fabs(bb.at(j).exp_cb)<490.0 && fabs(bb.at(j).pre_cb)<490.0)
		{
			t=bb.at(j).exp_cb-bb.at(j).pre_cb;
			e+=t*t;
			w+=1.0;
		}
		if(fabs(bb.at(j).exp_co)<490.0 && fabs(bb.at(j).pre_c)<490.0)
		{
			t=bb.at(j).exp_co-bb.at(j).pre_c;
			e+=t*t;
			w+=1.0;
		}
		if(fabs(bb.at(j).exp_h)<490.0 && fabs(bb.at(j).pre_h)<490.0)
		{
			t=bb.at(j).exp_h-bb.at(j).pre_h;
			e+=t*t*4;
			w+=4.0;
		}
		if(fabs(bb.at(j).exp_n)<490.0 && fabs(bb.at(j).pre_n)<490.0)
		{
			t=bb.at(j).exp_n-bb.at(j).pre_n;
			e+=t*t*0.4;
			w+=0.4;
		}
		if(fabs(bb.at(j).exp_ha)<490.0 && fabs(bb.at(j).pre_ha)<490.0)
		{
			t=bb.at(j).exp_ha-bb.at(j).pre_ha;
			e+=t*t*4;
			w+=4;
		}
		if(w>0.0)
		{
			e=sqrt(e/w);
			fprintf(fp,"%10d%10s%10.2f\n",bb.at(j).id0,Sequence::code2name(bb.at(j).code).c_str(),e);
		}
		else
			e=0.0;

	}		
	fclose(fp);		

}


double CMainbody::hill(double contact, double v1, double v2)
{
	double hill;
	
	if(contact<v1)	
		hill=0.0;
	else if(contact>v2)
		hill=1.0;
	else
	{
		hill=(contact-v1)/(v2-v1)*6-3;
		hill=(1-exp(-2*hill))/(1+exp(-2*hill));
		hill=hill/2+0.5;
	}

	return hill;

}


void CMainbody::predict_bb2()
{
	int i,j,j2,jj,n;
	int id;
	char code;
	vector< vector<double> > out;
	vector<double> in,in2;
	vector< vector<struct double_five> > ring_effect;
	vector< vector<struct ehbond> > hbond_effect;
	vector< vector<struct double_four> > ani_effect;
	vector<struct index_two> index;
	vector<double> cs_ca,cs_cb,cs_c,cs_h,cs_n;
	double pre[6];
	double temp;
	vector<double> eca,ecb,eco,eh,en;
	char name[4];
	FILE *fp;

	fp=fopen("bb_details.dat","w");
	


	traj->gethbond(&hbond,&hbond_effect);
	traj->getani(&anistropy,&bbnh,&ani_effect);
	traj->getring(&ring_index,&bbnh,&ring_effect);

	

	index.resize(pdb->getnres());	
	for(i=0;i<(int)index.size();i++)
		index.at(i).x1=index.at(i).x2=-1;
	for(i=0;i<(int)bb.size();i++)
		index.at(bb.at(i).id-1).x1=i+1;
	for(i=0;i<(int)bbnh.size();i++)
		index.at(bbnh.at(i).id-1).x2=i+1;


	

	for(i=0+1;i<(int)index.size()-1;i++)
	{	//cout<<"i is "<<i<<endl;
		if(index.at(i).x1<0)
			continue;
		id=i+1;
		code=pdb->code(id);
		pdb->name(id,name);

		if(strcmp(name,"CYS")==0)
			continue;

		dihe_process.ca2(id); 
		out=dihe_process.output2();

		cs_ca.clear();
		cs_cb.clear();
		cs_c.clear();
		cs_h.clear();
		cs_n.clear();

		for(n=0;n<nconf;n++)
		{		
			id=i+1;
			code=pdb->code(id);

			for(jj=0;jj<5;jj++)
				pre[jj]=0.0;
	 
			for(j=0;j<30;j++)
			{
				for(jj=0;jj<5;jj++)
					pre[jj]+=c_c[jj][j]*out.at(j).at(n);
			}

			for(j=42;j<60;j++)
			{
				for(jj=0;jj<5;jj++)
					pre[jj]+=c_c[jj][j-12]*out.at(j).at(n);
			}
				

			in.clear();
			for(j=0;j<12;j++)
			{
				in.push_back(out.at(j+30).at(n));
			}
			in2=Sequence::expand(code,&in);


			for(j=0;j<(int)in2.size();j++)
			{
				for(jj=0;jj<5;jj++)
					pre[jj]+=c_c[jj][j+48]*in2.at(j);
			}
			
		
			for(jj=0;jj<5;jj++)
			{
				for(j=-2;j<=0;j++)
				{
					j2=(j+2)*6;
					pre[jj]+=c_c[jj][288+j2]*hbond_effect.at(id+j).at(n).c_length;
					pre[jj]+=c_c[jj][289+j2]*hbond_effect.at(id+j).at(n).c_phi;
					pre[jj]+=c_c[jj][290+j2]*hbond_effect.at(id+j).at(n).c_psi;
					pre[jj]+=c_c[jj][291+j2]*hbond_effect.at(id+j).at(n).n_length;
					pre[jj]+=c_c[jj][292+j2]*hbond_effect.at(id+j).at(n).n_phi;
					pre[jj]+=c_c[jj][293+j2]*hbond_effect.at(id+j).at(n).n_psi;
				}
			}

			//sequence information
			code=pdb->code(id-1);
			Sequence::code2array(code,buffer);
			for(j=0;j<20;j++)		
			{
				for(jj=0;jj<5;jj++)
					pre[jj]+=c_c[jj][j+306]*buffer[j];
			}

			code=pdb->code(id);
			Sequence::code2array(code,buffer);
			for(j=0;j<20;j++)		
			{
				for(jj=0;jj<5;jj++)
					pre[jj]+=c_c[jj][j+326]*buffer[j];
			}

			code=pdb->code(id+1);
			Sequence::code2array(code,buffer);
			for(j=0;j<20;j++)		
			{
				for(jj=0;jj<5;jj++)
					pre[jj]+=c_c[jj][j+346]*buffer[j];
			}


			if(index.at(i).x2>0)
			{
				temp=0.0;
				for(j=0;j<5;j++)
					temp+=ring_effect.at(index.at(i).x2-1).at(n).x[j];
				temp*=c_h_add[4];
				for(j=0;j<4;j++)
					temp+=ani_effect.at(index.at(i).x2-1).at(n).x[j]*c_h_add[j];
				pre[3]+=temp;
			}

			cs_ca.push_back(pre[0]);
			cs_cb.push_back(pre[1]);
			cs_c.push_back(pre[2]);
			cs_n.push_back(pre[4]);
			cs_h.push_back(pre[3]);
		}
	
		fprintf(fp,"%8d%8s",id,name);
		fprintf(fp,"   CA   ");
		fprintf(fp,"%8.3f",bb.at(index.at(i).x1-1).exp_ca);
		for(n=0;n<nconf;n++)
				fprintf(fp,"%8.3f",cs_ca.at(n));

		fprintf(fp,"\n%8d%8s",id,name);
		fprintf(fp,"   CB   ");
		fprintf(fp,"%8.3f",bb.at(index.at(i).x1-1).exp_cb);
		for(n=0;n<nconf;n++)
				fprintf(fp,"%8.3f",cs_cb.at(n));

		fprintf(fp,"\n%8d%8s",id,name);
		fprintf(fp,"   C    ");	
		fprintf(fp,"%8.3f",bb.at(index.at(i).x1-1).exp_co);
		for(n=0;n<nconf;n++)
				fprintf(fp,"%8.3f",cs_c.at(n));

		fprintf(fp,"\n%8d%8s",id,name);
		fprintf(fp,"   H    ");
		if(index.at(i).x2>1)
			fprintf(fp,"%8.3f",bbnh.at(index.at(i).x2-1).exp_h);
		else
			fprintf(fp,"%8.3f",999.0);
		for(n=0;n<nconf;n++)
				fprintf(fp,"%8.3f",cs_h.at(n));

		fprintf(fp,"\n%8d%8s",id,name);
		fprintf(fp,"   N    ");
		if(index.at(i).x2>1)
			fprintf(fp,"%8.3f",bbnh.at(index.at(i).x2-1).exp_n);
		else
			fprintf(fp,"%8.3f",999.0);
		for(n=0;n<nconf;n++)
				fprintf(fp,"%8.3f",cs_n.at(n));
		fprintf(fp,"\n");


		for(n=0;n<5;n++)
			pre[n]=0;
		for(n=0;n<nconf;n++)
		{
			pre[0]+=cs_ca.at(n);
			pre[1]+=cs_cb.at(n);
			pre[2]+=cs_c.at(n);
			pre[3]+=cs_h.at(n);
			pre[4]+=cs_n.at(n);
		}
		for(n=0;n<5;n++)
			pre[n]/=nconf;
		pre[5]=999.0;
		pdb->attach_bbprediction(id,pre);
	}
	return;
}






void CMainbody::predict_proton()
{
	int i,j;
	int id;
	char code;
	int type;

	double c_ring[5]={-0.1939,-0.1629,-0.1620,-0.1872,-0.1987};
	double c_ani[4]={0.004478932,0.000925861,0.001692256,0.00040848288};
	double c_rand[10]={1.3876,1.2671,2.0763,0.0000,1.0064,0.9611,0.9102,0.8890,0.9976,0.8584};

	
	double cs_ring,cs_ani,cs_rand,cs;


	vector< vector<double> > hs;
	hs.resize(2);


	vector<struct double_five> ring_effect;
	vector<struct double_four> ani_effect;

	traj->getani(&anistropy,&protons,&ani_effect);
	traj->getring(&ring_index,&protons,&ring_effect);
		


	for(i=0;i<(int)protons.size();i++)
	{	
		id=protons.at(i).id;
		code=protons.at(i).code;	
		type=protons.at(i).type;


		cs_ring=0;
		for(j=0;j<5;j++)
			cs_ring+=ring_effect.at(i).x[j]*c_ring[j];

		cs_ani=0;
		for(j=0;j<4;j++)
			cs_ani+=ani_effect.at(i).x[j]*c_ani[j];
		cs_rand=c_rand[type-1];
		cs=cs_ring+cs_ani+cs_rand;
		pdb->attach_protonprediction(id,protons.at(i).name,cs);

		hs.at(0).push_back(protons.at(i).exp);
		hs.at(1).push_back(cs);
	}

	compare("Methyl 1H",hs);
	return;
}

void CMainbody::predict_proton2()
{
	int i,j,ii;
	int id;
	char code;
	int type;
	double cs_rand;
	vector<double> cs_ring,cs_ani,cs;
	double ccs_ring,ccs_ani,ccs;
	FILE *fp;

	fp=fopen("proton_details.dat","w");

	double c_ring[5]={-0.1939,-0.1629,-0.1620,-0.1872,-0.1987};
	double c_ani[4]={0.004478932,0.000925861,0.001692256,0.00040848288};
	double c_rand[10]={1.3876,1.2671,2.0763,0.0000,1.0064,0.9611,0.9102,0.8890,0.9976,0.8584};

	vector< vector<struct double_five> > ring_effect;
	vector< vector<struct double_four> > ani_effect;
	traj->getani(&anistropy,&protons,&ani_effect);
	traj->getring(&ring_index,&protons,&ring_effect);

	for(i=0;i<(int)protons.size();i++)
	{	
		id=protons.at(i).id;
		code=protons.at(i).code;	
		type=protons.at(i).type;

		for(ii=0;ii<nconf;ii++)
		{
			cs_ring.push_back(0.0);
			cs_ani.push_back(0.0);
			cs.push_back(0.0);
		}
		for(ii=0;ii<nconf;ii++)
		{
			cs_ring[ii]=0;
			for(j=0;j<5;j++)
				cs_ring[ii]+=ring_effect.at(i).at(ii).x[j]*c_ring[j];
			cs_ani[ii]=0;
			for(j=0;j<4;j++)
				cs_ani[ii]+=ani_effect.at(i).at(ii).x[j]*c_ani[j];
			cs_rand=c_rand[type-1];
			cs[ii]=cs_ring[ii]+cs_ani[ii]+cs_rand;
		}
		fprintf(fp,"%8d %8s %8s",id,Sequence::code2name(code).c_str(),protons.at(i).name.c_str());
		if(fabs(protons.at(i).exp)>0.00001)
			fprintf(fp," %8.3f",protons.at(i).exp);
		else
			fprintf(fp," %8.3f",999.9);
		for(ii=0;ii<nconf;ii++)
			fprintf(fp," %8.3f",cs[ii]);
		fprintf(fp,"\n");


		ccs=ccs_ring=ccs_ani=0.0;
		for(ii=0;ii<nconf;ii++)
		{	
			ccs+=cs[ii];
			ccs_ring+=cs_ring[ii];
			ccs_ani+=cs_ani[ii];
		}
		ccs/=nconf;
		ccs_ring/=nconf;
		ccs_ani/=nconf;
		pdb->attach_protonprediction(id,protons.at(i).name,ccs);
	}
	fclose(fp);
	return;
}


void CMainbody::predict_proton_static_new(void)
{
	double st = omp_get_wtime();
	int i,j;
	int id;
	int type;
	vector<double> out;
	int jump;
	float pre;
	double *c;
	vector< vector<double> > hs;

	vector<struct double_five> ring_effect(allprotons3_size);
	vector<struct double_four> ani_effect(allprotons3_size);

	allprotons=allprotons3;
	double_four* tmp = ani_effect.data();

	traj->getani_acc(anistropy_new,anistropy_size,allprotons3_new,allprotons3_size,&ani_effect);
	traj->getring_acc(ring_index_new, ring_index_size, allprotons3_new, allprotons3_size, &ring_effect);
	hs.resize(2);
	double_five *ring_effect_arr = ring_effect.data();
	int ring_effect_size = ring_effect.size();
	double_four *ani_effect_arr = ani_effect.data();
	int ani_effect_size = ani_effect.size();

	#pragma acc exit data copyout(ani_effect_arr[0:ani_effect_size],ring_effect_arr[0:ring_effect_size])

	for(i=0;i<(int)allprotons.size();i++)
	{		
		type=allprotons.at(i).type;
		id=allprotons.at(i).id;
		
		//too few data point and bad fitting !
		if(type==18 || type==19 || type==28 || type==49 || type==64  || type==82  || type==84  || type==97 )
			continue;
		//bad fitting !
		if(type==36 || type==48 || type==62 )
			continue;

		//use gobal or individual fitting?
		if(sep_table[type]>=0)
			c=c_sep[type-1];
		else
			c=c_all;

		jump=0;
		pre=c[type-1];
		jump+=98;

		for(j=0;j<5;j++)
			pre+=ring_effect.at(i).x[j]*c[j+jump];
		jump+=5;

		for(j=0;j<4;j++)
			//pre+=vec_arr[i*4+j]*c[j+jump];
			pre+=ani_effect.at(i).x[j]*c[j+jump];
		jump+=4;

		dihe_process.allproton(allprotons.at(i).id); 
		dihe_process.hb_expand(type);
		out=dihe_process.output();
		for(j=0;j<8*18;j++)
			pre+=out.at(j)*c[j+jump];
		jump+=8*18;

		pdb->attach_protonprediction(id,allprotons.at(i).name,pre);
		if(allprotons.at(i).name2!="")
			pdb->attach_protonprediction(id,allprotons.at(i).name2,pre);

		hs.at(0).push_back(allprotons.at(i).exp);
		hs.at(1).push_back(pre);
	}

	compare("Side chain protons",hs);
	cout << "predict_proton_static_new: " << omp_get_wtime() - st << " seconds" << endl;
	return;
}



void CMainbody::compare(char * buff, vector< vector<double> > t)
{
	int i,j,n;
	double e,ee;
	bool b;
	vector<bool> bs;

	for(i=0;i<(int)t.at(0).size();i++)
	{
		b=0;
		for(j=0;j<(int)t.size();j++)
		{
			if(t.at(j).at(i)<-400.0 || t.at(j).at(i)>400.0 )
				b=1;
		}
		bs.push_back(b);
		
	}

	cout<<buff<<": ";
	for(i=1;i<(int)t.size();i++)
	{
		ee=0;
		n=0;
		for(j=0;j<(int)t.at(0).size();j++)
		{
			if(bs.at(j)==0)
			{
				e=t.at(0).at(j)-t.at(i).at(j);
				ee+=e*e;
				n++;
			}
		}
		if(n>0)
		{
			ee/=n;
			ee=sqrt(ee);
			cout<<" "<<ee;
		}
		else
			cout<<" N.A.";
	}
	cout<<endl;
}

//This file include all trained parameters in PPM.
#include "data.h"


