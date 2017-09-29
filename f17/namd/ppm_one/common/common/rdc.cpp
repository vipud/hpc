#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>

using namespace std;

#include "rdc.h"



float CNmr::pythag(float p, float q)
{
   float r;
   if (p == 0) return (q);
   if (q == 0) return (p);
   if (p < 0) p = -p;
   if (q < 0) q = -q;
   if (p < q) {r = p; p = q; q = r;}
   r = q / p;
   return p * sqrt (1 + r * r);
}

float * CNmr::myvector(int nl, int nh)
{
	float *v;
	v=new float[nh-nl+1];
	return v-nl;
}

void CNmr::free_vector(float *v, int nl, int nh)
{
	delete [] (v+nl);
	return;
}

float ** CNmr::matrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	float **m;

	m=new float * [nrh-nrl+1];
	m -= nrl;

	for(i=nrl;i<=nrh;i++)
	{
		m[i]=new float[nch-ncl+1];
		m[i] -= ncl;
	}
	return m;
}

void CNmr::free_matrix(float **m,int nrl, int nrh, int ncl, int nch)
{
	int i;
	for(i=nrh;i>=nrl;i--)
	{
		delete [] (m[i]+ncl);
	}
	delete [] (m+nrl);
	return;
}

void CNmr::svdcmp(float **a, int m, int n, float w[], float **v)
{
	int flag,i,its,j,jj,k,l,nm;
	float anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=myvector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) 
				scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				if(f>0)
					g=-sqrt(s);
				else
					g=+sqrt(s);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				if(f>0)
					g=-sqrt(s);
				else
					g=+sqrt(s);	
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		if((fabs(w[i])+fabs(rv1[i]))>anorm)
			anorm=(fabs(w[i])+fabs(rv1[i]));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((float)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) printf("no convergence in 30 svdcmp iterations\n");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			if(f>0) f=((x-z)*(x+z)+h*((y/(f+g))-h))/x;
			else f=((x-z)*(x+z)+h*((y/(f-g))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
}

void CNmr::svbksb(float **u, int m, int n,float w[], float **v, float b[], float x[])
{
	int jj,j,i;
	float s,*tmp;

	tmp=myvector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	free_vector(tmp,1,n);
}

//a[1,n][1,n] real, systematic matrix --> orthogonal matrix
//d[1,n] --> is the diagonal elements.
//e[1,n] --> is the off-diagonal elements.
void CNmr::tred2(float **a, int n, float *d, float *e)
{
	int l,k,j,i;
	float scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[1]=0.0;
	e[1]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
		l=i-1;
		if (d[i]) {
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}


//d[1,n] is the diagonal elements.  -> eigen values
//e[1,n] is the off-diagonal elements
//z[1,n][1,n] is the output of tred2 -> eigen vectors
void CNmr::tqli(float d[], float e[], int n, float **z)
{
          int m,l,iter,i,k;
          float s,r,p,g,f,dd,c,b;
         
          for(i=2;i<=n;i++)
          	e[i-1]=e[i];
          e[n]=0.0;
          for(l=1;l<=n;l++)
          {
          	iter=0;
            do
            {
            	for (m=l;m<=n-1;m++)
            	{
                	dd=fabs(d[m])+fabs(d[m+1]);
                	if ((float)(fabs(e[m])+dd) == dd) break;
              	}
              	if (m != l)
              	{
	                if (iter++ == 30) 
	                  {printf("Too many iterations in tqli");exit;}
	                g=(d[l+1]-d[l])/(2.0*e[l]);
	                r=pythag(g,1.0);
	                if(g>=0) g=d[m]-d[l]+e[l]/(g+fabs(r));
	                else g=d[m]-d[l]+e[l]/(g-fabs(r));
	                s=c=1.0;
	                p=0.0;
                	for (i=m-1;i>=l;i--)
                	{
	                  f=s*e[i];
	                  b=c*e[i];
	                  e[i+1]=(r=pythag(f,g));
	                  if (r == 0.0)
	                  {
	                    d[i+1] -= p;
	                    e[m]=0.0;
	                    break;
	                  }
	                  s=f/r;
	                  c=g/r;
	                  g=d[i+1]-p;
	                  r=(d[i]-g)*s+2.0*c*b;
	                  d[i+1]=g+(p=s*r);
	                  g=c*r-b;
					  for (k=1;k<=n;k++)
                      {
                        f=z[k][i+1];
                        z[k][i+1]=s*z[k][i]+c*f;
                        z[k][i]=c*z[k][i]-s*f;
                      }
                    
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l]=g;
                e[m]=0.0;
              }
            } while (m != l);
          }
        }


int CNmr::loadpdb(CPdb *p_pdb, CTraj * p_traj)
{
	pdb=p_pdb;
	traj=p_traj;

	natom=pdb->getnatom();
	nres=pdb->getnres();
	nconf=traj->getnframe();

	return natom;
}



int CNmr::loadpdb(string name)
{
	pdb=new CPdb;
	traj=new CTraj;
	bnew=1;

	natom=pdb->loadpdb(name);
	nres=pdb->getnres();
	traj->setnres(nres);
	traj->setnatom(natom);
	nconf=traj->loadcoor(name);
	return natom;
}

CNmr::CNmr() 
{
	bnew=0;
};
CNmr::~CNmr() 
{
	if(bnew)
	{
		delete pdb;
		delete traj;
	}
	
};

//////////////////////////////


float CRdc::backcal(float *q,float **u,int m,int n,float b[], float b2[], float x[])
{
	int i,j;
	float err,e;
	float bb;

	
	err=0.0;
	for(i=1;i<=m;i++)
	{
		bb=0;
		for(j=1;j<=n;j++)
		{
			bb+=u[i][j]*x[j];
		}
		b2[i]=bb;
		e=bb-b[i];
		err+=e*e;
	}

    bb=0;
    for(i=1;i<=m;i++)
    {
		bb+=b[i]*b[i];
    }
    *q=sqrt(err/bb);
	err/=m;
	err=sqrt(err);
	return err;
}

float CRdc::error()
{
    return error(1);
}


float CRdc::error(int flag)
{
	int ii,i,j,m,nn;
	int shift;
	double x1,x2,y1,y2,z1,z2;
	float r,q;
	float t1,t2,t3,t4,t5;
	float **u,**uu,*w,**v; 
	float *b,*b2,*x;
	float *xy,*yz,*zx,*x_z,*y_z;
	vector<int> h,n;
	float qq;
    ofstream ftensor;
    ifstream ftensorin;
    int nexp;
    float *tensor;
    
    if(flag==1)
    {
        ftensor.open("rdc_tensor.dat");
        ftensor<<exp.size()<<endl;
        nexp=exp.size();
    }
    else
    {
        ftensorin.open("rdc_tensor.dat");
        ftensorin>>nexp;
        tensor=new float(nexp*5);
        for(i=0;i<nexp*5;i++)
            ftensorin>>tensor[i];
        ftensorin.close();
    }

	qq=0;
	for(ii=0;ii<exp.size() && ii<nexp;ii++)
	{
		nn=exp.at(ii).size();

		xy=myvector(0,nconf*nn-1);
		yz=myvector(0,nconf*nn-1);
		zx=myvector(0,nconf*nn-1);
		x_z=myvector(0,nconf*nn-1);
		y_z=myvector(0,nconf*nn-1);


		h.clear();
		n.clear();
		for(j=0;j<nn;j++)
		{
			h.push_back(exp.at(ii).at(j).atomid1);
			n.push_back(exp.at(ii).at(j).atomid2);
		}

		for(i=0;i<nconf;i++)
		{	
			m=i*nn;		
			for(j=0;j<nn;j++)
			{
				traj->getcoor(h.at(j),i,&x1,&y1,&z1);
				traj->getcoor(n.at(j),i,&x2,&y2,&z2);

				x1-=x2;	y1-=y2;z1-=z2;
				r=sqrt(x1*x1+y1*y1+z1*z1);
				x1/=r;y1/=r;z1/=r;
					
				t1=x1*x1-z1*z1;
				t2=y1*y1-z1*z1;
				t3=2*x1*y1;
				t4=2*x1*z1;
				t5=2*y1*z1;				
				x_z[m+j]=t1;
				y_z[m+j]=t2;
				xy[m+j]=t3;
				zx[m+j]=t4;
				yz[m+j]=t5;				
			}
		}

		u=matrix(1,nn,1,5);
		uu=matrix(1,nn,1,5);
		v=matrix(1,5,1,5);
		w=myvector(1,nn);	
		b=myvector(1,nn);
		b2=myvector(1,nn);
		x=myvector(1,5);

		for(i=1;i<=nn;i++)
		{
			for(j=1;j<=5;j++)
				u[i][j]=0.0;
		}		
	
		for(i=0;i<nconf;i++)
		{
			shift=i*nn;
			for(j=0;j<nn;j++)
			{
				u[j+1][1]+=x_z[shift+j];
				u[j+1][2]+=y_z[shift+j];
				u[j+1][3]+=xy[shift+j];
				u[j+1][4]+=yz[shift+j];
				u[j+1][5]+=zx[shift+j];
			}						
		}

		for(j=0;j<nn;j++)
		{
			u[j+1][1]/=nconf;
			u[j+1][2]/=nconf;
			u[j+1][3]/=nconf;
			u[j+1][4]/=nconf;
			u[j+1][5]/=nconf;
			uu[j+1][1]=u[j+1][1];
			uu[j+1][2]=u[j+1][2];
			uu[j+1][3]=u[j+1][3];
			uu[j+1][4]=u[j+1][4];
			uu[j+1][5]=u[j+1][5];
		}
		
		for(i=0;i<nn;i++)
			b[i+1]=exp.at(ii).at(i).vrdc;

		q=0;
		r=0;
		
            if(flag==1 && nn>10)
            {
                svdcmp(u,nn,5,w,v);
                svbksb(u,nn,5,w,v,b,x);
                ftensor<<x[1]<<" "<<x[2]<<" "<<x[3]<<" "<<x[4]<<" "<<x[5]<<endl;
                r=backcal(&q,uu,nn,5,b,b2,x);
            }
            else
            {
                x[1]=tensor[ii*5];x[2]=tensor[ii*5+1];x[3]=tensor[ii*5+2];x[4]=tensor[ii*5+3];x[5]=tensor[ii*5+4];
                r=backcal(&q,uu,nn,5,b,b2,x);
            }
			
			//printf("Set %d, r=%f q=%f\n",ii,r,q);
		

		fout<<"Set "<<ii<<" r="<<r<<" q="<<q<<endl;
		qq+=q;
		for(i=0;i<nn;i++)
		{
			exp.at(ii).at(i).pre_rdc=b2[i+1];
			fout<<exp.at(ii).at(i).resid1<<" "<<exp.at(ii).at(i).name1<<" ";
			fout<<exp.at(ii).at(i).resid2<<" "<<exp.at(ii).at(i).name2<<" ";
			fout<<exp.at(ii).at(i).vrdc<<" "<<exp.at(ii).at(i).pre_rdc<<endl;
		}
		fout<<endl;


		free_vector(w,1,nn);
		free_vector(b,1,nn);
		free_vector(b2,1,nn);
		free_vector(x,1,5);
		free_matrix(u,1,nn,1,5);
		free_matrix(uu,1,nn,1,5);
		free_matrix(v,1,5,1,5);
		free_vector(x_z,0,nn*nconf-1);
		free_vector(y_z,0,nn*nconf-1);
		free_vector(xy,0,nn*nconf-1);
		free_vector(yz,0,nn*nconf-1);
		free_vector(zx,0,nn*nconf-1);
	}
	
    if(flag==1) ftensor.close();
    
    return qq/ii;
}

void CRdc::fillin(void)
{
	int i,j;
	struct noeatoms noe1,noe2;
	struct exp_rdc rdc_line;

	for(i=0;i<exp.size();i++)
	{
		for(j=exp.at(i).size()-1;j>=0;j--)
		{
			rdc_line=exp.at(i).at(j);

			if(abs(rdc_line.resid1-rdc_line.resid2)<2 && rdc_line.resid1-1>=0 && rdc_line.resid2-1>=0 && rdc_line.resid1-1<nres && rdc_line.resid2-1<nres) 
			{
				noe1=pdb->query(rdc_line.resid1,rdc_line.name1);
				noe2=pdb->query(rdc_line.resid2,rdc_line.name2);
				if(noe1.atoms.size()>0 && noe2.atoms.size()>0)
				{
					exp.at(i).at(j).atomid1=noe1.atoms.at(0).at(0);
					exp.at(i).at(j).atomid2=noe2.atoms.at(0).at(0);
				}
				else
				{
					fout<<"In CRdc, exp set "<<i<<" res "<<rdc_line.resid1<<" name "<<rdc_line.name1<<" and ";
					fout<<" res "<<rdc_line.resid2<<" name "<<rdc_line.name2<<" cannot match pdb"<<endl;
					exp.at(i).erase(exp.at(i).begin()+j);
				}
			}
			else
			{
				fout<<"In CRdc, exp set "<<i<<" res "<<rdc_line.resid1<<" name "<<rdc_line.name1<<" and ";
				fout<<" res "<<rdc_line.resid2<<" name "<<rdc_line.name2<<" out of range of pdb"<<endl;
				exp.at(i).erase(exp.at(i).begin()+j);
			}
		}
	}
	return;
}


int CRdc::load(string filename)
{
	int i,j;
	ifstream fin(filename.c_str());
	string line;
	vector<string> block;
	bool bstart,inloop,begin;

	exp.clear();

	bstart=0;
	inloop=0;
	begin=0;

	block.clear();
	while(getline(fin,line))
	{
		if(line.size()<4)
			continue;

		if(begin==1 &&  line.find("stop_")!=string::npos)
		{
			bstart=begin=0;
			actualload(&block);
			block.clear();
		}

		if(begin==1)
			block.push_back(line);


		if(line.find("loop_")!=string::npos)
			inloop=1;
		if(line.find("stop_")!=string::npos)
			inloop=0;

		if(inloop==0 && line.find("save_CNS/XPLOR_dipolar_coupling_")!=string::npos)
			bstart=1;

		if(bstart==1 && inloop==1 && line.find("_RDC_constraint.RDC_constraint_list_ID")!=string::npos)
			begin=1;
	}


	string nmrseq;
	char c1,c2;
	vector<int> out;
	int adj,adj2;

	for(i=0;i<exp.size();i++)
	{
		for(j=0;j<exp.at(i).size();j++)
		{
			c1=Sequence::name2code(exp.at(i).at(j).resname1);
			c2=Sequence::name2code(exp.at(i).at(j).resname2);
			exp.at(i).at(j).resid1--;
			exp.at(i).at(j).resid2--;
			if(exp.at(i).at(j).resid1+1>nmrseq.size())
				nmrseq.resize(exp.at(i).at(j).resid1+1,'U');
			nmrseq.at(exp.at(i).at(j).resid1)=c1;
			if(exp.at(i).at(j).resid2+1>nmrseq.size())
				nmrseq.resize(exp.at(i).at(j).resid2+1,'U');
			nmrseq.at(exp.at(i).at(j).resid2)=c2;
		}
	}

	out=Sequence::align(pdb->getseq(),nmrseq);
	adj=0;adj2=0;
	for(i=out.size()/5;i<out.size()-out.size()/5;i++)
	{
		if(out.at(i)!=0)
		{
			adj+=(out.at(i)-i);
			adj2++;
		}
	}
	adj=((float)adj/adj2+0.5);
	adj--;
	
	for(i=0;i<exp.size();i++)
	for(j=0;j<exp.at(i).size();j++)
	{
		exp.at(i).at(j).resid1-=adj;
		exp.at(i).at(j).resid2-=adj;
	}

	fillin();
	return adj+1;
}


void CRdc::actualload(vector<string> *block)
{
	istringstream iss;
	int i;
	string p;
	struct exp_rdc rdc_line;
	vector<struct exp_rdc> rdcs;
	string oldname1,oldname2;
	bool bfirst;

	oldname1="non";
	oldname2="non";


	bfirst=1;
	rdcs.clear();
	for(i=0;i<block->size();i++)
	{
		iss.clear();
		iss.str(block->at(i));
		iss>>p>>p>>p;
		iss>>rdc_line.resid1;
		iss>>rdc_line.resname1;
		iss>>rdc_line.name1;
		iss>>p>>p;
		iss>>rdc_line.resid2;
		iss>>rdc_line.resname2;
		iss>>rdc_line.name2;
		iss>>rdc_line.vrdc;
		
		if(bfirst==1)
		{
			bfirst=0;
			oldname1=rdc_line.name1;
			oldname2=rdc_line.name2;
		}

		if(rdc_line.name1.compare(oldname1)!=0 || rdc_line.name2.compare(oldname2)!=0)
		{
			if(rdcs.size()>0)
				exp.push_back(rdcs);
			rdcs.clear();
			oldname1=rdc_line.name1;
			oldname2=rdc_line.name2;
		}
		
		rdcs.push_back(rdc_line);
	}
	if(rdcs.size()>0)
		exp.push_back(rdcs);

	return;
}


void CRdc::load(int adjin)
{

	ifstream fin("exp_rdc.dat");

	if (!fin.is_open())
		return;

	istringstream iss;
	string line,p;
	int set;
	int nexp;
	struct exp_rdc rdc_line;
	struct noeatoms noe1,noe2;


	exp.clear();
	adj=adjin;

	getline(fin,line);
	iss.clear();
	iss.str(line);
	iss>>nexp;
	exp.resize(nexp);

	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		iss>>p;

		if(p.compare("set")==0)
		{
			iss>>set;
			continue;
		}

		while(p.compare("resid")!=0)
		{
			iss>>p;
		}
		iss>>rdc_line.resid1;

		while(p.compare("name")!=0)
		{
			iss>>p;
		}
		iss>>rdc_line.name1;

		getline(fin,line);
		iss.clear();
		iss.str(line);
		iss>>p;
		while(p.compare("resid")!=0)
		{
			iss>>p;
		}
		iss>>rdc_line.resid2;

		while(p.compare("name")!=0)
		{
			iss>>p;
		}
		iss>>rdc_line.name2;

		iss>>p;
		iss>>rdc_line.vrdc;

		rdc_line.resid1+=adj;
		rdc_line.resid2+=adj;

		if(abs(rdc_line.resid1-rdc_line.resid2)<2 && rdc_line.resid1-1>=0 && rdc_line.resid2-1>=0 && rdc_line.resid1-1<nres && rdc_line.resid2-1<nres) 
		{
			noe1=pdb->query(rdc_line.resid1,rdc_line.name1);
			noe2=pdb->query(rdc_line.resid2,rdc_line.name2);
			if(noe1.atoms.size()>0 && noe2.atoms.size()>0)
			{
				rdc_line.atomid1=noe1.atoms.at(0).at(0);
				rdc_line.atomid2=noe2.atoms.at(0).at(0);
				exp.at(set-1).push_back(rdc_line);
			}
			else
			{
				fout<<"In CRdc, exp set "<<set<<" res "<<rdc_line.resid1<<" name "<<rdc_line.name1<<" and ";
				fout<<" res "<<rdc_line.resid2<<" name "<<rdc_line.name2<<" cannot match pdb"<<endl;
			}
		}
		else
		{
				fout<<"In CRdc, exp set "<<set<<" res "<<rdc_line.resid1<<" name "<<rdc_line.name1<<" and ";
				fout<<" res "<<rdc_line.resid2<<" name "<<rdc_line.name2<<" outof pdb range"<<endl; 
		}

	}

	

	return;
}



CRdc::CRdc() 
{	
	fout.open("rdc_out.dat");
};
CRdc::~CRdc() 
{
	fout.close();
};



////S2 order parameter parts
//contact model based prediciton of methyl S2 oder parameters
void CS2::methyl_contact(string outfile)
{
	float a,b,c;
	float temp;
	int i;
	int begin,stop;
	vector<int> heavy,boundary,heavy2;
	vector<int> tempatom;
	vector<struct proton> protons;
	vector<float> result;
	int id;
	float table[11]={0,1,2,4,0,2,2,3,3,2,3};
	ofstream fout(outfile.c_str());

	a=0.26;
	b=2.2;
	c=0.125;

	pdb->proton(&protons,0);
	heavy=pdb->getheavy();
	boundary=pdb->getboundary();

	for(i=0;i<(int)protons.size();i++)
	{
		id=protons.at(i).id;

		if(id>=2)
			begin=boundary.at(id-2);
		else
			begin=0;
		stop=boundary.at(id-1);

		heavy2=heavy;
		heavy2.erase(heavy2.begin()+begin,heavy2.begin()+stop);

		result.clear();
		tempatom.clear();
		tempatom.push_back(protons.at(i).cpos);
		traj->get_contact(3.40,0.0,tempatom,heavy2,&result);
		temp=result.at(0)*a/pow(table[protons.at(i).type],b);
		temp=(1-exp(-2.0*temp))/(1+exp(-2.0*temp))-c;
		fout<<protons.at(i).id<<" "<<protons.at(i).code<<" "<<protons.at(i).name<<" "<<temp<<endl;
	}
	fout.clear();
}

//contact model based prediction of backbone NH order parameters.
void CS2::contact_model(string outfile, string expfile)
{
	unsigned int i,j;
	int id;
	vector<struct nh_group> tt;
	vector<struct co_group> tt2;
	vector<int> index;
	vector<int> heavy;
	vector<int> heavy2;
	vector<int> boundary;
	vector<float> result,result2;
	istringstream iss;
	bool bexp;
	string line,p;
	vector<string> part;
	float temp;
	struct bbs2_pre pre;
	vector<struct bbs2_pre> pres;
	vector<int> tempatom;

	pdb->bbnh(&tt);
	pdb->bbco(&tt2);
	heavy=pdb->getheavy();
	boundary=pdb->getboundary();


	pres.clear();
	for(i=0;i<tt.size();i++)
	{
		id=tt.at(i).id;
		if(id<2 || id>boundary.size()-1)
			continue;
		pre.id=id;
		pre.h=tt.at(i).hpos;
		for(j=0;j<tt2.size();j++)
		{
			if(tt2.at(j).id+1==pre.id)
				pre.o=tt2.at(j).opos;
		}	
		if(id>=3) pre.self_begin=boundary.at(id-3);
		else pre.self_begin=0;
		pre.self_stop=boundary.at(id-1);
		pres.push_back(pre);
	}
	

	for(i=0;i<pres.size();i++)
	{
		pre=pres.at(i);
		tempatom.clear();
		tempatom.push_back(pre.h);
		tempatom.push_back(pre.o);
		heavy2=heavy;
		heavy2.erase(heavy2.begin()+pre.self_begin,heavy2.begin()+pre.self_stop);
		result.clear();
		traj->get_contact(1.00,0.0,tempatom,heavy2,&result);

		temp=(result.at(0)*0.8+result.at(1))*2.656;
		temp=(1-exp(-2.0*temp))/(1+exp(-2.0*temp))-0.1;
		pres.at(i).s2=temp;
		pres.at(i).experiment=0.0;
	}

	//read in experimental data
	ifstream fin(expfile.c_str());
	if(fin.is_open())
	{
		bexp=1;
		while(getline(fin,line))
		{
			p=line.substr(0,1);
			if(p=="#")
				continue;

			iss.clear();
			iss.str(line);
			part.clear();
			while(iss>>p)
			{
				part.push_back(p);
			}

			id=atoi(part.at(0).c_str());
			temp=atof(part.at(1).c_str());

			for(i=0;i<pres.size();i++)
			{
				if(pres.at(i).id==id)
					pres.at(i).experiment=temp;
			}
		}
		fin.close();
	}
	else
		bexp=0;



	ofstream fout(outfile.c_str());
	for(i=0;i<pres.size();i++)
	{
		fout<<pres.at(i).id<<" ";
		fout<<pres.at(i).s2;
		if(bexp)
			fout<<" "<<pres.at(i).experiment;
		fout<<endl;
	}
	fout.close();

	return;
};




void CS2::init(int site)
{
	int i;
	string name;
	FILE *fexp;

	ns2=0;
	nh.clear();
	xx.clear();
	yy.clear();
	zz.clear();

	flag=site;

	if(site==1)  //backbone NH
	{
		vector<struct nh_group> tt;
		pdb->bbnh(&tt);
		ns2=tt.size();
		nh.resize(ns2);
		for(i=0;i<ns2;i++)
		{
			nh.at(i).x1=tt.at(i).id;
			nh.at(i).x2=tt.at(i).npos;
			nh.at(i).x3=tt.at(i).hpos;
		}
	}
	else if(site==2) //backbone Ca-Ha
	{
		pdb->caha(&nh);
		ns2=nh.size();
	}

	else
	{
		flag=3;
		pdb->ired(&red);
		ns2=red.size();
		nh.resize(ns2);
		for(i=0;i<ns2;i++)
		{
			nh.at(i).x1=red.at(i).id;
			nh.at(i).x2=red.at(i).index1;
			nh.at(i).x3=red.at(i).index2;
		}
	}

	traj->getvector(nh,&xx,&yy,&zz);
	

	name="block.dat";
	fexp=fopen(name.c_str(),"rt");
	if(fexp!=NULL)
	{
		fscanf(fexp,"%d",&block_size);
		fclose(fexp);
	}
	else
		block_size=nconf;
	block_number=nconf/block_size;
	printf("block_size is %d, block_number is %d\n",block_size,block_number);
	
	
	return;
}


void CS2::doit()		
{
	int i,j,jj;
	int b,block_b,block_s;
	float t;
	float **u; 
	float *d,*e;
	int shift;
		
	u=matrix(1,ns2,1,ns2);
	d=myvector(1,ns2);	
	e=myvector(1,ns2);

	pre.clear();
	pre.resize(ns2);


	for(b=0;b<block_number;b++) //average within preset time window !
	{
		block_b=b*block_size;
		block_s=block_b+block_size;
		
		for(i=1;i<=ns2;i++)
		{
			for(j=1;j<=ns2;j++)
				u[i][j]=0.0;
		}		
				
		for(i=block_b;i<block_s;i++)
		{
			shift=i*ns2;
			for(j=0;j<ns2;j++)
			for(jj=0;jj<ns2;jj++)
			{
				t=xx[shift+j]*xx[shift+jj];
				t+=yy[shift+j]*yy[shift+jj];
				t+=zz[shift+j]*zz[shift+jj];
				u[j+1][jj+1]+=t*t*1.5-0.5;
			}	
		}
		for(j=0;j<ns2;j++)
		for(jj=0;jj<ns2;jj++)
		{
			u[j+1][jj+1]/=block_size;
		}

		/*
		ofstream fout("temp_matrix.dat");
		for(j=0;j<ns2;j++)
		{
			for(jj=0;jj<ns2;jj++)
				fout<<u[j+1][jj+1]<<" ";
			fout<<endl;
		}
		fout.close();
		*/
		
		
		tred2(u,ns2,d,e);
		tqli(d,e,ns2,u); //at this point, d is eigen values and u is the eigen vectors.

		/*
		ofstream fout2("temp_value.dat");
		for(jj=1;jj<=ns2;jj++)
			fout<<d[jj]<<endl;
		fout2.close();
		*/

		for(j=1;j<=ns2;j++)
		{
			t=0;
			for(jj=1;jj<=ns2-5;jj++)
				t+=d[jj]*u[j][jj]*u[j][jj];
			pre.at(j-1)+=(1-t);
		}		
	}

	for(j=0;j<ns2;j++)
	{
		pre.at(j)/=block_number;
		if(flag==3) red.at(j).s2.pre=pre.at(j);
	}


	free_vector(d,1,ns2);
	free_vector(e,1,ns2);
	free_matrix(u,1,ns2,1,ns2);

	if(flag==3)  pdb->loadred(&red);

	return;
}

void CS2::print(string name)
{
	int i;
	ofstream fout(name.c_str());

	for(i=0;i<ns2;i++)
	{
		if(flag==1) fout<<nh.at(i).x1<<" N H "<<pre.at(i)<<endl;
		else if(flag==2) fout<<nh.at(i).x1<<" CA HA "<<pre.at(i)<<endl;
		else if(flag==3) fout<<red.at(i).id<<" "<<red.at(i).code<<" "<<red.at(i).s2.name1<<" "<<red.at(i).s2.name2<<" "<<red.at(i).s2.pre<<endl;
	}
	fout.close();

	return;
}


CS2::CS2() {};
CS2::~CS2() {};
