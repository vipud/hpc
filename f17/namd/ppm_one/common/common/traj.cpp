#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
using namespace std;

#include "traj.h"
#include <omp.h>
#include <openacc.h>
using namespace ldw_math;

#pragma acc routine seq
	int my_cross(double z[3],double x[3],double y[3])
	{
			z[0]=x[1]*y[2]-x[2]*y[1];
			z[1]=-x[0]*y[2]+x[2]*y[0];
			z[2]=x[0]*y[1]-x[1]*y[0];
			return 0;
	}

#pragma acc routine seq
	double my_dot(double x[3],double y[3])
	{
			return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
	}

#pragma acc routine seq
	double my_mysign(double a,double b)
	{
			double result;

			a=fabs(a);
			if(b>=0)
			{
					result=a;
			}
			else
			{
					result=-a;
			}

			return result;
	}

#pragma acc routine seq
	double my_mymax(double a, double b)
	{
			double result;
			if(a>b)
			{
					result=a;
			}
			else
			{
					result=b;
			}
			return result;
	}

#pragma acc routine seq 
	double my_PYTHAG(double a, double b)
	{
		double at = fabs(a), bt = fabs(b), ct, result;

		if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
		else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
		else result = 0.0;
		return(result);
	}

#pragma acc routine seq
int my_dsvd(double a[6][3], int m, int n, double w[3], double v[3][3])
	{

		int flag, i, its, j, jj, k, l, nm;
		double c, f, h, s, x, y, z;
		double anorm = 0.0, g = 0.0, scale = 0.0;
		double rv1[3];
	  
		//if (m < n) 
		//{
			//fprintf(stderr, "#rows must be > #cols \n");
		//	return 0;
		//}
	  
		//rv1 = (double *)malloc((unsigned int) n*sizeof(double));
		//rv1 = new double[n];

	
		for (i = 0; i < n; i++) 
		{
			
			l = i + 1;
			rv1[i] = scale * g;
			g = 0.0;
			s = 0.0;
			scale = 0.0;
			if (i < m) 
			{
				for (k = i; k < m; k++) 
					scale += fabs((double)a[k][i]);
				if (scale) 
				{
					for (k = i; k < m; k++) 
					{
						a[k][i] = (double)((double)a[k][i]/scale);
						s += ((double)a[k][i] * (double)a[k][i]);
					}
					f = (double)a[i][i];
	                
					g = -my_mysign(sqrt(s), f);
	                
					h = f * g - s;
					a[i][i] = (double)(f - g);
					if (i != n - 1) 
					{
						for (j = l; j < n; j++) 
						{
							for (s = 0.0, k = i; k < m; k++) 
								s += ((double)a[k][i] * (double)a[k][j]);
							f = s / h;
							for (k = i; k < m; k++) 
								a[k][j] += (double)(f * (double)a[k][i]);
						}
					}
					for (k = i; k < m; k++) 
						a[k][i] = (double)((double)a[k][i]*scale);
				}
			}
			w[i] = (double)(scale * g);
	    
	
			g = 0.0;
			s = 0.0;
			scale = 0.0;
			if (i < m && i != n - 1) 
			{
				for (k = l; k < n; k++) 
					scale += fabs((double)a[i][k]);
				if (scale) 
				{
					for (k = l; k < n; k++) 
					{
						a[i][k] = (double)((double)a[i][k]/scale);
						s += ((double)a[i][k] * (double)a[i][k]);
					}
					f = (double)a[i][l];
					g = -my_mysign(sqrt(s), f);
					h = f * g - s;
					a[i][l] = (double)(f - g);
					for (k = l; k < n; k++) 
						rv1[k] = (double)a[i][k] / h;
					if (i != m - 1) 
					{
						for (j = l; j < m; j++) 
						{
							for (s = 0.0, k = l; k < n; k++) 
								s += ((double)a[j][k] * (double)a[i][k]);
							for (k = l; k < n; k++) 
								a[j][k] += (double)(s * rv1[k]);
						}
					}
					for (k = l; k < n; k++) 
						a[i][k] = (double)((double)a[i][k]*scale);
				}
			}
			anorm = my_mymax(anorm, (fabs((double)w[i]) + fabs(rv1[i])));


		}
	

		for (i = n - 1; i >= 0; i--) 
		{
			if (i < n - 1) 
			{
				if (g) 
				{
					for (j = l; j < n; j++)
						v[j][i] = (double)(((double)a[i][j] / (double)a[i][l]) / g);
						
					for (j = l; j < n; j++) 
					{
						for (s = 0.0, k = l; k < n; k++) 
							s += ((double)a[i][k] * (double)v[k][j]);
						for (k = l; k < n; k++) 
							v[k][j] += (double)(s * (double)v[k][i]);
					}
				}
				for (j = l; j < n; j++) 
					v[i][j] = v[j][i] = 0.0;
			}
			v[i][i] = 1.0;
			g = rv1[i];
			l = i;
		}

		for (i = n - 1; i >= 0; i--) 
		{
			l = i + 1;
			g = (double)w[i];
			if (i < n - 1) 
				for (j = l; j < n; j++) 
					a[i][j] = 0.0;
			if (g) 
			{
				g = 1.0 / g;
				if (i != n - 1) 
				{
					for (j = l; j < n; j++) 
					{
						for (s = 0.0, k = l; k < m; k++) 
							s += ((double)a[k][i] * (double)a[k][j]);
						f = (s / (double)a[i][i]) * g;
						for (k = i; k < m; k++) 
							a[k][j] += (double)(f * (double)a[k][i]);
					}
				}
				for (j = i; j < m; j++) 
					a[j][i] = (double)((double)a[j][i]*g);
			}
			else 
			{
				for (j = i; j < m; j++) 
					a[j][i] = 0.0;
			}
			++a[i][i];
		}


		for (k = n - 1; k >= 0; k--) 
		{                             
			for (its = 0; its < 30; its++) 
			{                        
				flag = 1;
				for (l = k; l >= 0; l--) 
				{                 
					nm = l - 1;
					if (fabs(rv1[l]) + anorm == anorm) 
					{
						flag = 0;
						break;
					}
					if (fabs((double)w[nm]) + anorm == anorm) 
						break;
				}
				if (flag) 
				{
					c = 0.0;
					s = 1.0;
					for (i = l; i <= k; i++) 
					{
						f = s * rv1[i];
						if (fabs(f) + anorm != anorm) 
						{
							g = (double)w[i];
							h = my_PYTHAG(f, g);
							w[i] = (double)h; 
							h = 1.0 / h;
							c = g * h;
							s = (- f * h);
							for (j = 0; j < m; j++) 
							{
								y = (double)a[j][nm];
								z = (double)a[j][i];
								a[j][nm] = (double)(y * c + z * s);
								a[j][i] = (double)(z * c - y * s);
							}
						}
					}
				}
		
				z = (double)w[k];
				if (l == k) 
				{                  
					if (z < 0.0) 
					{              
						w[k] = (double)(-z);
						for (j = 0; j < n; j++) 
							v[j][k] = (-v[j][k]);
					}
					//break;
					its = 30;
				} else {
				//if (its >= 30) {
				//	free((void*) rv1);
				//	fprintf(stderr, "No convergence after 30,000! iterations \n");
				//	return(0);
				//}
	    
				
				x = (double)w[l];
				nm = k - 1;
				y = (double)w[nm];
				g = rv1[nm];
				h = rv1[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = my_PYTHAG(f, 1.0);
				f = ((x - z) * (x + z) + h * ((y / (f + my_mysign(g, f))) - h)) / x;
	          
				
				c = s = 1.0;
				for (j = l; j <= nm; j++) 
				{
	
					i = j + 1;
					g = rv1[i];
					y = (double)w[i];
					h = s * g;
					g = c * g;
					z = my_PYTHAG(f, h);
					rv1[j] = z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y = y * c;
					for (jj = 0; jj < n; jj++) 
					{
						x = (double)v[jj][j];
						z = (double)v[jj][i];
						v[jj][j] = (double)(x * c + z * s);
						v[jj][i] = (double)(z * c - x * s);
					}
					z = my_PYTHAG(f, h);
					w[j] = (double)z;
					if (z) 
					{
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = (c * g) + (s * y);
					x = (c * y) - (s * g);
					for (jj = 0; jj < m; jj++) 
					{
						y = (double)a[jj][j];
						z = (double)a[jj][i];
						a[jj][j] = (double)(y * c + z * s);
						a[jj][i] = (double)(z * c - y * s);
					}
		
				}
		
				rv1[l] = 0.0;
				rv1[k] = f;
				w[k] = (double)x;
				} // end shitty else
			}
		}
		//delete(rv1);
		//free((void*) rv1);

		return 1;
	}


#pragma acc routine seq
	void my_ring(double x[6][3], int m, double ori[3])
	{
			int i,j;
			int id;
			double w[3];
			double v[3][3];
			double d;
			double xx[6][3];

			for(i=0;i<m;i++)
			for(j=0;j<3;j++)
					xx[i][j]=x[i][j];


			my_dsvd(xx, m, 3, w,v);

			d=w[0];id=0;
			if(w[1]<d)
			{
					d=w[1];
					id=1;
			}
			if(w[2]<d)
			{
					d=w[2];
					id=2;
			}

			for(i=0;i<3;i++)
					ori[i]=v[i][id];

			return;
	}

#pragma acc routine seq
	double my_area( double a, double b, double c )
	{
	  double s;
	  double y;
	  s = (a + b + c)/2;
	  s =  s * (s - a)*(s - b)*(s - c);
	  if(s<0)
		y=0;
	  else
		y = sqrt( s );
	  return y;
	}

#pragma acc routine seq
	double my_veclength(double x[3])
	{
			return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	}

#pragma acc routine seq
	void my_project(double ori[3], double p1[3], double p2[3])
	{
			double d;
			int i;

			d=ori[0]*p1[0]+ori[1]*p1[1]+ori[2]*p1[2];
			for(i=0;i<3;i++)
			{
					p2[i]=p1[i]-d*ori[i];
			}
			return;
	}

#pragma acc routine seq
	double my_effect(double x[6][3], int m, double ori[3], double p1[3])
	{
			int i,j;
			double t1[3];
			double t2[3];
			double t3[3];
			double t4[3];
			double tt,d;
			double s,ss;
			double leg1,leg2,leg3;
			double p2[3];

			ss=0;
			my_project(ori,p1,p2);
			for(i=0;i<m-1;i++)
			{
					for(j=0;j<3;j++)
					{
							t1[j]=x[i][j]-p1[j];
							t2[j]=x[i+1][j]-x[i][j];
							t3[j]=x[i+1][j]-p1[j];
					}
					leg1=my_veclength(t1);
					leg3=my_veclength(t3);
					d=1/(leg1*leg1*leg1)+1/(leg3*leg3*leg3);
					my_cross(t4,t2,t1);
					tt=my_dot(t4,ori);

					for(j=0;j<3;j++)
					{
							t1[j]=x[i][j]-p2[j];
							t3[j]=x[i+1][j]-p2[j];
					}
					leg1=my_veclength(t1);
					leg2=my_veclength(t2);
					leg3=my_veclength(t3);
					s=my_area(leg1,leg2,leg3);
					if(tt<0)
							s=-s;
					ss+=s*d;
			}

			//special case
			{
					i=m-1;
					for(j=0;j<3;j++)
					{
							t1[j]=x[i][j]-p1[j];
							t2[j]=x[0][j]-x[i][j];
							t3[j]=x[0][j]-p1[j];
					}
					leg1=my_veclength(t1);
					leg3=my_veclength(t3);
					d=1/(leg1*leg1*leg1)+1/(leg3*leg3*leg3);
					my_cross(t4,t2,t1);
					tt=my_dot(t4,ori);

					for(j=0;j<3;j++)
					{
							t1[j]=x[i][j]-p2[j];
							t3[j]=x[0][j]-p2[j];
					}
					leg1=my_veclength(t1);
					leg2=my_veclength(t2);
					leg3=my_veclength(t3);
					s=my_area(leg1,leg2,leg3);
					if(tt<0)
							s=-s;
					ss+=s*d;
			}

			return ss;
	}

	#pragma acc routine seq
	double my_coor_to_angle(double x2,double y2,double z2,double x3,double y3,double z3,double x4,double y4,double z4)
	{
			double angle;
			double b[3],c[3];

	        
			b[0]=x3-x2;b[1]=y3-y2;b[2]=z3-z2;
			c[0]=x4-x3;c[1]=y4-y3;c[2]=z4-z3;
			angle=my_dot(b,c)/my_veclength(b)/my_veclength(c);
			if(angle>1.0) angle=1.0;
			if(angle<-1.0) angle=-1.0;
	/*        angle=sqrt(1-angle*angle);*/
			return angle;
	}




// ///////////////////////////////////////////////////////////////////////////////////////////////////////////

void CTraj::clear()
{
	x.clear();
	y.clear();
	z.clear();
	nframe=0;
}




int CTraj::appendcoor(string filename)
{
	double xx,yy,zz;
	string line,part;

	ifstream fin(filename.c_str());
	
	while(getline(fin,line))
	{
		if(part=="ENDMDL" || part=="END")
		{
			if(x.size()%natom!=0)
				cout<<"In traj reading, suppose to read "<<natom<<" coors but actually read in "<<x.size()<<endl;
		}

		part=line.substr(0,6);
		if(part!="ATOM  " && part!="HETATM")
			continue;

		part=line.substr(30,8);
		xx=atof(part.c_str());
		x.push_back(xx);
		part=line.substr(38,8);
		yy=atof(part.c_str());
		y.push_back(yy);
		part=line.substr(46,8);
		zz=atof(part.c_str());
		z.push_back(zz);
	};


	nframe=x.size()/natom;
	return nframe;
}


// Updated function for OpenACC
// Includes data directives
// Copy x,y,z arrays into GPU
int CTraj::loadcoor(string filename)
{
	double xx,yy,zz;
	string line,part;
	bool bend;

	bend=0;


	ifstream fin(filename.c_str());

	
	while(getline(fin,line))
	{
		if(part=="ENDMDL" || part=="END")
		{
			//finished read first molecule 
			bend=1;
			if(x.size()%natom!=0)
				cout<<"In traj reading, suppose to read "<<natom<<" coors but actually read in "<<x.size()<<endl;
		}

		part=line.substr(0,6);
		if(part!="ATOM  " && part!="HETATM")
			continue;

		part=line.substr(30,8);
		xx=atof(part.c_str());
		x.push_back(xx);
		part=line.substr(38,8);
		yy=atof(part.c_str());
		y.push_back(yy);
		part=line.substr(46,8);
		zz=atof(part.c_str());
		z.push_back(zz);
		if(bend==0)
		{
			part=line.substr(12,4);
			atomname.push_back(part);
		}
	};

	if(natom==0)
		nframe=0;
	else
		nframe=x.size()/natom;

	// New variables to storing underlying x,y,z arrays for better GPU support
	x_arr = x.data();
	x_size = x.size();
	y_arr = y.data();
	y_size = y.size();
	z_arr = z.data();
	z_size = z.size();

#pragma acc enter data copyin(this)
#pragma acc enter data copyin( x_arr[0:x_size], y_arr[0:y_size], z_arr[0:z_size])

	return nframe;
}

int CTraj::set_range(int begin,int stop)
{
	if(begin<0)
	{
		begin=0;
		cout<<"Reset begin to "<<begin<<endl;
	}

	if(stop>nframe)
	{
		stop=nframe;
		cout<<"Reset stop to "<<stop<<endl;
	}

	if(stop==0)
		stop=nframe;

	nframe=stop-begin;


	if(stop<nframe)
	{
		x.erase(x.begin()+stop*natom,x.begin()+x.size());
		y.erase(y.begin()+stop*natom,y.begin()+y.size());
		z.erase(z.begin()+stop*natom,z.begin()+z.size());
	}
	if(begin>0)
	{
		x.erase(x.begin(),x.begin()+begin*natom);
		y.erase(y.begin(),y.begin()+begin*natom);
		z.erase(z.begin(),z.begin()+begin*natom);
	}

	cout<<"Using frames from "<<begin<<" to "<<stop<<endl;
	return nframe;
}



void CTraj::getschbond2(vector<struct proton> *protons, vector<struct bbhbond_group> *bb, vector< vector<ehbond> > *effect, vector< vector<eschbond> > *effect_sc)
{
	int i,j,k;
	int base;
	int n,h[3],c,o;
	int nh;
	int nid,cid;
	double u[3];
	float x2,x3,x4,x5;
	float y2,y3,y4,y5;
	float z2,z3,z4,z5;
	float phi,psi,d;

	int kk,good_kk;
	float min_d;

	effect->resize(nres);
	for(i=0;i<(int)effect->size();i++)
	{
		effect->at(i).resize(nframe);
		for(k=0;k<nframe;k++)
		{
			effect->at(i).at(k).n_length=0;
			effect->at(i).at(k).c_length=0;
			effect->at(i).at(k).n_phi=0;
			effect->at(i).at(k).c_phi=0;
			effect->at(i).at(k).n_psi=0;
			effect->at(i).at(k).c_psi=0;
		}
	}

	effect_sc->resize(protons->size());
	for(i=0;i<(int)effect_sc->size();i++)
	{
		effect_sc->at(i).resize(nframe);
		for(k=0;k<nframe;k++)
		{
			effect_sc->at(i).at(k).n_length=0;
			effect_sc->at(i).at(k).c_length=0;
			effect_sc->at(i).at(k).n_phi=0;
			effect_sc->at(i).at(k).c_phi=0;
			effect_sc->at(i).at(k).n_psi=0;
			effect_sc->at(i).at(k).c_psi=0;
		}
	}


	
	for(i=0;i<(int)protons->size();i++)
	{
		for(j=0;j<(int)bb->size();j++)
		{
			nh=protons->at(i).nh;
			for(k=0;k<nh;k++)
				h[k]=protons->at(i).hpos[k];
			n=protons->at(i).cpos;
			if(h[0]<=-1 || n<=-1 || (nh==2 && h[1]<0) || (nh==3 && h[2]<0) ) 
				continue;
			c=bb->at(j).cpos;
			o=bb->at(j).opos;
			if(o<=-1 || c<=-1)
				continue;
			n--;c--;o--;
			for(k=0;k<nh;k++)
				h[k]--;

			for(k=0;k<nframe;k++)
			{
				//cout<<"i,j,k is "<<i<<" "<<j<<" "<<k<<endl;
				min_d=100000.0;
				for(kk=0;kk<nh;kk++)
				{
					base=k*natom;
					u[0]=x[h[kk]+base]-x[o+base];
					u[1]=y[h[kk]+base]-y[o+base];
					u[2]=z[h[kk]+base]-z[o+base];
					d=veclength(u);
					if(d<min_d)
					{
						min_d=d;
						good_kk=kk;
					}
				}
				d=min_d;
				x2=x[n+base];y2=y[n+base];z2=z[n+base];
				x3=x[h[good_kk]+base];y3=y[h[good_kk]+base];z3=z[h[good_kk]+base];
				x4=x[o+base];y4=y[o+base];z4=z[o+base];
				x5=x[c+base];y5=y[c+base];z5=z[c+base];
				phi=coor_to_angle(x2,y2,z2,x3,y3,z3,x4,y4,z4);
				psi=coor_to_angle(x3,y3,z3,x4,y4,z4,x5,y5,z5);
				if(d<3 && phi>0.5 && psi>0.5)
				{
					d=1/(d-1);
					phi*=phi;
					psi*=psi;
					nid=protons->at(i).id-1;
					cid=bb->at(j).id-1; /* C start from 0*/
					
					//effect->at(nid).at(k).n_length=d;
					//effect->at(nid).at(k).n_phi=phi;
					//effect->at(nid).at(k).n_psi=psi;

					effect_sc->at(i).at(k).n_length=d;
					effect_sc->at(i).at(k).n_phi=phi;
					effect_sc->at(i).at(k).n_psi=psi;
					effect_sc->at(i).at(k).id=protons->at(i).id;
					effect_sc->at(i).at(k).type=protons->at(i).type;

					if(bb->at(j).type==1)
					{
						effect->at(cid).at(k).c_length=d;
						effect->at(cid).at(k).c_phi=phi;
						effect->at(cid).at(k).c_psi=psi;
					}
					
				}
			}//for k
		}//for j
	}//for i

	return;
}


// New funtion for OpenACC
void CTraj::gethbond_acc(bbhbond_group *hbond, int _hbond_size, ehbond *effect_arr, int effect_size)
{
	double st = omp_get_wtime();
	int i,j,k;
	int base;
	int nid,cid;
	int n,h,c,o;
	//double u[3];
	double x2,x3,x4,x5,y2,y3,y4,y5,z2,z3,z4,z5;
	double d,phi,psi;

	//effect->resize(nres);
	//ehbond *effect_arr = effect->data();
	//int effect_size = effect->size();

	//double *x_arr = x_arr;
	//double *y_arr = y_arr;
	//double *z_arr = z_arr;
	//int x_size = x_size;
	//int y_size = y_size;
	//int z_size = z_size;

	//#pragma acc enter data copyin(effect_arr[0:effect_size])

#pragma acc data present(x_arr[0:x_size],y_arr[0:y_size],z_arr[0:z_size],effect_arr[0:effect_size],hbond[0:_hbond_size])
{

#pragma acc parallel
{
#pragma acc loop gang independent
	for(i=0;i<_hbond_size;i++)
	{
#pragma acc loop vector independent private(x2,x3,x4,x5,y2,y3,y4,y5,z2,z3,z4,z5,n,h,c,o,nid,cid,k,phi,psi)
		for(j=0;j<_hbond_size;j++)
		{
			k=j-i;
			if(k<3 && k>-3)
				continue;
			n=hbond[i].npos;
			h=hbond[i].hpos;
			if(h<=-1 || n<=-1) 
				continue;
			c=hbond[j].cpos;
			o=hbond[j].opos;
			if(o<=-1 || c<=-1)
				continue;
			n--;h--;c--;o--;
			if(h<0 || n<0 || o<0 || c<0)
				continue;

			double u[3];
#pragma acc loop seq independent
			for(k=0;k<nframe;k++)
			{
				//cout<<"i,j,k is "<<i<<" "<<j<<" "<<k<<endl;
				base=k*natom;
				u[0]=x_arr[h+base]-x_arr[o+base];
				u[1]=y_arr[h+base]-y_arr[o+base];
				u[2]=z_arr[h+base]-z_arr[o+base];
				d=my_veclength(u);
				x2=x_arr[n+base];y2=y_arr[n+base];z2=z_arr[n+base];
				x3=x_arr[h+base];y3=y_arr[h+base];z3=z_arr[h+base];
				x4=x_arr[o+base];y4=y_arr[o+base];z4=z_arr[o+base];
				x5=x_arr[c+base];y5=y_arr[c+base];z5=z_arr[c+base];
				phi=my_coor_to_angle(x2,y2,z2,x3,y3,z3,x4,y4,z4);
				psi=my_coor_to_angle(x3,y3,z3,x4,y4,z4,x5,y5,z5);
				if(d<3 && phi>0.5 && psi>0.5)
				{
					d=1/(d-1);
					phi*=phi;
					psi*=psi;
					nid=hbond[i].id-1;
					cid=hbond[j].id-1; /* C start from 0*/
					if(hbond[i].type==1)
					{
						#pragma acc atomic update
						effect_arr[nid].n_length+=d;
						#pragma acc atomic update
						effect_arr[nid].n_phi+=phi;
						#pragma acc atomic update
						effect_arr[nid].n_psi+=psi;
					}
					if(hbond[j].type==1)
					{
						#pragma acc atomic update
						effect_arr[cid].c_length+=d;
						#pragma acc atomic update
						effect_arr[cid].c_phi+=phi;
						#pragma acc atomic update
						effect_arr[cid].c_psi+=psi;
					}
				}
			}
		}
	}
} // end parallel region
#pragma acc parallel loop
	for(i=0;i<effect_size;i++)
	{
		effect_arr[i].n_length/=nframe;
		effect_arr[i].c_length/=nframe;
		effect_arr[i].n_phi/=nframe;
		effect_arr[i].c_phi/=nframe;
		effect_arr[i].n_psi/=nframe;
		effect_arr[i].c_psi/=nframe;
	}

} // END DATA REGION

	cout << "gethbond: " << omp_get_wtime() - st << " seconds" << endl;
	return;
}



//hbond type (1, bb) (12, sc OH) (13,sc NH) (22 or 23, sc CO)
void CTraj::gethbond(vector<bbhbond_group> *hbond,vector<ehbond> *effect)
{
	double st = omp_get_wtime();
	int i,j,k;
	int base;
	int nid,cid;
	int n,h,c,o;
	double u[3];
	double x2,x3,x4,x5,y2,y3,y4,y5,z2,z3,z4,z5;
	double d,phi,psi;

	
	effect->resize(nres);
	for(i=0;i<(int)hbond->size();i++)
	{

		for(j=0;j<(int)hbond->size();j++)
		{
			k=j-i;
			if(k<3 && k>-3)
				continue;
			n=hbond->at(i).npos;
			h=hbond->at(i).hpos;
			if(h<=-1 || n<=-1) 
				continue;
			c=hbond->at(j).cpos;
			o=hbond->at(j).opos;
			if(o<=-1 || c<=-1)
				continue;
			n--;h--;c--;o--;
			if(h<0 || n<0 || o<0 || c<0)
				continue;



			for(k=0;k<nframe;k++)
			{
				//cout<<"i,j,k is "<<i<<" "<<j<<" "<<k<<endl;
				base=k*natom;
				u[0]=x[h+base]-x[o+base];
				u[1]=y[h+base]-y[o+base];
				u[2]=z[h+base]-z[o+base];
				d=veclength(u);
				x2=x[n+base];y2=y[n+base];z2=z[n+base];
				x3=x[h+base];y3=y[h+base];z3=z[h+base];
				x4=x[o+base];y4=y[o+base];z4=z[o+base];
				x5=x[c+base];y5=y[c+base];z5=z[c+base];
				phi=coor_to_angle(x2,y2,z2,x3,y3,z3,x4,y4,z4);
				psi=coor_to_angle(x3,y3,z3,x4,y4,z4,x5,y5,z5);
				if(d<3 && phi>0.5 && psi>0.5)
				{
					d=1/(d-1);
					phi*=phi;
					psi*=psi;
					nid=hbond->at(i).id-1;
					cid=hbond->at(j).id-1; /* C start from 0*/
					if(hbond->at(i).type==1)
					{
						effect->at(nid).n_length+=d;
						effect->at(nid).n_phi+=phi;
						effect->at(nid).n_psi+=psi;
					}
					if(hbond->at(j).type==1)
					{
						effect->at(cid).c_length+=d;
						effect->at(cid).c_phi+=phi;
						effect->at(cid).c_psi+=psi;
					}
				}
			}
		}
	}

	for(i=0;i<(int)effect->size();i++)
	{
		effect->at(i).n_length/=nframe;
		effect->at(i).c_length/=nframe;
		effect->at(i).n_phi/=nframe;
		effect->at(i).c_phi/=nframe;
		effect->at(i).n_psi/=nframe;
		effect->at(i).c_psi/=nframe;
	}
	cout << "gethbond: " << omp_get_wtime() - st << " seconds" << endl;
	return;
}


//hbond type (1, bb) (12, sc OH) (13,sc NH) (22 or 23, sc CO)
void CTraj::gethbond(vector<bbhbond_group> *hbond,vector<ehbond> *effect, double cutoff)
{
	int i,j,k;
	int base;
	int nid,cid;
	int n,h,c,o;
	double u[3];
	double x2,x3,x4,x5,y2,y3,y4,y5,z2,z3,z4,z5;
	double d,phi,psi;

	
	effect->resize(nres);
	for(i=0;i<(int)hbond->size();i++)
	{

		for(j=0;j<(int)hbond->size();j++)
		{
			k=j-i;
			if(k<3 && k>-3)
				continue;
			n=hbond->at(i).npos;
			h=hbond->at(i).hpos;
			if(h<=-1 || n<=-1) 
				continue;
			c=hbond->at(j).cpos;
			o=hbond->at(j).opos;
			if(o<=-1 || c<=-1)
				continue;
			n--;h--;c--;o--;
			if(h<0 || n<0 || o<0 || c<0)
				continue;



			for(k=0;k<nframe;k++)
			{
				//cout<<"i,j,k is "<<i<<" "<<j<<" "<<k<<endl;
				base=k*natom;
				u[0]=x[h+base]-x[o+base];
				u[1]=y[h+base]-y[o+base];
				u[2]=z[h+base]-z[o+base];
				d=veclength(u);
				x2=x[n+base];y2=y[n+base];z2=z[n+base];
				x3=x[h+base];y3=y[h+base];z3=z[h+base];
				x4=x[o+base];y4=y[o+base];z4=z[o+base];
				x5=x[c+base];y5=y[c+base];z5=z[c+base];
				phi=coor_to_angle(x2,y2,z2,x3,y3,z3,x4,y4,z4);
				psi=coor_to_angle(x3,y3,z3,x4,y4,z4,x5,y5,z5);
				if(d<3 && phi>0.5 && psi>0.5)
				{
					d=1/(d-1);
					phi*=phi;
					psi*=psi;
					nid=hbond->at(i).id-1;
					cid=hbond->at(j).id-1; /* C start from 0*/

					if(d>cutoff)
						d=cutoff;

					if(hbond->at(i).type==1)
					{
						effect->at(nid).n_length+=d;
						effect->at(nid).n_phi+=phi;
						effect->at(nid).n_psi+=psi;
					}
					if(hbond->at(j).type==1)
					{
						effect->at(cid).c_length+=d;
						effect->at(cid).c_phi+=phi;
						effect->at(cid).c_psi+=psi;
					}
				}
			}
		}
	}

	for(i=0;i<(int)effect->size();i++)
	{
		effect->at(i).n_length/=nframe;
		effect->at(i).c_length/=nframe;
		effect->at(i).n_phi/=nframe;
		effect->at(i).c_phi/=nframe;
		effect->at(i).n_psi/=nframe;
		effect->at(i).c_psi/=nframe;
	}
	return;
}



void CTraj::gethbond(vector<bbhbond_group> *hbond,vector< vector<ehbond> > *effect)
{
	int i,j,k;
	int base;
	int nid,cid;
	int n,h,c,o;
	double u[3];
	double x2,x3,x4,x5,y2,y3,y4,y5,z2,z3,z4,z5;
	double d,phi,psi;

	effect->resize(nres);
	for(i=0;i<(int)effect->size();i++)
	{
		effect->at(i).resize(nframe);
		for(k=0;k<nframe;k++)
		{
			effect->at(i).at(k).n_length=0;
			effect->at(i).at(k).c_length=0;
			effect->at(i).at(k).n_phi=0;
			effect->at(i).at(k).c_phi=0;
			effect->at(i).at(k).n_psi=0;
			effect->at(i).at(k).c_psi=0;
		}
	}

	for(i=0;i<(int)hbond->size();i++)
	{

		for(j=0;j<(int)hbond->size();j++)
		{
			k=j-i;
			if(k<3 && k>-3)
				continue;
			n=hbond->at(i).npos;
			h=hbond->at(i).hpos;
			if(h<0) 
				continue;
			c=hbond->at(j).cpos;
			o=hbond->at(j).opos;
			if(o<0)
				continue;
			n--;h--;c--;o--;


			for(k=0;k<nframe;k++)
			{
				//cout<<"i,j,k is "<<i<<" "<<j<<" "<<k<<endl;
				base=k*natom;
				u[0]=x[h+base]-x[o+base];
				u[1]=y[h+base]-y[o+base];
				u[2]=z[h+base]-z[o+base];
				d=veclength(u);
				x2=x[n+base];y2=y[n+base];z2=z[n+base];
				x3=x[h+base];y3=y[h+base];z3=z[h+base];
				x4=x[o+base];y4=y[o+base];z4=z[o+base];
				x5=x[c+base];y5=y[c+base];z5=z[c+base];
				phi=coor_to_angle(x2,y2,z2,x3,y3,z3,x4,y4,z4);
				psi=coor_to_angle(x3,y3,z3,x4,y4,z4,x5,y5,z5);
				if(d<3 && phi>0.5 && psi>0.5)
				{
					d=1/(d-1);
					phi*=phi;
					psi*=psi;
					nid=hbond->at(i).id-1;
					cid=hbond->at(j).id-1; /* C start from 0*/
					if(hbond->at(i).type==1)
					{
						effect->at(nid).at(k).n_length+=d;
						effect->at(nid).at(k).n_phi+=phi;
						effect->at(nid).at(k).n_psi+=psi;
					}
					if(hbond->at(j).type==1)
					{
						effect->at(cid).at(k).c_length+=d;
						effect->at(cid).at(k).c_phi+=phi;
						effect->at(cid).at(k).c_psi+=psi;
					}
				}
			}
		}
	}
	return;
}


void CTraj::gethbond2(vector<bbhbond_group> *hbond,vector< vector<ehbond> > *effect)
{
	int i,j,k;
	int base;
	int nid,cid;
	int n,h,c,o;
	double u[3];
	double x2,x3,x4,x5,y2,y3,y4,y5,z2,z3,z4,z5;
	double d,phi,psi;

	effect->resize(nres);
	for(i=0;i<(int)effect->size();i++)
	{
		effect->at(i).resize(nframe);
		for(k=0;k<nframe;k++)
		{
			effect->at(i).at(k).n_length=0;
			effect->at(i).at(k).c_length=0;
			effect->at(i).at(k).n_phi=0;
			effect->at(i).at(k).c_phi=0;
			effect->at(i).at(k).n_psi=0;
			effect->at(i).at(k).c_psi=0;
		}
	}

	for(i=0;i<(int)hbond->size();i++)
	{

		for(j=0;j<(int)hbond->size();j++)
		{
			k=j-i;
			if(k<=3 && k>=-3)
				continue;
			n=hbond->at(i).npos;
			h=hbond->at(i).hpos;
			if(h<0) 
				continue;
			c=hbond->at(j).cpos;
			o=hbond->at(j).opos;
			if(o<0)
				continue;
			n--;h--;c--;o--;


			for(k=0;k<nframe;k++)
			{
				//cout<<"i,j,k is "<<i<<" "<<j<<" "<<k<<endl;
				base=k*natom;
				u[0]=x[h+base]-x[o+base];
				u[1]=y[h+base]-y[o+base];
				u[2]=z[h+base]-z[o+base];
				d=veclength(u);
				x2=x[n+base];y2=y[n+base];z2=z[n+base];
				x3=x[h+base];y3=y[h+base];z3=z[h+base];
				x4=x[o+base];y4=y[o+base];z4=z[o+base];
				x5=x[c+base];y5=y[c+base];z5=z[c+base];
				phi=coor_to_angle(x2,y2,z2,x3,y3,z3,x4,y4,z4);
				psi=coor_to_angle(x3,y3,z3,x4,y4,z4,x5,y5,z5);
				if(d<3 && phi>0.5 && psi>0.5)
				{
					d=1/(d-1);
					phi*=phi;
					psi*=psi;
					nid=hbond->at(i).id-1;
					cid=hbond->at(j).id-1; /* C start from 0*/
					if(hbond->at(i).type==1)
					{
						effect->at(nid).at(k).n_length=d;
						effect->at(nid).at(k).n_phi=phi;
						effect->at(nid).at(k).n_psi=psi;
					}
					if(hbond->at(i).type==1)
					{
						effect->at(cid).at(k).c_length=d;
						effect->at(cid).at(k).c_phi=phi;
						effect->at(cid).at(k).c_psi=psi;
					}
				}
			}
		}
	}
	return;
}




void CTraj::dis_matrix(vector<int> *index,vector<double> *dis)
{
	int i,j,n1,n2;
	double d1,d2,d3,d;
	for(i=0;i<(int)index->size();i++)
	{
		n1=index->at(i)-1;
		for(j=0;j<(int)index->size();j++)
		{		
			n2=index->at(j)-1;
			d1=x.at(n1)-x.at(n2);
			d2=y.at(n1)-y.at(n2);
			d3=z.at(n1)-z.at(n2);
			d=sqrt(d1*d1+d2*d2+d3*d3);
			dis->push_back(d);
		}
	}
	return;
}


void CTraj::getangle(vector<struct dihe_group> *index, vector<double> * angle)
{
	int i,j;
	int x1,x2,x3,x4;
	int base;
	double t;


	for(i=0;i<nframe;i++)
	{
		base=i*natom-1;  // -1 because of C starts from 0 but not 1 !
		for(j=0;j<(int)index->size();j++)
		{
			x1=base+index->at(j).x1;
			x2=base+index->at(j).x2;
			x3=base+index->at(j).x3;
			x4=base+index->at(j).x4;
			t=coor_to_angle(x[x1],y[x1],z[x1],x[x2],y[x2],z[x2],x[x3],y[x3],z[x3]);
			angle->push_back(t);
			t=coor_to_angle(x[x2],y[x2],z[x2],x[x3],y[x3],z[x3],x[x4],y[x4],z[x4]);
			angle->push_back(t);
		}
	}
	return;
}



void CTraj::getdihe(vector<struct dihe_group> *index, vector<double> * dihe)
{	
	double st = omp_get_wtime();
	int i,j;
	int x1,x2,x3,x4;
	int base;
	double t;


	for(i=0;i<nframe;i++)
	{
		base=i*natom-1;  // -1 because of C starts from 0 but not 1 !
		for(j=0;j<(int)index->size();j++)
		{
			if(index->at(j).bgood==1)
			{
				x1=base+index->at(j).x1;
				x2=base+index->at(j).x2;
				x3=base+index->at(j).x3;
				x4=base+index->at(j).x4;
				t=coor_to_dihe(x[x1],y[x1],z[x1],x[x2],y[x2],z[x2],x[x3],y[x3],z[x3],x[x4],y[x4],z[x4]);
				dihe->push_back(t);
			}
			else
				dihe->push_back(-1000.0);
		}
	}
	cout << "traj::getdihe: " << omp_get_wtime()-st << " seconds" << endl;
	return;
}


void CTraj::getring(vector<struct ring_group> *index, vector<struct nh_group>* select, vector<struct double_five> *ring_effect)
{	
	double st = omp_get_wtime();
	int i,j,ii,jj,m;
	int base;
	int t[6];
	double p1[3],t1[3],t2[3],t3[3];
	double u[6][3];
	double sum[3];
	double ori[3];
	double e;
	struct double_five temp;
	
	for(i=0;i<5;i++)
		temp.x[i]=0;
	for(i=0;i<(int)select->size();i++)
		ring_effect->push_back(temp);


	for(i=0;i<nframe;i++)
	{
		base=i*natom;  
		for(j=0;j<(int)index->size();j++)
		{
			switch(index->at(j).x1)
			{
				case 1:
		        	case 2:
		       		case 5:
				m=6; break;
				case 4:
				case 3:
				m=5;break;
			}

			t[0]=index->at(j).x2;
			t[1]=index->at(j).x3;
			t[2]=index->at(j).x4;
			t[3]=index->at(j).x5;
			t[4]=index->at(j).x6;
			t[5]=index->at(j).x7;

           
            		for(ii=0;ii<m;ii++)
			{
				u[ii][0]=x[t[ii]+base-1];
                		u[ii][1]=y[t[ii]+base-1];
                		u[ii][2]=z[t[ii]+base-1];
            		}
           
            		for(jj=0;jj<3;jj++)
				sum[jj]=0; 
			for(ii=0;ii<m;ii++)
            		{
				for(jj=0;jj<3;jj++)
				{
					sum[jj]+=u[ii][jj];
				}
			}
			for(jj=0;jj<3;jj++)
				sum[jj]/=m; 
			for(ii=0;ii<m;ii++)
			{
				for(jj=0;jj<3;jj++)
					u[ii][jj]-=sum[jj];
			}

			ring(u,m,ori);

			for(jj=0;jj<3;jj++)
			{
				t1[jj]=u[0][jj]-u[1][jj];
				t2[jj]=u[2][jj]-u[1][jj];
			}
           		cross(t3,t1,t2);
            		if(dot(t3,ori)<0)
            		{
				for(jj=0;jj<3;jj++)
					ori[jj]=-ori[jj];
			}

			for(ii=0;ii<(int)select->size();ii++)
			{
				if(select->at(ii).hpos>=1)
				{
					e=0;
					p1[0]=x[base+select->at(ii).hpos-1]-sum[0];
					p1[1]=y[base+select->at(ii).hpos-1]-sum[1];
					p1[2]=z[base+select->at(ii).hpos-1]-sum[2];
					e+=effect(u,m,ori,p1); 
					e=e*10;
					ring_effect->at(ii).x[index->at(j).x1-1]+=e;
				}
			}
		}
	}

	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<5;jj++)
		{
			ring_effect->at(ii).x[jj]/=nframe;
		}
	}
	cout << "getring1: " << omp_get_wtime() - st << " seconds" << endl;
	return;
}


// New function for OpenACC
void CTraj::getring_acc(ring_group *index, int index_size, nh_group *select, int select_size, double_five *ring_effect_arr, int ring_effect_size)
{	
	double st = omp_get_wtime();
	int i,j,ii,jj,m;
	int base;

	//double_five *ring_effect_arr = ring_effect->data();
	//#pragma acc enter data copyin(ring_effect_arr[0:select_size])

	/*double *x_arr = x_arr;
	double *y_arr = y_arr;
	double *z_arr = z_arr;
	int x_size = x_size;
	int y_size = y_size;
	int z_size = z_size;*/

	#pragma acc data present(index[0:index_size],select[0:select_size], \
		x_arr[0:x_size],y_arr[0:y_size],z_arr[0:z_size],ring_effect_arr[0:select_size])
	{

	for(i=0;i<nframe;i++)
	{
		base=i*natom;  
		#pragma acc parallel loop independent gang private(i,j,ii,jj,m)
		for(j=0;j<index_size;j++)
		{
			int t_p[6];
			double t1_p[3],t2_p[3],t3_p[3];
			double u_p[6][3];
			double sum_p[3];
			double ori_p[3];

			if(index[j].x1 == 1 || index[j].x1 == 2 || index[j].x1 == 5){
				m=6;
			} else if(index[j].x1 == 4 || index[j].x1 == 3) {
				m=5;
			}

			t_p[0]=index[j].x2;
			t_p[1]=index[j].x3;
			t_p[2]=index[j].x4;
			t_p[3]=index[j].x5;
			t_p[4]=index[j].x6;
			t_p[5]=index[j].x7;

           		/*#pragma acc loop seq
            		for(ii=0;ii<m;ii++)
			{
				u_p[ii][0]=x_arr[t_p[ii]+base-1];
                		u_p[ii][1]=y_arr[t_p[ii]+base-1];
                		u_p[ii][2]=z_arr[t_p[ii]+base-1];
            		}*/
			u_p[0][0]=x_arr[t_p[0]+base-1]; u_p[0][1]=y_arr[t_p[0]+base-1]; u_p[0][2]=z_arr[t_p[0]+base-1];
			u_p[1][0]=x_arr[t_p[1]+base-1]; u_p[1][1]=y_arr[t_p[1]+base-1]; u_p[1][2]=z_arr[t_p[1]+base-1];
			u_p[2][0]=x_arr[t_p[2]+base-1]; u_p[2][1]=y_arr[t_p[2]+base-1]; u_p[2][2]=z_arr[t_p[2]+base-1];
			u_p[3][0]=x_arr[t_p[3]+base-1]; u_p[3][1]=y_arr[t_p[3]+base-1]; u_p[3][2]=z_arr[t_p[3]+base-1];
			u_p[4][0]=x_arr[t_p[4]+base-1]; u_p[4][1]=y_arr[t_p[4]+base-1]; u_p[4][2]=z_arr[t_p[4]+base-1];
			if(m >= 6) {
				u_p[5][0]=x_arr[t_p[5]+base-1]; u_p[5][1]=y_arr[t_p[5]+base-1]; u_p[5][2]=z_arr[t_p[5]+base-1];
			} else {
				u_p[5][0]=0; u_p[5][1]=0; u_p[5][2]=0;
			}
           
           		/*#pragma acc loop seq
            		for(jj=0;jj<3;jj++)
				sum_p[jj]=0;

           		#pragma acc loop seq
			for(ii=0;ii<m;ii++)
            		{
				sum_p[ii] = 
           			#pragma acc loop seq
				for(jj=0;jj<3;jj++)
				{
					sum_p[jj]+=u_p[ii][jj];
				}
			}*/

			sum_p[0] = (u_p[0][0] + u_p[1][0] + u_p[2][0] + u_p[3][0] + u_p[4][0] + u_p[5][0]) / m;
			sum_p[1] = (u_p[0][1] + u_p[1][1] + u_p[2][1] + u_p[3][1] + u_p[4][1] + u_p[5][1]) / m;
			sum_p[2] = (u_p[0][2] + u_p[1][2] + u_p[2][2] + u_p[3][2] + u_p[4][2] + u_p[5][2]) / m;

           		/*#pragma acc loop seq
			for(jj=0;jj<3;jj++)
				sum_p[jj]/=m; */

           		/*#pragma acc loop seq
			for(ii=0;ii<m;ii++)
			{
           			#pragma acc loop seq
				for(jj=0;jj<3;jj++)
					u_p[ii][jj]-=sum_p[jj];
			}*/

			u_p[0][0] -= sum_p[0]; u_p[0][1] -= sum_p[1]; u_p[0][2] -= sum_p[2];
			u_p[1][0] -= sum_p[0]; u_p[1][1] -= sum_p[1]; u_p[1][2] -= sum_p[2];
			u_p[2][0] -= sum_p[0]; u_p[2][1] -= sum_p[1]; u_p[2][2] -= sum_p[2];
			u_p[3][0] -= sum_p[0]; u_p[3][1] -= sum_p[1]; u_p[3][2] -= sum_p[2];
			u_p[4][0] -= sum_p[0]; u_p[4][1] -= sum_p[1]; u_p[4][2] -= sum_p[2];
			u_p[5][0] -= sum_p[0]; u_p[5][1] -= sum_p[1]; u_p[5][2] -= sum_p[2];
			

			my_ring(u_p,m,ori_p);

           		/*#pragma acc loop seq
			for(jj=0;jj<3;jj++)
			{
				t1_p[jj]=u_p[0][jj]-u_p[1][jj];
				t2_p[jj]=u_p[2][jj]-u_p[1][jj];
			}*/

			t1_p[0] = u_p[0][0] - u_p[1][0];
			t1_p[1] = u_p[0][1] - u_p[1][1];
			t1_p[2] = u_p[0][2] - u_p[1][2];
			t2_p[0] = u_p[2][0] - u_p[1][0];
			t2_p[1] = u_p[2][1] - u_p[1][1];
			t2_p[2] = u_p[2][2] - u_p[1][2];

           		my_cross(t3_p,t1_p,t2_p);
            		if(my_dot(t3_p,ori_p)<0)
            		{
	           		/*#pragma acc loop seq
				for(jj=0;jj<3;jj++)
					ori_p[jj]=-ori_p[jj];*/
				ori_p[0] = -ori_p[0];
				ori_p[1] = -ori_p[1];
				ori_p[2] = -ori_p[2];
			}

			#pragma acc loop vector independent
			for(ii=0;ii<select_size;ii++)
			{
				double e_p;
				double p1_p[3];
				if(select[ii].hpos>=1)
				{
					p1_p[0]=x_arr[base+select[ii].hpos-1]-sum_p[0];
					p1_p[1]=y_arr[base+select[ii].hpos-1]-sum_p[1];
					p1_p[2]=z_arr[base+select[ii].hpos-1]-sum_p[2];
					e_p=my_effect(u_p,m,ori_p,p1_p); 
					e_p*=10;
					#pragma acc atomic update
					ring_effect_arr[ii].x[index[j].x1-1]+=e_p;
				}
			}
		}
	}

	#pragma acc parallel loop independent
	for(ii=0;ii<select_size;ii++)
	{
		/*#pragma acc loop seq
		for(jj=0;jj<5;jj++)
		{
			ring_effect_arr[ii].x[jj]/=nframe;
		}*/
		ring_effect_arr[ii].x[0]/=nframe;
		ring_effect_arr[ii].x[1]/=nframe;
		ring_effect_arr[ii].x[2]/=nframe;
		ring_effect_arr[ii].x[3]/=nframe;
		ring_effect_arr[ii].x[4]/=nframe;
	}

	} // END DATA REGION

	cout << "getring(nh group): " << omp_get_wtime() - st << " seconds" << endl;
	return;
}


void CTraj::getring(vector<struct ring_group> *index, vector<struct nh_group>* select, vector< vector<struct double_five> > *ring_effect)
{	
	int i,j,ii,jj,m;
	int base;
	int t[6];
	double p1[3],t1[3],t2[3],t3[3];
	double u[6][3];
	double sum[3];
	double ori[3];
	double e;
	struct double_five temp;
	
	for(i=0;i<5;i++)
		temp.x[i]=0;
	ring_effect->resize(select->size());
	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<nframe;jj++)
		{
			ring_effect->at(ii).push_back(temp);
		}
	}


	for(i=0;i<nframe;i++)
	{
		base=i*natom;  
		for(j=0;j<(int)index->size();j++)
		{
			switch(index->at(j).x1)
			{
				case 1:
                case 2:
                case 5:
					m=6;
					break;
				case 4:
				case 3:
					m=5;
					break;
			}
			t[0]=index->at(j).x2;
			t[1]=index->at(j).x3;
			t[2]=index->at(j).x4;
			t[3]=index->at(j).x5;
			t[4]=index->at(j).x6;
			t[5]=index->at(j).x7;

           
            for(ii=0;ii<m;ii++)
			{
				u[ii][0]=x[t[ii]+base-1];
                u[ii][1]=y[t[ii]+base-1];
                u[ii][2]=z[t[ii]+base-1];
            }
           
            for(jj=0;jj<3;jj++)
				sum[jj]=0; 
			for(ii=0;ii<m;ii++)
            {
				for(jj=0;jj<3;jj++)
				{
					sum[jj]+=u[ii][jj];
				}
			}
			for(jj=0;jj<3;jj++)
				sum[jj]/=m; 
			for(ii=0;ii<m;ii++)
			{
				for(jj=0;jj<3;jj++)
					u[ii][jj]-=sum[jj];
			}

			ring(u,m,ori);

            for(jj=0;jj<3;jj++)
            {
                   t1[jj]=u[0][jj]-u[1][jj];
                   t2[jj]=u[2][jj]-u[1][jj];
            }
            cross(t3,t1,t2);
            if(dot(t3,ori)<0)
            {
				for(jj=0;jj<3;jj++)
					ori[jj]=-ori[jj];
			}

			for(ii=0;ii<(int)select->size();ii++)
			{
				if(select->at(ii).hpos>=1)
				{
					e=0;
					p1[0]=x[base+select->at(ii).hpos-1]-sum[0];
					p1[1]=y[base+select->at(ii).hpos-1]-sum[1];
					p1[2]=z[base+select->at(ii).hpos-1]-sum[2];
					e+=effect(u,m,ori,p1); 
					e=e*10;
					ring_effect->at(ii).at(i).x[index->at(j).x1-1]+=e;
				}
			}
		}
	}

	return;
}



// New function for OpenACC
void CTraj::getring_acc(ring_group *index, int index_size, proton *select, int select_size, double_five *ring_effect_arr, int ring_effect_size)
{	
	double st = omp_get_wtime();
	int i,j,ii,jj,m,k;
	int base;

	//double_five *ring_effect_arr = ring_effect->data();
	//#pragma acc enter data copyin(ring_effect_arr[0:select_size])

	/*double *x_arr = x_arr;
	double *y_arr = y_arr;
	double *z_arr = z_arr;
	int x_size = x_size;
	int y_size = y_size;
	int z_size = z_size;*/

	#pragma acc data present(index[0:index_size],select[0:select_size], \
		x_arr[0:x_size],y_arr[0:y_size],z_arr[0:z_size],ring_effect_arr[0:select_size])
	{

	for(i=0;i<nframe;i++)
	{
		base=i*natom;  
		#pragma acc parallel loop independent gang private(i,j,ii,jj,m)
		for(j=0;j<index_size;j++)
		{
			int t_p[6];
			double t1_p[3],t2_p[3],t3_p[3];
			double u_p[6][3];
			double sum_p[3];
			double ori_p[3];

			if(index[j].x1 == 1 || index[j].x1 == 2 || index[j].x1 == 5){
				m=6;
			} else if(index[j].x1 == 4 || index[j].x1 == 3) {
				m=5;
			}

			t_p[0]=index[j].x2;
			t_p[1]=index[j].x3;
			t_p[2]=index[j].x4;
			t_p[3]=index[j].x5;
			t_p[4]=index[j].x6;
			t_p[5]=index[j].x7;

           		/*#pragma acc loop seq
            		for(ii=0;ii<m;ii++)

			{

				u_p[ii][0]=x_arr[t_p[ii]+base-1];

                		u_p[ii][1]=y_arr[t_p[ii]+base-1];

                		u_p[ii][2]=z_arr[t_p[ii]+base-1];

            		}*/
			u_p[0][0]=x_arr[t_p[0]+base-1]; u_p[0][1]=y_arr[t_p[0]+base-1]; u_p[0][2]=z_arr[t_p[0]+base-1];
			u_p[1][0]=x_arr[t_p[1]+base-1]; u_p[1][1]=y_arr[t_p[1]+base-1]; u_p[1][2]=z_arr[t_p[1]+base-1];
			u_p[2][0]=x_arr[t_p[2]+base-1]; u_p[2][1]=y_arr[t_p[2]+base-1]; u_p[2][2]=z_arr[t_p[2]+base-1];
			u_p[3][0]=x_arr[t_p[3]+base-1]; u_p[3][1]=y_arr[t_p[3]+base-1]; u_p[3][2]=z_arr[t_p[3]+base-1];
			u_p[4][0]=x_arr[t_p[4]+base-1]; u_p[4][1]=y_arr[t_p[4]+base-1]; u_p[4][2]=z_arr[t_p[4]+base-1];
			if(m >= 6) {
				u_p[5][0]=x_arr[t_p[5]+base-1]; u_p[5][1]=y_arr[t_p[5]+base-1]; u_p[5][2]=z_arr[t_p[5]+base-1];
			} else {
				u_p[5][0]=0; u_p[5][1]=0; u_p[5][2]=0;
			}
           
           		/*#pragma acc loop seq
            		for(jj=0;jj<3;jj++)

				sum_p[jj]=0;

           		#pragma acc loop seq
			for(ii=0;ii<m;ii++)
            		{
				sum_p[ii] = 
           			#pragma acc loop seq
				for(jj=0;jj<3;jj++)
				{
					sum_p[jj]+=u_p[ii][jj];
				}
			}*/

			sum_p[0] = (u_p[0][0] + u_p[1][0] + u_p[2][0] + u_p[3][0] + u_p[4][0] + u_p[5][0]) / m;
			sum_p[1] = (u_p[0][1] + u_p[1][1] + u_p[2][1] + u_p[3][1] + u_p[4][1] + u_p[5][1]) / m;
			sum_p[2] = (u_p[0][2] + u_p[1][2] + u_p[2][2] + u_p[3][2] + u_p[4][2] + u_p[5][2]) / m;

           		/*#pragma acc loop seq
			for(jj=0;jj<3;jj++)

				sum_p[jj]/=m; */

           		/*#pragma acc loop seq
			for(ii=0;ii<m;ii++)

			{

           			#pragma acc loop seq

				for(jj=0;jj<3;jj++)

					u_p[ii][jj]-=sum_p[jj];

			}*/

			u_p[0][0] -= sum_p[0]; u_p[0][1] -= sum_p[1]; u_p[0][2] -= sum_p[2];
			u_p[1][0] -= sum_p[0]; u_p[1][1] -= sum_p[1]; u_p[1][2] -= sum_p[2];
			u_p[2][0] -= sum_p[0]; u_p[2][1] -= sum_p[1]; u_p[2][2] -= sum_p[2];
			u_p[3][0] -= sum_p[0]; u_p[3][1] -= sum_p[1]; u_p[3][2] -= sum_p[2];
			u_p[4][0] -= sum_p[0]; u_p[4][1] -= sum_p[1]; u_p[4][2] -= sum_p[2];
			u_p[5][0] -= sum_p[0]; u_p[5][1] -= sum_p[1]; u_p[5][2] -= sum_p[2];
			

			my_ring(u_p,m,ori_p);

           		/*#pragma acc loop seq
			for(jj=0;jj<3;jj++)

			{

				t1_p[jj]=u_p[0][jj]-u_p[1][jj];

				t2_p[jj]=u_p[2][jj]-u_p[1][jj];

			}*/

			t1_p[0] = u_p[0][0] - u_p[1][0];
			t1_p[1] = u_p[0][1] - u_p[1][1];
			t1_p[2] = u_p[0][2] - u_p[1][2];
			t2_p[0] = u_p[2][0] - u_p[1][0];
			t2_p[1] = u_p[2][1] - u_p[1][1];
			t2_p[2] = u_p[2][2] - u_p[1][2];

           		my_cross(t3_p,t1_p,t2_p);
            		if(my_dot(t3_p,ori_p)<0)
            		{
	           		/*#pragma acc loop seq
				for(jj=0;jj<3;jj++)

					ori_p[jj]=-ori_p[jj];*/
				ori_p[0] = -ori_p[0];
				ori_p[1] = -ori_p[1];
				ori_p[2] = -ori_p[2];
			}

			#pragma acc loop vector
			for(ii=0;ii<select_size;ii++)
			{
				double p1_p[3];
				double e_pp=0;
				#pragma acc loop seq reduction(+:e_pp)
				for(k=0;k<select[ii].nh;k++)
				{
					p1_p[0]=x_arr[base+select[ii].hpos[k]-1]-sum_p[0];
					p1_p[1]=y_arr[base+select[ii].hpos[k]-1]-sum_p[1];
					p1_p[2]=z_arr[base+select[ii].hpos[k]-1]-sum_p[2];
					e_pp+=my_effect(u_p,m,ori_p,p1_p); 
				}
				e_pp*=(10*3/select[ii].nh);
				#pragma acc atomic update
				ring_effect_arr[ii].x[index[j].x1-1]+=e_pp;
			}
		}
	}

	#pragma acc parallel loop independent
	for(ii=0;ii<select_size;ii++)
	{
		/*#pragma acc loop seq
		for(jj=0;jj<5;jj++)

		{

			ring_effect_arr[ii].x[jj]/=nframe;

		}*/
		ring_effect_arr[ii].x[0]/=nframe;
		ring_effect_arr[ii].x[1]/=nframe;
		ring_effect_arr[ii].x[2]/=nframe;
		ring_effect_arr[ii].x[3]/=nframe;
		ring_effect_arr[ii].x[4]/=nframe;
	}

	} // END DATA REGION

	cout << "getring(proton): " << omp_get_wtime() - st << " seconds" << endl;
	return;
}


void CTraj::getring(vector<struct ring_group> *index, vector<struct proton>* select, vector<struct double_five> *ring_effect)
{	
	int i,j,ii,jj,m,k;
	int base;
	int t[6];
	double p1[3],t1[3],t2[3],t3[3];
	double u[6][3];
	double sum[3];
	double ori[3];
	double e;
	struct double_five temp;
	
	for(i=0;i<5;i++)
		temp.x[i]=0;
	for(i=0;i<(int)select->size();i++)
		ring_effect->push_back(temp);


	for(i=0;i<nframe;i++)
	{
		base=i*natom;  
		for(j=0;j<(int)index->size();j++)
		{
			switch(index->at(j).x1)
			{
				case 1:
                case 2:
                case 5:
					m=6;
					break;
				case 4:
				case 3:
					m=5;
					break;
			}
			t[0]=index->at(j).x2;
			t[1]=index->at(j).x3;
			t[2]=index->at(j).x4;
			t[3]=index->at(j).x5;
			t[4]=index->at(j).x6;
			t[5]=index->at(j).x7;

           
            for(ii=0;ii<m;ii++)
			{
				u[ii][0]=x[t[ii]+base-1];
                u[ii][1]=y[t[ii]+base-1];
                u[ii][2]=z[t[ii]+base-1];
            }
           
            for(jj=0;jj<3;jj++)
				sum[jj]=0; 
			for(ii=0;ii<m;ii++)
            {
				for(jj=0;jj<3;jj++)
				{
					sum[jj]+=u[ii][jj];
				}
			}
			for(jj=0;jj<3;jj++)
				sum[jj]/=m; 
			for(ii=0;ii<m;ii++)
			{
				for(jj=0;jj<3;jj++)
					u[ii][jj]-=sum[jj];
			}

			ring(u,m,ori);

            for(jj=0;jj<3;jj++)
            {
                   t1[jj]=u[0][jj]-u[1][jj];
                   t2[jj]=u[2][jj]-u[1][jj];
            }
            cross(t3,t1,t2);
            if(dot(t3,ori)<0)
            {
				for(jj=0;jj<3;jj++)
					ori[jj]=-ori[jj];
			}

			for(ii=0;ii<(int)select->size();ii++)
			{
				e=0;	
				for(k=0;k<(int)select->at(ii).nh;k++)
				{
					p1[0]=x[base+select->at(ii).hpos[k]-1]-sum[0];
					p1[1]=y[base+select->at(ii).hpos[k]-1]-sum[1];
					p1[2]=z[base+select->at(ii).hpos[k]-1]-sum[2];
					e+=effect(u,m,ori,p1); 
				}
				e=e*10*3/select->at(ii).nh;
				ring_effect->at(ii).x[index->at(j).x1-1]+=e;
			}
		}
	}

	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<5;jj++)
		{
			ring_effect->at(ii).x[jj]/=nframe;
		}
	}
	return;
}





void CTraj::getring(vector<struct ring_group> *index, vector<struct proton>* select, vector< vector<struct double_five> > *ring_effect)
{	
	int i,j,ii,jj,m,k;
	int base;
	int t[6];
	double p1[3],t1[3],t2[3],t3[3];
	double u[6][3];
	double sum[3];
	double ori[3];
	double e;
	struct double_five temp;
	
	for(i=0;i<5;i++)
		temp.x[i]=0;
	ring_effect->resize(select->size());
	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<nframe;jj++)
		{
			ring_effect->at(ii).push_back(temp);
		}
	}


	for(i=0;i<nframe;i++)
	{
		base=i*natom;  
		for(j=0;j<(int)index->size();j++)
		{
			switch(index->at(j).x1)
			{
				case 1:
                case 2:
                case 5:
					m=6;
					break;
				case 4:
				case 3:
					m=5;
					break;
			}
			t[0]=index->at(j).x2;
			t[1]=index->at(j).x3;
			t[2]=index->at(j).x4;
			t[3]=index->at(j).x5;
			t[4]=index->at(j).x6;
			t[5]=index->at(j).x7;

           
            for(ii=0;ii<m;ii++)
			{
				u[ii][0]=x[t[ii]+base-1];
                u[ii][1]=y[t[ii]+base-1];
                u[ii][2]=z[t[ii]+base-1];
            }
           
            for(jj=0;jj<3;jj++)
				sum[jj]=0; 
			for(ii=0;ii<m;ii++)
            {
				for(jj=0;jj<3;jj++)
				{
					sum[jj]+=u[ii][jj];
				}
			}
			for(jj=0;jj<3;jj++)
				sum[jj]/=m; 
			for(ii=0;ii<m;ii++)
			{
				for(jj=0;jj<3;jj++)
					u[ii][jj]-=sum[jj];
			}

			ring(u,m,ori);

            for(jj=0;jj<3;jj++)
            {
                   t1[jj]=u[0][jj]-u[1][jj];
                   t2[jj]=u[2][jj]-u[1][jj];
            }
            cross(t3,t1,t2);
            if(dot(t3,ori)<0)
            {
				for(jj=0;jj<3;jj++)
					ori[jj]=-ori[jj];
			}

			for(ii=0;ii<(int)select->size();ii++)
			{
				e=0;	
				for(k=0;k<(int)select->at(ii).nh;k++)
				{
					p1[0]=x[base+select->at(ii).hpos[k]-1]-sum[0];
					p1[1]=y[base+select->at(ii).hpos[k]-1]-sum[1];
					p1[2]=z[base+select->at(ii).hpos[k]-1]-sum[2];
					e+=effect(u,m,ori,p1); 
				}
				e=e*10*3/select->at(ii).nh;
				ring_effect->at(ii).at(i).x[index->at(j).x1-1]+=e;
			}
		}
	}
	return;
}


void CTraj::getring_bb(vector<struct ring_group> *index, vector<struct bb_group>* select, vector<struct double_five> *ring_effect, enum bb_carbon c)
{	
	int i,j,ii,jj,m;
	int base;
	int t[6];
	double p1[3],t1[3],t2[3],t3[3];
	double u[6][3];
	double sum[3];
	double ori[3];
	double e;
	struct double_five temp;
	
	for(i=0;i<5;i++)
		temp.x[i]=0;
	for(i=0;i<(int)select->size();i++)
		ring_effect->push_back(temp);


	for(i=0;i<nframe;i++)
	{
		base=i*natom;  
		for(j=0;j<(int)index->size();j++)
		{
			switch(index->at(j).x1)
			{
				case 1:
                case 2:
                case 5:
					m=6;
					break;
				case 4:
				case 3:
					m=5;
					break;
			}
			t[0]=index->at(j).x2;
			t[1]=index->at(j).x3;
			t[2]=index->at(j).x4;
			t[3]=index->at(j).x5;
			t[4]=index->at(j).x6;
			t[5]=index->at(j).x7;

           
            for(ii=0;ii<m;ii++)
			{
				u[ii][0]=x[t[ii]+base-1];
                u[ii][1]=y[t[ii]+base-1];
                u[ii][2]=z[t[ii]+base-1];
            }
           
            for(jj=0;jj<3;jj++)
				sum[jj]=0; 
			for(ii=0;ii<m;ii++)
            {
				for(jj=0;jj<3;jj++)
				{
					sum[jj]+=u[ii][jj];
				}
			}
			for(jj=0;jj<3;jj++)
				sum[jj]/=m; 
			for(ii=0;ii<m;ii++)
			{
				for(jj=0;jj<3;jj++)
					u[ii][jj]-=sum[jj];
			}

			ring(u,m,ori);

            for(jj=0;jj<3;jj++)
            {
                   t1[jj]=u[0][jj]-u[1][jj];
                   t2[jj]=u[2][jj]-u[1][jj];
            }
            cross(t3,t1,t2);
            if(dot(t3,ori)<0)
            {
				for(jj=0;jj<3;jj++)
					ori[jj]=-ori[jj];
			}

			for(ii=0;ii<(int)select->size();ii++)
			{
				e=0;

				if((select->at(ii).capos>=1 && c==bb_ca) || (select->at(ii).cbpos>=1 && c==bb_cb)||(select->at(ii).copos>=1 && c==bb_co))
				switch (c)
				{
				case bb_ca:
					p1[0]=x[base+select->at(ii).capos-1]-sum[0];
					p1[1]=y[base+select->at(ii).capos-1]-sum[1];
					p1[2]=z[base+select->at(ii).capos-1]-sum[2];
					break;
				case bb_cb:
					p1[0]=x[base+select->at(ii).cbpos-1]-sum[0];
					p1[1]=y[base+select->at(ii).cbpos-1]-sum[1];
					p1[2]=z[base+select->at(ii).cbpos-1]-sum[2];
					break;
				case bb_co:
					p1[0]=x[base+select->at(ii).copos-1]-sum[0];
					p1[1]=y[base+select->at(ii).copos-1]-sum[1];
					p1[2]=z[base+select->at(ii).copos-1]-sum[2];
					break;
				}
				e+=effect(u,m,ori,p1); 
				e=e*10;
				ring_effect->at(ii).x[index->at(j).x1-1]+=e;
			}
		}
	}

	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<5;jj++)
		{
			ring_effect->at(ii).x[jj]/=nframe;
		}
	}
	return;
}



void CTraj::getani(vector<struct ani_group> *index, vector<struct methyl_group>* select, vector<struct double_four> *ani_effect, enum methyl c)
{
	int i,j,ii,jj,k;
	int i1,i2,i3;
	int base;
	double center[3];
	double v1[3];
	double v2[3];
	double ori[3];
	double cosa;
	double length;
	double e;


	struct double_four temp;
	
	for(i=0;i<4;i++)
		temp.x[i]=0;
	for(i=0;i<(int)select->size();i++)
		ani_effect->push_back(temp);
	
	for(i=0;i<nframe;i++)
	{
		base=i*natom;
		for(j=0;j<(int)index->size();j++)
		{
			i1=index->at(j).pos[0]+base-1;
			i2=index->at(j).pos[1]+base-1;
			i3=index->at(j).pos[2]+base-1;
			center[0]=(x[i1]+x[i2]+x[i3])/3;
			center[1]=(y[i1]+y[i2]+y[i3])/3;
			center[2]=(z[i1]+z[i2]+z[i3])/3;

			v1[0]=x[i1]-x[i2];
			v1[1]=y[i1]-y[i2];
			v1[2]=z[i1]-z[i2];

			v2[0]=x[i3]-x[i2];
			v2[1]=y[i3]-y[i2];
			v2[2]=z[i3]-z[i2];

			cross(ori,v1,v2);

			for(jj=0;jj<(int)select->size();jj++)
			{
				e=0;				
				switch(c)
				{
				case hydrogen:
					for(k=0;k<3;k++)
					{
						i1=base+select->at(jj).cpos-1+k+1;
						v1[0]=center[0]-x[i1];
						v1[1]=center[1]-y[i1];
						v1[2]=center[2]-z[i1];
						length=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
						cosa=dot(v1,ori);
						cosa/=sqrt(ori[0]*ori[0]+ori[1]*ori[1]+ori[2]*ori[2]);
						cosa/=sqrt(length);					 
						e+=(1-3*cosa*cosa)/(length*sqrt(length));
					}
					ani_effect->at(jj).x[index->at(j).type-1]+=e/3.0*1000;
					break;

				case carbon:
					i1=base+select->at(jj).cpos-1;
					v1[0]=center[0]-x[i1];
					v1[1]=center[1]-y[i1];
					v1[2]=center[2]-z[i1];
					length=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
					cosa=dot(v1,ori);
					cosa/=sqrt(ori[0]*ori[0]+ori[1]*ori[1]+ori[2]*ori[2]);
					cosa/=sqrt(length);					 
					e+=(1-3*cosa*cosa)/(length*sqrt(length));
					ani_effect->at(jj).x[index->at(j).type-1]+=e*1000;
					break;
				}
			}
		}
	}

	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<4;jj++)
		{
			ani_effect->at(ii).x[jj]/=nframe;
		}
	}
	return;
}


// New function for OpenACC
void CTraj::getani_acc(ani_group *index, int index_size, proton *select, int select_size, double_four *ani_effect_arr)
{
	if(index_size <= 0 || select_size <= 0)
		return;
	double st = omp_get_wtime();
	int i,j,ii,jj,k;
	int i1,i2,i3;
	int base;
	double cosa;
	double length;
	double e;

	//double_four *ani_effect_arr = ani_effect->data();
	//#pragma acc enter data copyin(ani_effect_arr[0:select_size])

	/*double *x_arr = x_arr;
	double *y_arr = y_arr;
	double *z_arr = z_arr;
	int x_size = x_size;
	int y_size = y_size;
	int z_size = z_size;*/

	#pragma acc data present(index[0:index_size],select[0:select_size], \
	x_arr[0:x_size],y_arr[0:y_size],z_arr[0:z_size], \
	ani_effect_arr[0:select_size])
	{

	for(i=0;i<nframe;i++)
	{

		base=i*natom;
		#pragma acc parallel
		{
		#pragma acc loop independent gang private(i1,i2,i3,e,cosa,length,jj,k)
		for(j=0;j<index_size;j++)
		{
			double center_p[3];
			double v1_p[3];
			double v2_p[3];
			double ori_p[3];
			i1=index[j].pos[0]+base-1;
			i2=index[j].pos[1]+base-1;
			i3=index[j].pos[2]+base-1;

			center_p[0]=(x_arr[i1]+x_arr[i2]+x_arr[i3])/3;
			center_p[1]=(y_arr[i1]+y_arr[i2]+y_arr[i3])/3; 
			center_p[2]=(z_arr[i1]+z_arr[i2]+z_arr[i3])/3;

			v1_p[0]=x_arr[i1]-x_arr[i2];
			v1_p[1]=y_arr[i1]-y_arr[i2];
			v1_p[2]=z_arr[i1]-z_arr[i2];

			v2_p[0]=x_arr[i3]-x_arr[i2];
			v2_p[1]=y_arr[i3]-y_arr[i2];
			v2_p[2]=z_arr[i3]-z_arr[i2];

			my_cross(ori_p,v1_p,v2_p);

			#pragma acc loop vector
			for(jj=0;jj<select_size;jj++) 
			{
				double e_pp = 0;
				int i1_pp;
				double length_pp;
				double cosa_pp;
				double v1_pp[3];


				// This loop is actually only either size 1 or 2
				#pragma acc loop seq reduction(+:e)				
				for(k=0;k<select[jj].nh;k++)
				{	
					i1_pp=base+select[jj].hpos[k]-1;

					v1_pp[0]=center_p[0]-x_arr[i1_pp];
					v1_pp[1]=center_p[1]-y_arr[i1_pp];
					v1_pp[2]=center_p[2]-z_arr[i1_pp];

					length_pp=v1_pp[0]*v1_pp[0]+v1_pp[1]*v1_pp[1]+v1_pp[2]*v1_pp[2];

					cosa_pp=my_dot(v1_pp,ori_p);

					cosa_pp/=sqrt(ori_p[0]*ori_p[0]+ori_p[1]*ori_p[1]+ori_p[2]*ori_p[2]);
					cosa_pp/=sqrt(length_pp);					 
					e_pp+=(1-3*cosa_pp*cosa_pp)/(length_pp*sqrt(length_pp));
				} // for(k=0;k<select[jj].nh;k++)

				#pragma acc atomic update
				ani_effect_arr[jj].x[index[j].type-1] += e_pp/select[jj].nh*1000;
			} // for(jj=0;jj<select_size;jj++)
		} // for(j=0;j<index_size;j++)

		} // END PARALLEL REGION

	} // for(i=0;i<nframe;i++)

	#pragma acc parallel loop independent
	for(j=0; j<select_size; j++)
	{
		ani_effect_arr[j].x[0] /= nframe;
		ani_effect_arr[j].x[1] /= nframe;
		ani_effect_arr[j].x[2] /= nframe;
		ani_effect_arr[j].x[3] /= nframe;
	}
	} // END DATA REGION


	cout << "getani(proton): " << omp_get_wtime() - st << " seconds" << endl;
	return;
}



// New function for OpenACC
void CTraj::getani_acc(ani_group *index, int index_size, proton *select, int select_size, double_four *ani_effect_arr, int ani_effect_size)
{
	if(index_size <= 0 || select_size <= 0)
		return;
	double st = omp_get_wtime();
	int i,j,ii,jj,k;
	int i1,i2,i3;
	int base;
	double cosa;
	double length;
	double e;

	//double_four *ani_effect_arr = ani_effect->data();
	//#pragma acc enter data copyin(ani_effect_arr[0:select_size])

	/*double *x_arr = x_arr;
	double *y_arr = y_arr;
	double *z_arr = z_arr;
	int x_size = x_size;
	int y_size = y_size;
	int z_size = z_size;*/

	#pragma acc data present(index[0:index_size],select[0:select_size], \
	x_arr[0:x_size],y_arr[0:y_size],z_arr[0:z_size], \
	ani_effect_arr[0:select_size])
	{

	for(i=0;i<nframe;i++)
	{

		base=i*natom;
		#pragma acc parallel
		{
		#pragma acc loop independent gang private(i1,i2,i3,e,cosa,length,jj,k)
		for(j=0;j<index_size;j++)
		{
			double center_p[3];
			double v1_p[3];
			double v2_p[3];
			double ori_p[3];
			i1=index[j].pos[0]+base-1;
			i2=index[j].pos[1]+base-1;
			i3=index[j].pos[2]+base-1;

			center_p[0]=(x_arr[i1]+x_arr[i2]+x_arr[i3])/3;
			center_p[1]=(y_arr[i1]+y_arr[i2]+y_arr[i3])/3; 
			center_p[2]=(z_arr[i1]+z_arr[i2]+z_arr[i3])/3;

			v1_p[0]=x_arr[i1]-x_arr[i2];
			v1_p[1]=y_arr[i1]-y_arr[i2];
			v1_p[2]=z_arr[i1]-z_arr[i2];

			v2_p[0]=x_arr[i3]-x_arr[i2];
			v2_p[1]=y_arr[i3]-y_arr[i2];
			v2_p[2]=z_arr[i3]-z_arr[i2];

			my_cross(ori_p,v1_p,v2_p);

			#pragma acc loop vector
			for(jj=0;jj<select_size;jj++) 
			{
				double e_pp = 0;
				int i1_pp;
				double length_pp;
				double cosa_pp;
				double v1_pp[3];


				// This loop is actually only either size 1 or 2
				#pragma acc loop seq reduction(+:e)				
				for(k=0;k<select[jj].nh;k++)
				{	
					i1_pp=base+select[jj].hpos[k]-1;

					v1_pp[0]=center_p[0]-x_arr[i1_pp];
					v1_pp[1]=center_p[1]-y_arr[i1_pp];
					v1_pp[2]=center_p[2]-z_arr[i1_pp];

					length_pp=v1_pp[0]*v1_pp[0]+v1_pp[1]*v1_pp[1]+v1_pp[2]*v1_pp[2];

					cosa_pp=my_dot(v1_pp,ori_p);

					cosa_pp/=sqrt(ori_p[0]*ori_p[0]+ori_p[1]*ori_p[1]+ori_p[2]*ori_p[2]);
					cosa_pp/=sqrt(length_pp);					 
					e_pp+=(1-3*cosa_pp*cosa_pp)/(length_pp*sqrt(length_pp));
				} // for(k=0;k<select[jj].nh;k++)

				#pragma acc atomic update
				ani_effect_arr[jj].x[index[j].type-1] += e_pp/select[jj].nh*1000;
			} // for(jj=0;jj<select_size;jj++)
		} // for(j=0;j<index_size;j++)

		} // END PARALLEL REGION

	} // for(i=0;i<nframe;i++)

	#pragma acc parallel loop independent
	for(j=0; j<select_size; j++)
	{
		ani_effect_arr[j].x[0] /= nframe;
		ani_effect_arr[j].x[1] /= nframe;
		ani_effect_arr[j].x[2] /= nframe;
		ani_effect_arr[j].x[3] /= nframe;
	}
	} // END DATA REGION


	cout << "getani(proton): " << omp_get_wtime() - st << " seconds" << endl;
	return;
}

void CTraj::getani(vector<struct ani_group> *index, vector<struct proton>* select, vector<struct double_four> *ani_effect)
{
	int i,j,ii,jj,k;
	int i1,i2,i3;
	int base;
	double center[3];
	double v1[3];
	double v2[3];
	double ori[3];
	double cosa;
	double length;
	double e;


	//struct double_four temp;
	
	//for(i=0;i<4;i++)
	//	temp.x[i]=0;
	//for(i=0;i<(int)select->size();i++)
	//	ani_effect->push_back(temp);

	double_four *ani_effect_arr = ani_effect->data();
	
	for(i=0;i<nframe;i++)
	{

		base=i*natom;
		for(j=0;j<(int)index->size();j++)
		{
			double ori_p[3];
			double center_p[3];
			double v1_p[3];
			double v2_p[3];
			i1=index->at(j).pos[0]+base-1;
			i2=index->at(j).pos[1]+base-1;
			i3=index->at(j).pos[2]+base-1;
			center_p[0]=(x[i1]+x[i2]+x[i3])/3;
			center_p[1]=(y[i1]+y[i2]+y[i3])/3;
			center_p[2]=(z[i1]+z[i2]+z[i3])/3;

			v1_p[0]=x[i1]-x[i2];
			v1_p[1]=y[i1]-y[i2];
			v1_p[2]=z[i1]-z[i2];

			v2_p[0]=x[i3]-x[i2];
			v2_p[1]=y[i3]-y[i2];
			v2_p[2]=z[i3]-z[i2];

			cross(ori_p,v1_p,v2_p);

			for(jj=0;jj<(int)select->size();jj++)
			{
				//cout<<jj<<endl;
				double e_p =0;
				double v1_pp[3];
				double length_p;
				double cosa_p;
				double i1_p;
				//e=0;				
				for(k=0;k<(int)select->at(jj).nh;k++)
				{
					i1_p=base+select->at(jj).hpos[k]-1;
					v1_p[0]=center_p[0]-x[i1_p];
					v1_p[1]=center_p[1]-y[i1_p];
					v1_p[2]=center_p[2]-z[i1_p];
					length_p=v1_p[0]*v1_p[0]+v1_p[1]*v1_p[1]+v1_p[2]*v1_p[2];
					cosa_p=dot(v1_p,ori_p);
					cosa_p/=sqrt(ori_p[0]*ori_p[0]+ori_p[1]*ori_p[1]+ori_p[2]*ori_p[2]);
					cosa_p/=sqrt(length_p);					 
					e_p+=(1-3*cosa_p*cosa_p)/(length_p*sqrt(length_p));
				}
				ani_effect_arr[jj].x[index->at(j).type-1]+=e_p/select->at(jj).nh*1000;
			}
		}
	}

	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<4;jj++)
		{
			ani_effect_arr[ii].x[jj]/=nframe;
		}
	}
	return;
}



void CTraj::getani(vector<struct ani_group> *index, vector<struct proton>* select, vector< vector<struct double_four>  > *ani_effect)
{
	int i,j,ii,jj,k;
	int i1,i2,i3;
	int base;
	double center[3];
	double v1[3];
	double v2[3];
	double ori[3];
	double cosa;
	double length;
	double e;


	struct double_four temp;
	
	for(i=0;i<4;i++)
		temp.x[i]=0;
	ani_effect->resize(select->size());
	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<nframe;jj++)
		{
			ani_effect->at(ii).push_back(temp);
		}
	}

	
	for(i=0;i<nframe;i++)
	{
		base=i*natom;
		for(j=0;j<(int)index->size();j++)
		{
			i1=index->at(j).pos[0]+base-1;
			i2=index->at(j).pos[1]+base-1;
			i3=index->at(j).pos[2]+base-1;
			center[0]=(x[i1]+x[i2]+x[i3])/3;
			center[1]=(y[i1]+y[i2]+y[i3])/3;
			center[2]=(z[i1]+z[i2]+z[i3])/3;

			v1[0]=x[i1]-x[i2];
			v1[1]=y[i1]-y[i2];
			v1[2]=z[i1]-z[i2];

			v2[0]=x[i3]-x[i2];
			v2[1]=y[i3]-y[i2];
			v2[2]=z[i3]-z[i2];

			cross(ori,v1,v2);

			for(jj=0;jj<(int)select->size();jj++)
			{
				e=0;				
				for(k=0;k<(int)select->at(jj).nh;k++)
				{
					i1=base+select->at(jj).hpos[k]-1;
					v1[0]=center[0]-x[i1];
					v1[1]=center[1]-y[i1];
					v1[2]=center[2]-z[i1];
					length=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
					cosa=dot(v1,ori);
					cosa/=sqrt(ori[0]*ori[0]+ori[1]*ori[1]+ori[2]*ori[2]);
					cosa/=sqrt(length);					 
					e+=(1-3*cosa*cosa)/(length*sqrt(length));
				}
				ani_effect->at(jj).at(i).x[index->at(j).type-1]+=e/select->at(jj).nh*1000;
			}
		}
	}


	return;
}


// New function for OpenACC
void CTraj::getani_acc(ani_group *index, int index_size, nh_group *select, int select_size, double_four *ani_effect_arr, int ani_effect_size)
{
	if(index_size <= 0 || select_size <= 0)
		return;
	double st = omp_get_wtime();
	int i,j,ii,jj;
	int i1,i2,i3;
	int base;
	double cosa;
	double length;
	double e;

	//double_four *ani_effect_arr = ani_effect->data();
	//#pragma acc enter data copyin(ani_effect_arr[0:select_size])

	//double *x_arr = x_arr;
	//double *y_arr = y_arr;
	//double *z_arr = z_arr;
	//int x_size = x_size;
	//int y_size = y_size;
	//int z_size = z_size;

	#pragma acc data present(index[0:index_size],select[0:select_size],x_arr[0:x_size],y_arr[0:y_size],z_arr[0:z_size],ani_effect_arr[0:select_size])
	{

	for(i=0;i<nframe;i++)
	{
		base=i*natom;
		#pragma acc parallel
		{

		#pragma acc loop independent gang private(i1,i2,i3,e,cosa,length,jj)
		for(j=0;j<index_size;j++)
		{
			double center_p[3];
			double v1_p[3];
			double v2_p[3];
			double ori_p[3];
			i1=index[j].pos[0]+base-1;
			i2=index[j].pos[1]+base-1;
			i3=index[j].pos[2]+base-1;
			center_p[0]=(x_arr[i1]+x_arr[i2]+x_arr[i3])/3;
			center_p[1]=(y_arr[i1]+y_arr[i2]+y_arr[i3])/3;
			center_p[2]=(z_arr[i1]+z_arr[i2]+z_arr[i3])/3;

			v1_p[0]=x_arr[i1]-x_arr[i2];
			v1_p[1]=y_arr[i1]-y_arr[i2];
			v1_p[2]=z_arr[i1]-z_arr[i2];

			v2_p[0]=x_arr[i3]-x_arr[i2];
			v2_p[1]=y_arr[i3]-y_arr[i2];
			v2_p[2]=z_arr[i3]-z_arr[i2];

			my_cross(ori_p,v1_p,v2_p);

			#pragma acc loop independent vector
			for(jj=0;jj<select_size;jj++)
			{
				double e_pp;
				int i1_pp;
				double v1_pp[3];
				double length_pp;
				double cosa_pp;
				if(select[jj].hpos>=1)
				{			
					i1_pp=base+select[jj].hpos-1;
					v1_pp[0]=center_p[0]-x_arr[i1_pp];
					v1_pp[1]=center_p[1]-y_arr[i1_pp];
					v1_pp[2]=center_p[2]-z_arr[i1_pp];
					length_pp=v1_pp[0]*v1_pp[0]+v1_pp[1]*v1_pp[1]+v1_pp[2]*v1_pp[2];
					cosa_pp=my_dot(v1_pp,ori_p);
					cosa_pp/=sqrt(ori_p[0]*ori_p[0]+ori_p[1]*ori_p[1]+ori_p[2]*ori_p[2]);
					cosa_pp/=sqrt(length_pp);					 
					e_pp=(1-3*cosa_pp*cosa_pp)/(length_pp*sqrt(length_pp));
					#pragma acc atomic update
					ani_effect_arr[jj].x[index[j].type-1]+=e_pp*1000;
				}
			}
		}
		} // end parallel region
	} // end frame loop

	#pragma acc parallel loop independent
	for(j=0; j<select_size; j++)
	{

		ani_effect_arr[j].x[0] /= nframe;
		ani_effect_arr[j].x[1] /= nframe;
		ani_effect_arr[j].x[2] /= nframe;
		ani_effect_arr[j].x[3] /= nframe;
	}

	} // END DATA REGION

	cout << "getani(nh group): " << omp_get_wtime() - st << " seconds" << endl;
	return;
}



void CTraj::getani(vector<struct ani_group> *index, vector<struct nh_group>* select, vector< vector<struct double_four>  > *ani_effect)
{
	int i,j,ii,jj;
	int i1,i2,i3;
	int base;
	double center[3];
	double v1[3];
	double v2[3];
	double ori[3];
	double cosa;
	double length;
	double e;


	struct double_four temp;
	
	for(i=0;i<4;i++)
		temp.x[i]=0;
	ani_effect->resize(select->size());
	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<nframe;jj++)
		{
			ani_effect->at(ii).push_back(temp);
		}
	}
	
	for(i=0;i<nframe;i++)
	{
		base=i*natom;
		for(j=0;j<(int)index->size();j++)
		{
			i1=index->at(j).pos[0]+base-1;
			i2=index->at(j).pos[1]+base-1;
			i3=index->at(j).pos[2]+base-1;
			center[0]=(x[i1]+x[i2]+x[i3])/3;
			center[1]=(y[i1]+y[i2]+y[i3])/3;
			center[2]=(z[i1]+z[i2]+z[i3])/3;

			v1[0]=x[i1]-x[i2];
			v1[1]=y[i1]-y[i2];
			v1[2]=z[i1]-z[i2];

			v2[0]=x[i3]-x[i2];
			v2[1]=y[i3]-y[i2];
			v2[2]=z[i3]-z[i2];

			cross(ori,v1,v2);

			for(jj=0;jj<(int)select->size();jj++)
			{
				if(select->at(jj).hpos>=1)
				{
					e=0;				
					i1=base+select->at(jj).hpos-1;
					v1[0]=center[0]-x[i1];
					v1[1]=center[1]-y[i1];
					v1[2]=center[2]-z[i1];
					length=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
					cosa=dot(v1,ori);
					cosa/=sqrt(ori[0]*ori[0]+ori[1]*ori[1]+ori[2]*ori[2]);
					cosa/=sqrt(length);					 
					e+=(1-3*cosa*cosa)/(length*sqrt(length));
					ani_effect->at(jj).at(i).x[index->at(j).type-1]+=e*1000;
				}
			}
		}
	}

	return;
}


void CTraj::getani(vector<struct ani_group> *index, vector<struct nh_group>* select, vector<struct double_four> *ani_effect)
{
	int i,j,ii,jj;
	int i1,i2,i3;
	int base;
	double center[3];
	double v1[3];
	double v2[3];
	double ori[3];
	double cosa;
	double length;
	double e;


	struct double_four temp;
	
	for(i=0;i<4;i++)
		temp.x[i]=0;
	for(i=0;i<(int)select->size();i++)
		ani_effect->push_back(temp);
	
	for(i=0;i<nframe;i++)
	{
		base=i*natom;
		for(j=0;j<(int)index->size();j++)
		{
			i1=index->at(j).pos[0]+base-1;
			i2=index->at(j).pos[1]+base-1;
			i3=index->at(j).pos[2]+base-1;
			center[0]=(x[i1]+x[i2]+x[i3])/3;
			center[1]=(y[i1]+y[i2]+y[i3])/3;
			center[2]=(z[i1]+z[i2]+z[i3])/3;

			v1[0]=x[i1]-x[i2];
			v1[1]=y[i1]-y[i2];
			v1[2]=z[i1]-z[i2];

			v2[0]=x[i3]-x[i2];
			v2[1]=y[i3]-y[i2];
			v2[2]=z[i3]-z[i2];

			cross(ori,v1,v2);

			for(jj=0;jj<(int)select->size();jj++)
			{
				if(select->at(jj).hpos>=1)
				{
					e=0;				
					i1=base+select->at(jj).hpos-1;
					v1[0]=center[0]-x[i1];
					v1[1]=center[1]-y[i1];
					v1[2]=center[2]-z[i1];
					length=v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
					cosa=dot(v1,ori);
					cosa/=sqrt(ori[0]*ori[0]+ori[1]*ori[1]+ori[2]*ori[2]);
					cosa/=sqrt(length);					 
					e+=(1-3*cosa*cosa)/(length*sqrt(length));
					ani_effect->at(jj).x[index->at(j).type-1]+=e*1000;
				}
			}
		}
	}

	for(ii=0;ii<(int)select->size();ii++)
	{
		for(jj=0;jj<4;jj++)
		{
			ani_effect->at(ii).x[jj]/=nframe;
		}
	}
	return;
}



double CTraj::noedistance_frame(vector<int> *att1, vector<int> *att2,int j)
{
	double x0,y0,z0;
	double r2,r6,sum,sumsum;
	int i1,j1;
	int n1,n2;
	int jump;

	sumsum=0;
	for(i1=0;i1<(int)att1->size();i1++)
	for(j1=0;j1<(int)att2->size();j1++)
	{
		n1=att1->at(i1)-1;
		n2=att2->at(j1)-1;
		
		jump=j*natom;
		x0=x.at(jump+n1)-x.at(jump+n2);
		y0=y.at(jump+n1)-y.at(jump+n2);
		z0=z.at(jump+n1)-z.at(jump+n2);
		r2=x0*x0+y0*y0+z0*z0;
		r6=r2*r2*r2;
		sum=1.0/r6;
		sumsum+=sum;
	}
	sumsum/=i1;
	sumsum/=j1;
	r2=pow(sumsum,-1/6.0);
	return r2;
}

double CTraj::noedistance(vector<int> *att1, vector<int> *att2)
{
	double x0,y0,z0;
	double r2,r6,sum,sumsum;
	int i1,j1;
	int n1,n2;
	int j;
	int jump;

	sumsum=0;
	for(i1=0;i1<(int)att1->size();i1++)
	for(j1=0;j1<(int)att2->size();j1++)
	{
		n1=att1->at(i1)-1;
		n2=att2->at(j1)-1;

		sum=0.0;
		for(j=0;j<nframe;j++)
		{
			jump=j*natom;
			x0=x.at(jump+n1)-x.at(jump+n2);
			y0=y.at(jump+n1)-y.at(jump+n2);
			z0=z.at(jump+n1)-z.at(jump+n2);
			r2=x0*x0+y0*y0+z0*z0;
			r6=r2*r2*r2;
			sum+=1.0/r6;
		}
		sum/=nframe;
		sumsum+=sum;
	}
	sumsum/=i1;
	sumsum/=j1;
	r2=pow(sumsum,-1/6.0);
	return r2;
}


void CTraj::evaluatenmrcons_frame(vector<struct noeline> *nmrcons, double cutoff)
{
	int ii,i,j;
	int pos1,pos2;
	int n1,n2;
	vector<int> att1,att2;
	double d;
	ofstream fout("noe_frame.dat");


	for(ii=0;ii<nframe;ii++)
	{
		for(i=0;i<(int)nmrcons->size();i++)
		{	

			if(fabs((float)(nmrcons->at(i).resid1-nmrcons->at(i).resid2))<=2)
				continue;

			nmrcons->at(i).obs.clear();
			nmrcons->at(i).pos1.clear();
			nmrcons->at(i).pos2.clear();
			for(n1=0;n1<(int)nmrcons->at(i).index1.atoms.size();n1++)
			for(n2=0;n2<(int)nmrcons->at(i).index2.atoms.size();n2++)
			{
				att1=nmrcons->at(i).index1.atoms.at(n1);
				att2=nmrcons->at(i).index2.atoms.at(n2);
				nmrcons->at(i).obs.push_back(noedistance_frame(&att1,&att2,ii));
				nmrcons->at(i).pos1.push_back(att1.at(0));
				nmrcons->at(i).pos2.push_back(att2.at(0));
			}

			d=100000.0;
			pos1=pos2=0;
			for(j=0;j<nmrcons->at(i).obs.size();j++)
			{
				if(nmrcons->at(i).obs.at(j)<d)
				{
					d=nmrcons->at(i).obs.at(j);
					pos1=nmrcons->at(i).pos1.at(j);
					pos2=nmrcons->at(i).pos2.at(j);
				}
			}
			nmrcons->at(i).bvio=0;
			nmrcons->at(i).d=d;
			
			if(nmrcons->at(i).a+cutoff<d)
				fout<<d<<" ";
			else
				fout<<d<<" ";
		}
		fout<<endl;

	}
	fout.close();

	ofstream fout2("noe_frame_head.dat");
	for(i=0;i<(int)nmrcons->size();i++)
	{
		if(fabs((float)(nmrcons->at(i).resid1-nmrcons->at(i).resid2))<=2)
			continue;

		fout2<<nmrcons->at(i).resid1<<" ";
		fout2<<nmrcons->at(i).resname1<<" ";
		fout2<<nmrcons->at(i).atomname1<<" ";
		fout2<<nmrcons->at(i).resid2<<" ";
		fout2<<nmrcons->at(i).resname2<<" ";
		fout2<<nmrcons->at(i).atomname2<<" ";
		fout2<<nmrcons->at(i).b<<" ";
		fout2<<nmrcons->at(i).c<<" ";
		fout2<<nmrcons->at(i).a<<endl;
	}
	fout2.close();
	return;
}

void CTraj::evulatenmrcons(vector<struct noeline> *nmrcons, double cutoff)
{
	int i,j;
	int pos1,pos2;
	int n1,n2;
	vector<int> att1,att2;
	double d;
	ofstream fout("violation.dat");
	ofstream fout2("fullfilled.dat");
	ofstream fout3("grouped.dat");
	ofstream fchimera("for_chimera.dat");
	ofstream fchimera2("for_chimera2.dat");


	
	for(i=0;i<(int)nmrcons->size();i++)
	{	
		for(n1=0;n1<(int)nmrcons->at(i).index1.atoms.size();n1++)
		for(n2=0;n2<(int)nmrcons->at(i).index2.atoms.size();n2++)
		{
			att1=nmrcons->at(i).index1.atoms.at(n1);
			att2=nmrcons->at(i).index2.atoms.at(n2);
			nmrcons->at(i).obs.push_back(noedistance(&att1,&att2));
			nmrcons->at(i).pos1.push_back(att1.at(0));
			nmrcons->at(i).pos2.push_back(att2.at(0));
		}
	}

	for(i=0;i<(int)nmrcons->size();i++)
	{
		d=100000.0;
		pos1=pos2=0;
		for(j=0;j<nmrcons->at(i).obs.size();j++)
		{
			if(nmrcons->at(i).obs.at(j)<d)
			{
				d=nmrcons->at(i).obs.at(j);
				pos1=nmrcons->at(i).pos1.at(j);
				pos2=nmrcons->at(i).pos2.at(j);
			}
		}
		nmrcons->at(i).bvio=0;
		nmrcons->at(i).d=d;

		if(nmrcons->at(i).a+cutoff<d)
		{
			nmrcons->at(i).bvio=1;
			fout<<nmrcons->at(i).group<<" ";
			fout<<i<<" "<<nmrcons->at(i).resid1<<" "<<nmrcons->at(i).resname1<<" "<<nmrcons->at(i).atomname1<<" ";
			fout<<nmrcons->at(i).resid2<<" "<<nmrcons->at(i).resname2<<" "<<nmrcons->at(i).atomname2<<" ";
			fout<<d<<" "<<nmrcons->at(i).c<<" "<<nmrcons->at(i).a<<endl;
			if(abs(nmrcons->at(i).resid1-nmrcons->at(i).resid2)>2 && d-nmrcons->at(i).a>cutoff)
			{
				fchimera<<"distance @/serialNumber="<<pos1<<" @/serialNumber="<<pos2<<endl;
				fchimera<<"display :"<<nmrcons->at(i).resid1<<endl;
				fchimera<<"display :"<<nmrcons->at(i).resid2<<endl;
				fchimera2<<"distance :"<<nmrcons->at(i).resid1<<"@CA :"<<nmrcons->at(i).resid2<<"@CA"<<endl;
			}
		}
		else
		{		
			fout2<<nmrcons->at(i).group<<" ";
			fout2<<i<<" "<<nmrcons->at(i).resid1<<" "<<nmrcons->at(i).resname1<<" "<<nmrcons->at(i).atomname1<<" ";
			fout2<<nmrcons->at(i).resid2<<" "<<nmrcons->at(i).resname2<<" "<<nmrcons->at(i).atomname2<<" ";
			fout2<<d<<" "<<nmrcons->at(i).c<<" "<<nmrcons->at(i).a<<endl;
		}
	}

	int ngroup=nmrcons->at(i-1).group;
	int i1=0;
	int i2=0;


	for(i=0;i<=ngroup;i++)
	{
		for(;i1<nmrcons->size() && nmrcons->at(i1).group<=i;i1++)
		{
			if(nmrcons->at(i1).a+cutoff<nmrcons->at(i1).d)
			{
				fout3<<nmrcons->at(i1).group<<" ";
				fout3<<i1<<" "<<nmrcons->at(i1).resid1<<" "<<nmrcons->at(i1).resname1<<" "<<nmrcons->at(i1).atomname1<<" ";
				fout3<<nmrcons->at(i1).resid2<<" "<<nmrcons->at(i1).resname2<<" "<<nmrcons->at(i1).atomname2<<" ";
				fout3<<nmrcons->at(i1).d<<" "<<nmrcons->at(i1).c<<" "<<nmrcons->at(i1).a<<" VIOLATED "<<endl;
			}
		}
		

		for(;i2<nmrcons->size() && nmrcons->at(i2).group<=i;i2++)
		{
			if(nmrcons->at(i2).a+cutoff>=nmrcons->at(i2).d)
			{
				fout3<<nmrcons->at(i2).group<<" ";
				fout3<<i2<<" "<<nmrcons->at(i2).resid1<<" "<<nmrcons->at(i2).resname1<<" "<<nmrcons->at(i2).atomname1<<" ";
				fout3<<nmrcons->at(i2).resid2<<" "<<nmrcons->at(i2).resname2<<" "<<nmrcons->at(i2).atomname2<<" ";
				fout3<<nmrcons->at(i2).d<<" "<<nmrcons->at(i2).c<<" "<<nmrcons->at(i2).a<<endl;
			}
		}
		fout3<<endl;
	}
	
	fout.close();
	fout2.close();
	fout3.close();
	fchimera.close();
	fchimera2.close();


	return;
}


void CTraj::rmsd_matrix(vector< vector<double> > *rmsd,vector<int> *ca, int skip)
{
	float *x1,*y1,*z1,*x2,*y2,*z2;
	int nca,jumpi,jumpj;
	int i,j,k;
	vector< double> t;
	class CRmsd rmsdf;

	nca=ca->size();
	
	x1=new float[nframe*nca];
	x2=new float[nframe*nca];
	y1=new float[nframe*nca];
	y2=new float[nframe*nca];
	z1=new float[nframe*nca];
	z2=new float[nframe*nca];

	for(i=0;i<nframe;i++)
	{
		if((i+1)%skip==0)
		{
			
			t.clear();
			jumpi=i*natom-1;
			for(j=0;j<i;j++)
				t.push_back(rmsd->at(j).at(i));
			t.push_back(0.0);
			for(j=i+1;j<nframe;j++)
			{
				if((j+1)%skip==0)
				{
					jumpj=j*natom-1;
					for(k=0;k<nca;k++)
					{
						x1[k]=x[jumpi+ca->at(k)];
						y1[k]=y[jumpi+ca->at(k)];
						z1[k]=z[jumpi+ca->at(k)];

						x2[k]=x[jumpj+ca->at(k)];
						y2[k]=y[jumpj+ca->at(k)];
						z2[k]=z[jumpj+ca->at(k)];
					}

					t.push_back(rmsdf.calculate_rotation_rmsd(x1,y1,z1,x2,y2,z2,nca));
				}
			}
			rmsd->push_back(t);
		}
	}
	return;
}



void CTraj::getvector(vector<struct index_three> nh,vector<double> *xx,vector<double> *yy,vector<double> *zz)
{
	int i,m,j;
	double r;
	double x1,y1,z1;
	for(i=0;i<nframe;i++)
	{	
		m=i*natom-1;		
		for(j=0;j<nh.size();j++)
		{
			x1=x.at(m+nh.at(j).x3)-x.at(m+nh.at(j).x2);
			y1=y.at(m+nh.at(j).x3)-y.at(m+nh.at(j).x2);
			z1=z.at(m+nh.at(j).x3)-z.at(m+nh.at(j).x2);
			r=sqrt(x1*x1+y1*y1+z1*z1);
			if(r>0)
			{
				x1/=r;y1/=r;z1/=r;
			}	
			else
			{
				x1=1.0;y1=z1=0.0;
			}	
			xx->push_back(x1);
			yy->push_back(y1);
			zz->push_back(z1);		
		}
	}
}


void CTraj:: getcoor(int ipos,int iframe,double *xx,double *yy,double *zz)
{
	ipos+=iframe*natom-1;
	*xx=x.at(ipos);
	*yy=y.at(ipos);
	*zz=z.at(ipos);
}


void CTraj::getcoor(vector<int> pos,int iframe,vector<double> *xx,vector<double> *yy,vector<double> *zz)
{
	int j;
	int adj;
	unsigned int i;

	adj=iframe*natom-1;
	for(i=0;i<pos.size();i++)
	{
		j=pos.at(i)+adj;
		xx->push_back(x.at(j));
		yy->push_back(y.at(j));
		zz->push_back(z.at(j));
	}
	return;
}


void CTraj::getcoor(vector<int> pos,vector<float> *xx,vector<float> *yy,vector<float> *zz)
{
	int j;
	int adj;
	unsigned int i,iframe;

	for(iframe=0;iframe<nframe;iframe++)
	{
		adj=iframe*natom-1;
		for(i=0;i<pos.size();i++)
		{
			j=pos.at(i)+adj;
			xx->push_back(x.at(j));
			yy->push_back(y.at(j));
			zz->push_back(z.at(j));
		}
	}
	return;
}



// Updated function for OpenACC
void CTraj::get_contact(float rc,float shift, vector<int> pos, vector<int> used, vector<float> * result)
{	
	//convert used to array
	int *used_arr = used.data();
	int used_size = used.size();
	//convert pos to array
	int *pos_arr = pos.data();
	int pos_size = pos.size();
	//convert result to array
	float *result_arr = result->data();
	int result_size = result->size();	

	double start_time = omp_get_wtime();
	int i,j;
	int ii,jj;
	float contact;
	float x0,y0,z0;
	float rr;

	double *this_x_arr = x_arr;
	double *this_y_arr = y_arr;
	double *this_z_arr = z_arr;
	
#pragma acc parallel copyin (used_arr[0:used_size], pos_arr[0:pos_size]) copyout (result_arr[0:result_size]) default(present)
{
	#pragma acc loop gang private(ii, x0, y0, z0)
	for(i=0;i<pos_size;i++) // pos.size() = 16 still worth unraveling?
	{
		contact=0.0;
		ii = pos_arr[i];
		//ii=pos.at(i);
		if(ii<0)
		{
			//result->push_back(-1.0);
			result_arr[i] = -1.0;
			continue;
		}
		ii--;

		//x0=x.at(ii);
                //y0=y.at(ii);
                //z0=z.at(ii);
               
		//x0 = x_arr[ii];
		//y0 = y_arr[ii];// do x0, y0 need to be pointers?
		//z0 = z_arr[ii]; 
		
		x0 = this_x_arr[ii];
		y0 = this_y_arr[ii];
		z0 = this_z_arr[ii];
		
		#pragma acc loop vector private(jj, rr) reduction(+: contact)
		for(j=0;j<used_size;j++)
		{
			jj=used_arr[j];
			if(jj<0)
				continue;
			jj--;
			//rr=(x.at(jj)-x0)*(x.at(jj)-x0)+(y.at(jj)-y0)*(y.at(jj)-y0)+(z.at(jj)-z0)*(z.at(jj)-z0);
			//rr=(x_arr[jj]-x0)*(x_arr[jj]-x0)+(y_arr[jj]-y0)*(y_arr[jj]-y0)+(z_arr[jj]-z0)*(z_arr[jj]-z0);	
			rr=(this_x_arr[jj]-x0)*(this_x_arr[jj]-x0)+(this_y_arr[jj]-y0)*(this_y_arr[jj]-y0)+(this_z_arr[jj]-z0)*(this_z_arr[jj]-z0);
			rr=sqrt(rr)-shift;
			contact+=exp(-rr/rc);				
		}

		//result->push_back(contact);
		result_arr[i] = contact;
	}
}
//	cout << (int)used.size() << endl; used size is around 1700 for kod and 50000 for tube_small, about half the atom size`
//	printf((int)pos.size()); //tomh
	cout << "get_contact: " << omp_get_wtime() - start_time << "seconds" << endl; //speeds: tubesmall: 4 seconds // 11.2 total seconds runtime
	return;
}


// New function for OpenACC
//void CTraj::get_all_contacts(vector<struct bb_group> *bb, vector<struct index_two> *index, int index_size, int *c2, int c2_size, float *results, int results_size)
//void CTraj::get_all_contacts(vector<struct bb_group> *bb, index_two *index, int index_size, int *c2, int c2_size, float *results, int results_size)
void CTraj::get_all_contacts(bb_group *bb, int bb_size, index_two *index, int index_size, int *c2, int c2_size, float *results, int results_size)
{

	double st = omp_get_wtime();

	//cout << "Get all contacts" << endl;
	// Variables
	int i, j;
	int ii1, ii2, ii3 ,jj;
	float contact1, contact2, contact3;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	float rr1,rr2,rr3;
	double xx,yy,zz;

	// Array containing all coordinates
	int c1_size = (index_size-2)*3;
	int *c1 = new int[c1_size];

	// Avoid copying "this" pointer
	//double *x_arr_this = x_arr;
	//double *y_arr_this = y_arr;
	//double *z_arr_this = z_arr;
	//int x_arr_size_this = x_size;
	//int y_arr_size_this = y_size;
	//int z_arr_size_this = z_size;
	//#pragma acc enter data copyin(results[0:results_size])
#pragma acc data create(c1[0:c1_size]) present(x_arr[0:x_size],y_arr[0:y_size],z_arr[0:z_size],c2[0:c2_size],results[0:results_size], \
	bb[0:bb_size], index[0:index_size], this)
{


	// Load up c1 same way as in predict_bb_static_ann
	//          index0           index1             index(index_size-2)
	// c1[ {coords at i=1}, {coords at i=2},... {coords at i=index_size-1} ]
	#pragma acc parallel loop independent
	for(i=0+1;i<index_size-1;i++)
	{
		if(index[i].x1 <= 0)
		{
			c1[((i-1)*3)+0]=-1;
			c1[((i-1)*3)+1]=-1;
			c1[((i-1)*3)+2]=-1;
		} else {
			c1[((i-1)*3)+0]=bb[index[i].x1-1].capos;
			c1[((i-1)*3)+1]=bb[index[i].x1-1].cbpos;
			c1[((i-1)*3)+2]=bb[index[i].x1-1].copos;
		}
	}


	//#pragma acc enter data copyin(c1[0:(index_size-2)*3],c2[0:c2_size],results[0:results_size])
	//#pragma acc enter data copyin(results[0:results_size])
	#pragma acc parallel loop independent private(ii1,ii2,ii3,x1,x2,x3,y1,y2,y3,z1,z2,z3)
	for(i=0+1;i<(int)index_size-1;i++)
	{
		contact1=0.0; contact2=0.0; contact3=0.0;

		ii1=c1[((i-1)*3)+0]; ii2=c1[((i-1)*3)+1]; ii3=c1[((i-1)*3)+2];

		if(ii1 < 0){ 
			x1=0; y1=0; z1=0;
		} else { 
			ii1--; x1=x_arr[ii1]; y1=y_arr[ii1]; z1=z_arr[ii1];
		}

		if(ii2 < 0){
			x2=0; y2=0; z2=0;
		} else {
			ii2--; x2=x_arr[ii2]; y2=y_arr[ii2]; z2=z_arr[ii2];
		}
		
		if(ii3 < 0){
			x3=0; y3=0; z3=0;
		} else {
			ii3--; x3=x_arr[ii3]; y3=y_arr[ii3]; z3=z_arr[ii3];
		}

		#pragma acc loop independent reduction(+:contact1) reduction(+:contact2) \
		reduction(+:contact3) private(jj,xx,yy,zz,rr1,rr2,rr3)
		for(j=0;j<c2_size;j++)
		{
			jj=c2[j];
			if(jj>=0){
				jj--;
				xx = x_arr[jj];
				yy = y_arr[jj];
				zz = z_arr[jj];
				rr1=(xx-x1)*(xx-x1)+(yy-y1)*(yy-y1)+(zz-z1)*(zz-z1);
				rr2=(xx-x2)*(xx-x2)+(yy-y2)*(yy-y2)+(zz-z2)*(zz-z2);
				rr3=(xx-x3)*(xx-x3)+(yy-y3)*(yy-y3)+(zz-z3)*(zz-z3);
				rr1=sqrt(rr1);
				rr2=sqrt(rr2);
				rr3=sqrt(rr3);
				contact1+=exp(-rr1/3.0);
				contact2+=exp(-rr2/3.0);
				contact3+=exp(-rr3/3.0);
			}				
		}
	
		if(ii1 < -1){
			results[((i-1)*3)+0]=-1.0;
		} else {
			results[((i-1)*3)+0]=contact1;
		}
		if(ii2 < -1){
			results[((i-1)*3)+1]=-1.0;
		} else {
			results[((i-1)*3)+1]=contact2;
		}if(ii3 < -1){
			results[((i-1)*3)+2]=-1.0;
		} else {
			results[((i-1)*3)+2]=contact3;
		}

	}

} // end data region
	//cout << "End get all contacts" << endl;
	delete(c1);
	cout << "get_all_contacts: " << omp_get_wtime() - st << " seconds" << endl;
}	


// Updated function for OpenACC
// Was replaced by better optmized version
void CTraj::get_contact(vector<int> pos, int* used, int used_size, vector<float> * result)
{
	int j;
	int ii1, ii2, ii3 ,jj;
	float contact1, contact2, contact3;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	float rr1,rr2,rr3;
	double xx,yy,zz;

//////////////////////////////////////////////////////////////


	int *pos_arr = pos.data();

	contact1=0.0;
	contact2=0.0;
	contact3=0.0;

	ii1=pos_arr[0];
	ii2=pos_arr[1];
	ii3=pos_arr[2];

	if(ii1 < 0){
		x1=0;
		y1=0;
		z1=0;
	} else {
		ii1--;
		x1=x[ii1];
		y1=y[ii1];
		z1=z[ii1];
	}
	if(ii2 < 0){
		x2=0;
		y2=0;
		z2=0;
	} else {
		ii2--;
		x2=x[ii2];
		y2=y[ii2];
		z2=z[ii2];
	}if(ii3 < 0){
		x3=0;
		y3=0;
		z3=0;
	} else {
		ii3--;
		x3=x[ii3];
		y3=y[ii3];
		z3=z[ii3];
	}

	double *x_arr_this = x_arr;
	double *y_arr_this = y_arr;
	double *z_arr_this = z_arr;

	#pragma acc parallel loop present(used[0:used_size],x_arr_this[0:x_size],y_arr_this[0:y_size],z_arr_this[0:z_size]) \
		reduction(+:contact1) reduction(+:contact2) reduction(+:contact3) private(jj,xx,yy,zz,rr1,rr2,rr3)
	for(j=0;j<used_size;j++)
	{
		jj=used[j];
		if(jj>=0){
			jj--;
			xx = x_arr_this[jj];
			yy = y_arr_this[jj];
			zz = z_arr_this[jj];
			rr1=(xx-x1)*(xx-x1)+(yy-y1)*(yy-y1)+(zz-z1)*(zz-z1);
			rr2=(xx-x2)*(xx-x2)+(yy-y2)*(yy-y2)+(zz-z2)*(zz-z2);
			rr3=(xx-x3)*(xx-x3)+(yy-y3)*(yy-y3)+(zz-z3)*(zz-z3);
			rr1=sqrt(rr1);
			rr2=sqrt(rr2);
			rr3=sqrt(rr3);
			contact1+=exp(-rr1/3.0);
			contact2+=exp(-rr2/3.0);
			contact3+=exp(-rr3/3.0);
		}				
	}
	
	if(ii1 < -1){
		result->push_back(-1.0);
	} else {
		result->push_back(contact1);
	}
	if(ii2 < -1){
		result->push_back(-1.0);
	} else {
		result->push_back(contact2);
	}if(ii3 < -1){
		result->push_back(-1.0);
	} else {
		result->push_back(contact3);
	}

//////////////////////////////////////////////////////////////

/*
	for(i=0;i<(int)pos.size();i++)
	{
		contact=0.0;
		ii=pos.at(i);
		if(ii<0)
		{
			result->push_back(-1.0);
			continue;
		}
		ii--;
		x0=x.at(ii);
		y0=y.at(ii);
		z0=z.at(ii);
		for(j=0;j<(int)used.size();j++)
		{
			jj=used.at(j);
			if(jj<0)
				continue;
			jj--;
			rr=(x.at(jj)-x0)*(x.at(jj)-x0)+(y.at(jj)-y0)*(y.at(jj)-y0)+(z.at(jj)-z0)*(z.at(jj)-z0);
			rr=sqrt(rr);
			contact+=exp(-rr/3.0);				
		}
		result->push_back(contact);
	}
*/
}



CTraj::CTraj()
{
	natom=0;
	nres=0;
	nframe=0;
};

CTraj::~CTraj()
{
	#pragma acc exit data delete(x_arr, y_arr, z_arr)
	#pragma acc exit data delete(this)
};


