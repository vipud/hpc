#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <cstring>


using namespace std;

#include "supply.h"
using namespace ldw_math;
 


void CCommandline::pharse(int argc, char** argv)
{
	int i,j;
	for(i=1;i<argc;i++)
	{
		for(j=0;j<narg;j++)
		{
			if(arguments.at(j).compare(argv[i])==0 && argv[i][0]=='-' && argv[i][1]=='-') //started with --, multi followings
			{
				parameters.at(j)=" ";			  // remove default one.
				while(i+1<argc && argv[i+1][0]!='-')  //store multi entries in to this value
				{
					parameters.at(j).append(argv[i+1]);
					parameters.at(j).append(" ");
					i++;
				}
			}
			else if(arguments.at(j).compare(argv[i])==0) // //started with -, normal one
			{
				if(i+1<argc && argv[i+1][0]!='-')
				{
					parameters.at(j)=argv[i+1];
					i++;
				}
				else
					parameters.at(j)="yes";
			}
		}
	}
	return;
}


void CCommandline::init(vector<string> in,vector<string> in2)
{
	int i;
	string t;
	
	narg=in.size();
	for(i=0;i<narg;i++)
	{
		if(in.at(i).at(0)=='-')
		{
			arguments.push_back(in.at(i));
			parameters.push_back(in2.at(i));
		}
		else
		{
			t="-";
			t+=in.at(i);
			arguments.push_back(t);
			parameters.push_back(in2.at(i));
		}
	}
	return;
};

void CCommandline::init(vector<string> in,vector<string> in2, vector<string> in3)
{
	int i;
	string t;
	
	narg=in.size();
	for(i=0;i<narg;i++)
	{
		if(in.at(i).at(0)=='-')
		{
			arguments.push_back(in.at(i));
			parameters.push_back(in2.at(i));
			informations.push_back(in3.at(i));
		}
		else
		{
			t="-";
			t+=in.at(i);
			arguments.push_back(t);
			parameters.push_back(in2.at(i));
			informations.push_back(in3.at(i));
		}
	}
	return;
};

string CCommandline::query(string in)
{
	string out;
	int i;

	out="no";
	for(i=0;i<narg;i++)
	{
		if(arguments.at(i).compare(in)==0)
			out=parameters.at(i);
	}

	return out;
}

void CCommandline::print(void)
{
	int i;
	printf("Command line arguments:\n");
	for(i=0;i<(int)arguments.size();i++)
	{
		printf("%-15s %15s",arguments.at(i).c_str(),parameters.at(i).c_str());
		if(i<(int)informations.size()) printf("     %-50s",informations.at(i).c_str());
		printf("\n");
	}
	printf("\n");
	return;
}


CCommandline::CCommandline() {};
CCommandline::~CCommandline() {};


void CRmsd::setup_rotation(float x[],float y[],float z[],
						float x0[],float y0[],float z0[], 
						int n_list)
	{
	  int n;
	  mov_com0 = mov_com1 = mov_com2 = 0.0f;
	  ref_com0 = ref_com1 = ref_com2 = 0.0f;
	  
	  for (n=0; n<n_list; n++) 
	  { 
		  mov_com0 += x0[n];
		  ref_com0 += x[n];
		  mov_com1 += y0[n];
		  ref_com1 += y[n];
		  mov_com2 += z0[n];
		  ref_com2 += z[n];
	  }
	    
	  mov_com0 /= n_list;
	  ref_com0 /= n_list;
	  mov_to_ref0 = ref_com0 - mov_com0;
	  mov_com1 /= n_list;
	  ref_com1 /= n_list;
	  mov_to_ref1 = ref_com1 - mov_com1;
	  mov_com2 /= n_list;
	  ref_com2 /= n_list;
	  mov_to_ref2 = ref_com2 - mov_com2;

	  for (n=0; n<n_list; n++) 
	  { 
		  x0[n] -= mov_com0;
		  x[n] -= ref_com0;
		  y0[n] -= mov_com1;
		  y[n] -= ref_com1;
		  z0[n] -= mov_com2;
		  z[n] -= ref_com2;
	  }

	  R00 = R01 = R02 = 0.0f;
	  R10 = R11 = R12 = 0.0f;
	  R20 = R21 = R22 = 0.0f;
	  E0 = 0.0f;

	  for (n=0; n<n_list; n++) 
	  {	
		  E0 +=x0[n]*x0[n]+x[n]*x[n]+y0[n]*y0[n]+y[n]*y[n]+z0[n]*z0[n]+z[n]*z[n];
		  R00 += x0[n] * x[n];
		  R01 += x0[n] * y[n];
		  R02 += x0[n] * z[n];
		  R10 += y0[n] * x[n];
		  R11 += y0[n] * y[n];
		  R12 += y0[n] * z[n];
		  R20 += z0[n] * x[n];
		  R21 += z0[n] * y[n];
		  R22 += z0[n] * z[n];
	  }
	  E0 *= 0.5f;

	  return;
	  };

int CRmsd::jacobi3(int* n_rot)
	{
	  int count;
	  float tresh, theta, tau, t, sum, s, h, g, c;
	  float b0,b1,b2;
	  float z0,z1,z2;
	  float gg,hh;

	  vec01=vec02=vec10=vec12=vec20=vec21=0.0f;
	  vec00=vec11=vec22=1.0f;
	  b0 = eigenval0 = RtR00;
	  b1 = eigenval1 = RtR11;
	  b2 = eigenval2 = RtR22;

	  z0=z1=z2=0.0f; 
	  *n_rot = 0;

	  /* 50 tries */
	  for (count=0; count<50; count++)     
	  {

		sum = 0.0f;
		sum += (float)fabs(RtR01);
		sum += (float)fabs(RtR02);
		sum += (float)fabs(RtR12);

		if (sum == 0.0) 
			return(1);

		if (count < 3) 
		  tresh = sum * 0.2f / 9.0f;    
		else       
		  tresh = 0.0f;      

		g = 100.0f * (float)fabs(RtR01);
		if ( count > 3  &&  fabs(eigenval0)+g == fabs(eigenval0) &&  fabs(eigenval1)+g == fabs(eigenval1) ) 
		{
			RtR01 = 0.0f;
		} 
		else if (fabs(RtR01) > tresh) 
		{
			h = eigenval1 - eigenval0;         
			if (fabs(h)+g == fabs(h))
				t = RtR01 / h; 
			else 
			{
				theta = 0.5f * h / (RtR01);
				t = 1.0f / ( (float)fabs(theta) + (float)sqrt(1.0f + theta*theta) );
				if (theta < 0.0f) 
					t = -t;
			}
			
			c = 1.0f / (float) sqrt(1 + t*t);
			s = t * c;
			tau = s / (1.0f + c);
			h = t * RtR01;

			z0 -= h;
			z1 += h;
			eigenval0 -= h;
			eigenval1 += h;

			RtR01 = 0.0f;
			gg=RtR02; hh=RtR12; RtR02=gg-s*(hh+gg*tau); RtR12=hh+s*(gg-hh*tau); //ROTATE(RtR, 0, 2, 1, 2)  
			gg=vec00; hh=vec01; vec00=gg-s*(hh+gg*tau); vec01=hh+s*(gg-hh*tau); //ROTATE(vec, 0, 0, 0, 1)
			gg=vec10; hh=vec11; vec10=gg-s*(hh+gg*tau); vec11=hh+s*(gg-hh*tau); //ROTATE(vec, 1, 0, 1, 1)
			gg=vec20; hh=vec21; vec20=gg-s*(hh+gg*tau); vec21=hh+s*(gg-hh*tau); //ROTATE(vec, 2, 0, 2, 1)
			++(*n_rot);
		}

		g = 100.0f * (float)fabs(RtR02);
		if ( count > 3  &&  fabs(eigenval0)+g == fabs(eigenval0) &&  fabs(eigenval2)+g == fabs(eigenval2) ) 
			RtR02 = 0.0f;
		else if (fabs(RtR02) > tresh) 
		{
			h = eigenval2 - eigenval0;
			if (fabs(h)+g == fabs(h))
				t = RtR02 / h;
			else
			{
				theta = 0.5f * h / (RtR02);
				t = 1.0f / ( (float)fabs(theta) + (float)sqrt(1.0f + theta*theta) );
				if (theta < 0.0f) 
					t = -t;
			}
			c = 1.0f / (float) sqrt(1 + t*t);
			s = t * c;
			tau = s / (1.0f + c);
			h = t * RtR02;

			z0 -= h;
			z2 += h;
			eigenval0 -= h;
			eigenval2 += h;

			RtR02 = 0.0f;
   			gg=RtR01; hh=RtR12; RtR01=gg-s*(hh+gg*tau); RtR12=hh+s*(gg-hh*tau);   
			gg=vec00; hh=vec02; vec00=gg-s*(hh+gg*tau); vec02=hh+s*(gg-hh*tau); 
			gg=vec10; hh=vec12; vec10=gg-s*(hh+gg*tau); vec12=hh+s*(gg-hh*tau); 
			gg=vec20; hh=vec22; vec20=gg-s*(hh+gg*tau); vec22=hh+s*(gg-hh*tau); 
			++(*n_rot);
		}

		g = 100.0f * (float)fabs(RtR12);
		if ( count > 3  &&  fabs(eigenval1)+g == fabs(eigenval1) &&  fabs(eigenval2)+g == fabs(eigenval2) )
			RtR12 = 0.0f;
		else if (fabs(RtR12) > tresh) 
		{
			h = eigenval2 - eigenval1;
			if (fabs(h)+g == fabs(h)) 
				t = RtR12 / h;
			else 
			{
				theta = 0.5f * h / (RtR12);
				t = 1.0f / ( (float)fabs(theta) + (float)sqrt(1.0f + theta*theta) );
				if (theta < 0.0f) 
					t = -t;
			}
	          
			c = 1.0f / (float) sqrt(1 + t*t);
			s = t * c;
			tau = s / (1.0f + c);
			h = t * RtR12;

			z1 -= h;
			z2 += h;
			eigenval1 -= h;
			eigenval2 += h;

			RtR12 = 0.0f;
   			gg=RtR01; hh=RtR02; RtR01=gg-s*(hh+gg*tau); RtR02=hh+s*(gg-hh*tau); //ROTATE(a, 0, 1, 0, 2)  
			gg=vec01; hh=vec02; vec01=gg-s*(hh+gg*tau); vec02=hh+s*(gg-hh*tau); //ROTATE(vec, 0, 1, 0, 2)
			gg=vec11; hh=vec12; vec11=gg-s*(hh+gg*tau); vec12=hh+s*(gg-hh*tau); //ROTATE(vec, 1, 1, 1, 2)
			gg=vec21; hh=vec22; vec21=gg-s*(hh+gg*tau); vec22=hh+s*(gg-hh*tau); //ROTATE(vec, 2, 1, 2, 2)
			++(*n_rot);
		}
		
		b0 += z0; eigenval0 = b0; z0 = 0.0f;
		b1 += z1; eigenval1 = b1; z1 = 0.0f;
		b2 += z2; eigenval2 = b2; z2 = 0.0f;
	  }

	  printf("Too many iterations in jacobi3\n");
	  return (0);
	};

int CRmsd::diagonalize_symmetric()
	{
		int k;
		int n_rot;
		float val; 
		
	  if (!jacobi3(&n_rot)) 
	  {
		  printf("convergence failed\n");
		  return (0);
	  }


	  k = 0;
	  val = eigenval0;
	  if (eigenval1 >= val)
	  { 
		  k = 1;
		  val = eigenval1;
	  }
	  if (eigenval2 >= val)
	  { 
		  k = 2;
		  val = eigenval2;
	  }

	  if (k != 0) 
	  {
		  if(k==1)
		  {
			  eigenval1 = eigenval0;
			  eigenval0 = val;
			  val = vec00;
			  vec00 = vec01;
			  vec01 = val;	
				  
			  val = vec10;
			  vec10 = vec11;
			  vec11 = val;
				  
			  val = vec20;
			  vec20 = vec21;
			  vec21 = val;
		  }

		  else
		  {
			  eigenval2 = eigenval0;
			  eigenval0 = val;
			  val = vec00;
			  vec00 = vec02;
			  vec02 = val;	
				  
			  val = vec10;
			  vec10 = vec12;
			  vec12 = val;
				  
			  val = vec20;
			  vec20 = vec22;
			  vec22 = val;
		  }

	  }


	  k = 1;
	  val = eigenval1;
	  if (eigenval2 >= val)
	  { 
		  k = 2;
		  val = eigenval2;
		  eigenval2 = eigenval1;
		  eigenval1 = val;
		  
		  val = vec01;
		  vec01 = vec02;
		  vec02 = val;
			  
		  val = vec11;
		  vec11 = vec12;
		  vec12 = val;
			  
		  val = vec21;
		  vec21 = vec22;
		  vec22 = val;
	  }

	  right_eigenvec00 = vec00;      right_eigenvec01 = vec10;      right_eigenvec02 = vec20;
	  right_eigenvec10 = vec01;      right_eigenvec11 = vec11;      right_eigenvec12 = vec21;
	  right_eigenvec20 = vec02;      right_eigenvec21 = vec12;      right_eigenvec22 = vec22;
	  
	  return (1);
	};

int CRmsd::calculate_rotation_matrix()
	{
	  float Rt00,Rt01,Rt02,Rt10,Rt11,Rt12,Rt20,Rt21,Rt22;
	  float v0,v1,v2;
	  float sigma;

	  Rt00=R00;Rt01=R10;Rt02=R20;Rt10=R01;Rt11=R11;Rt12=R21;Rt20=R02;Rt21=R12;Rt22=R22;

	  RtR00 = 0.0f;
	  RtR00 += Rt00 * R00;
	  RtR00 += Rt10 * R01;
	  RtR00 += Rt20 * R02;
	  RtR01 = 0.0f;
	  RtR01 += Rt00 * R10;
	  RtR01 += Rt10 * R11;
	  RtR01 += Rt20 * R12;
	  RtR02 = 0.0f;
	  RtR02 += Rt00 * R20;
	  RtR02 += Rt10 * R21;
	  RtR02 += Rt20 * R22;

	  RtR10 = 0.0f;
	  RtR10 += Rt01 * R00;
	  RtR10 += Rt11 * R01;
	  RtR10 += Rt21 * R02;
	  RtR11 = 0.0f;
	  RtR11 += Rt01 * R10;
	  RtR11 += Rt11 * R11;
	  RtR11 += Rt21 * R12;
	  RtR12 = 0.0f;
	  RtR12 += Rt01 * R20;
	  RtR12 += Rt11 * R21;
	  RtR12 += Rt21 * R22;

	  RtR20 = 0.0f;
	  RtR20 += Rt02 * R00;
	  RtR20 += Rt12 * R01;
	  RtR20 += Rt22 * R02;
	  RtR21 = 0.0f;
	  RtR21 += Rt02 * R10;
	  RtR21 += Rt12 * R11;
	  RtR21 += Rt22 * R12;
	  RtR22 = 0.0f;
	  RtR22 += Rt02 * R20;
	  RtR22 += Rt12 * R21;
	  RtR22 += Rt22 * R22;



	  if (!diagonalize_symmetric())
		return(0);

	 
	  right_eigenvec20 = -right_eigenvec11*right_eigenvec02 + right_eigenvec12*right_eigenvec01;
	  right_eigenvec21 = -right_eigenvec12*right_eigenvec00 + right_eigenvec10*right_eigenvec02;
	  right_eigenvec22 = -right_eigenvec10*right_eigenvec01 + right_eigenvec11*right_eigenvec00;
	  

	  left_eigenvec00 = right_eigenvec00*Rt00+right_eigenvec01*Rt01+right_eigenvec02*Rt02;
	  left_eigenvec01 = right_eigenvec00*Rt10+right_eigenvec01*Rt11+right_eigenvec02*Rt12;
	  left_eigenvec02 = right_eigenvec00*Rt20+right_eigenvec01*Rt21+right_eigenvec02*Rt22;
	      
	  left_eigenvec10 = right_eigenvec10*Rt00+right_eigenvec11*Rt01+right_eigenvec12*Rt02;
	  left_eigenvec11 = right_eigenvec10*Rt10+right_eigenvec11*Rt11+right_eigenvec12*Rt12;
	  left_eigenvec12 = right_eigenvec10*Rt20+right_eigenvec11*Rt21+right_eigenvec12*Rt22;
	      
	  left_eigenvec20 = right_eigenvec20*Rt00+right_eigenvec21*Rt01+right_eigenvec22*Rt02;
	  left_eigenvec21 = right_eigenvec20*Rt10+right_eigenvec21*Rt11+right_eigenvec22*Rt12;
	  left_eigenvec22 = right_eigenvec20*Rt20+right_eigenvec21*Rt21+right_eigenvec22*Rt22;
	      

	  float tteemmpp;
	  tteemmpp=(float)sqrt(left_eigenvec00*left_eigenvec00+left_eigenvec01*left_eigenvec01+left_eigenvec02*left_eigenvec02);	
	  left_eigenvec00/=tteemmpp;left_eigenvec01/=tteemmpp;left_eigenvec02/=tteemmpp;
	  tteemmpp=(float)sqrt(left_eigenvec10*left_eigenvec10+left_eigenvec11*left_eigenvec11+left_eigenvec12*left_eigenvec12);	
	  left_eigenvec10/=tteemmpp;left_eigenvec11/=tteemmpp;left_eigenvec12/=tteemmpp;
	  tteemmpp=(float)sqrt(left_eigenvec20*left_eigenvec20+left_eigenvec21*left_eigenvec21+left_eigenvec22*left_eigenvec22);	
	  left_eigenvec20/=tteemmpp;left_eigenvec21/=tteemmpp;left_eigenvec22/=tteemmpp;

	  
	  v0=left_eigenvec01*left_eigenvec12-left_eigenvec02*left_eigenvec11;
	  v1=left_eigenvec02*left_eigenvec10-left_eigenvec00*left_eigenvec12;
	  v2=left_eigenvec00*left_eigenvec11-left_eigenvec01*left_eigenvec10;

	  if(v0*left_eigenvec20+v1*left_eigenvec21+v2*left_eigenvec22< 0.0)
		sigma = -1.0f;
	  else 
		sigma = 1.0f;
	  
		
	 
	  left_eigenvec20 = v0; 
	  left_eigenvec21 = v1; 
	  left_eigenvec22 = v2; 

	  
	  
	  residual = E0 - (float) sqrt(fabs(eigenval0)) - (float) sqrt(fabs(eigenval1)) - sigma*(float) sqrt(fabs(eigenval2));
	  return (1);
	};


	void CRmsd::get_rotation_matrix()
	{
	  
	  U00 = 0.0f;
	  U00 += left_eigenvec00 * right_eigenvec00;
	  U00 += left_eigenvec10 * right_eigenvec10;
	  U00 += left_eigenvec20 * right_eigenvec20;
	  U01 = 0.0f;
	  U01 += left_eigenvec00 * right_eigenvec01;
	  U01 += left_eigenvec10 * right_eigenvec11;
	  U01 += left_eigenvec20 * right_eigenvec21;
	  U02 = 0.0f;
	  U02 += left_eigenvec00 * right_eigenvec02;
	  U02 += left_eigenvec10 * right_eigenvec12;
	  U02 += left_eigenvec20 * right_eigenvec22;

	  U10 = 0.0f;
	  U10 += left_eigenvec01 * right_eigenvec00;
	  U10 += left_eigenvec11 * right_eigenvec10;
	  U10 += left_eigenvec21 * right_eigenvec20;
	  U11 = 0.0f;
	  U11 += left_eigenvec01 * right_eigenvec01;
	  U11 += left_eigenvec11 * right_eigenvec11;
	  U11 += left_eigenvec21 * right_eigenvec21;
	  U12 = 0.0f;
	  U12 += left_eigenvec01 * right_eigenvec02;
	  U12 += left_eigenvec11 * right_eigenvec12;
	  U12 += left_eigenvec21 * right_eigenvec22;

	  U20 = 0.0f;
	  U20 += left_eigenvec02 * right_eigenvec00;
	  U20 += left_eigenvec12 * right_eigenvec10;
	  U20 += left_eigenvec22 * right_eigenvec20;
	  U21 = 0.0f;
	  U21 += left_eigenvec02 * right_eigenvec01;
	  U21 += left_eigenvec12 * right_eigenvec11;
	  U21 += left_eigenvec22 * right_eigenvec21;
	  U22 = 0.0f;
	  U22 += left_eigenvec02 * right_eigenvec02;
	  U22 += left_eigenvec12 * right_eigenvec12;
	  U22 += left_eigenvec22 * right_eigenvec22; 
		return;
	};

		

	float CRmsd::calculate_rotation_rmsd(float x[],float y[],float z[],float x0[],float y0[],float z0[],int n_list)
	{
		setup_rotation(x,y,z,x0,y0,z0,n_list);
		calculate_rotation_matrix();	  
		residual = (float)fabs(residual);
		return (float)sqrt( fabs((float) (residual)*2.0/((float)n_list)));
	};

	float CRmsd::direct_rmsd(float x[],float y[],float z[],float x0[],float y0[],float z0[],int n_list)
	{ 
		int n;
		float tx,ty,tz,d;

		//move to center
		for (n=0; n<n_list; n++) 
		{
			x0[n] -= mov_com0;
			x[n] -= ref_com0;
			y0[n] -= mov_com1;
			y[n] -= ref_com1;
			z0[n] -= mov_com2;
			z[n] -= ref_com2;
		}

		get_rotation_matrix();

		//rotate it
		for(n=0;n<n_list;n++)
		{
			tx=x0[n]*U00+y0[n]*U01+z0[n]*U02;
			ty=x0[n]*U10+y0[n]*U11+z0[n]*U12;
			tz=x0[n]*U20+y0[n]*U21+z0[n]*U22;
			x0[n]=tx;
			y0[n]=ty;
			z0[n]=tz;
		}

		//calculate RMSD directly
		d=0;
		for(n=0;n<n_list;n++)
		{
			tx=x0[n]-x[n];
			ty=y0[n]-y[n];
			tz=z0[n]-z[n];
			d+=tx*tx+ty*ty+tz*tz;
		}
		d=d/(n_list);
		d=sqrt(d);

		return d;
		
	};

	CRmsd::CRmsd() {};
	CRmsd::~CRmsd() {};




namespace ldw_math
{
	float veclengthf(float x[3])
	{
			return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	}

	double veclength(double x[3])
	{
			return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	}

	void cross(double z[3],double x[3],double y[3])
	{
			z[0]=x[1]*y[2]-x[2]*y[1];
			z[1]=-x[0]*y[2]+x[2]*y[0];
			z[2]=x[0]*y[1]-x[1]*y[0];
			return;
	}

	double dot(double x[3],double y[3])
	{
			return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
	}

	double coor_to_angle(double x2,double y2,double z2,double x3,double y3,double z3,double x4,double y4,double z4)
	{
			double angle;
			double b[3],c[3];

	        
			b[0]=x3-x2;b[1]=y3-y2;b[2]=z3-z2;
			c[0]=x4-x3;c[1]=y4-y3;c[2]=z4-z3;
			angle=dot(b,c)/veclength(b)/veclength(c);
			if(angle>1.0) angle=1.0;
			if(angle<-1.0) angle=-1.0;
	/*        angle=sqrt(1-angle*angle);*/
			return angle;
	}

	double coor_to_length(double x3,double y3,double z3,double x4,double y4,double z4)
	{
			return (x4-x3)*(x4-x3)+(y4-y3)*(y4-y3)+(z4-z3)*(z4-z3);
	}

	double coor_to_dihe(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double x4,double y4,double z4)
	{
			double t;
			double angle;
			double a[3],b[3],c[3];
			double d[3],e[3];
			double f[3];

			a[0]=x2-x1;a[1]=y2-y1;a[2]=z2-z1;
			b[0]=x3-x2;b[1]=y3-y2;b[2]=z3-z2;
			c[0]=x4-x3;c[1]=y4-y3;c[2]=z4-z3;
			cross(d,a,b);
			cross(e,b,c);
			cross(f,d,e);


			angle=dot(d,e)/veclength(d)/veclength(e);
			if(angle>1.0) angle=1.0;
			if(angle<-1.0) angle=-1.0;
			angle=acos(angle);
			if((t=dot(b,f))<0)
					angle=-angle;
			return angle;
	}

	double mysign(double a,double b)
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

	double mymax(double a, double b)
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


	 
	double PYTHAG(double a, double b)
	{
		double at = fabs(a), bt = fabs(b), ct, result;

		if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
		else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
		else result = 0.0;
		return(result);
	}


	double area( double a, double b, double c )
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


	void fit(vector<double> x,vector<double> y, vector<double> *z,double *a,double *b,double *rms,double *r)
	{
		int i;
		int ndata;
		double t,t2,sxoss,syoss,sx,sy,stx,sty,stz,ss;
		double chi2;

		ndata=x.size();
		
		sx=sy=0.0;
		stx=sty=stz=0.0;
		*b=0;
		*a=0;
		
		for(i=0;i<ndata;i++)
		{
			sx+=x[i];
			sy+=y[i];
		}
		ss=ndata;
		
		sxoss=sx/ss;
		syoss=sy/ss;
		
		for(i=0;i<ndata;i++)
		{
			t=x[i]-sxoss;
			stx+=t*t;
			*b+=t*y[i];
			t2=y[i]-syoss;
			sty+=t2*t2;
			stz+=t*t2;		
		}
		
		*b/=stx;	
		*a=(sy-sx*(*b))/ss;
		
		chi2=0.0;
		for(i=0;i<ndata;i++)
		{
			t=*a+(*b)*x[i];
			z->at(i)=t;
			t=y[i]-t;
			chi2+=t*t;		
		}
		chi2/=ss;
		*rms=sqrt(chi2);	
		*r=stz/sqrt(stx*sty);	
	}

	int dsvd(double a[6][3], int m, int n, double *w, double v[3][3])
	{
		int flag, i, its, j, jj, k, l, nm;
		double c, f, h, s, x, y, z;
		double anorm = 0.0, g = 0.0, scale = 0.0;
		double *rv1;
	  
		if (m < n) 
		{
			fprintf(stderr, "#rows must be > #cols \n");
			return(0);
		}
	  
		rv1 = (double *)malloc((unsigned int) n*sizeof(double));

	/* Householder reduction to bidiagonal form */
		for (i = 0; i < n; i++) 
		{
			/* left-hand reduction */
			l = i + 1;
			rv1[i] = scale * g;
			g = s = scale = 0.0;
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
	                
					g = -mysign(sqrt(s), f);
	                
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
	    
			/* right-hand reduction */
			g = s = scale = 0.0;
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
					g = -mysign(sqrt(s), f);
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
			anorm = mymax(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
		}
	  
		/* accumulate the right-hand transformation */
		for (i = n - 1; i >= 0; i--) 
		{
			if (i < n - 1) 
			{
				if (g) 
				{
					for (j = l; j < n; j++)
						v[j][i] = (double)(((double)a[i][j] / (double)a[i][l]) / g);
						/* double division to avoid underflow */
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
	  
		/* accumulate the left-hand transformation */
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

		/* diagonalize the bidiagonal form */
		for (k = n - 1; k >= 0; k--) 
		{                             /* loop over singular values */
			for (its = 0; its < 30; its++) 
			{                         /* loop over allowed iterations */
				flag = 1;
				for (l = k; l >= 0; l--) 
				{                     /* test for splitting */
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
							h = PYTHAG(f, g);
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
				{                  /* convergence */
					if (z < 0.0) 
					{              /* make singular value nonnegative */
						w[k] = (double)(-z);
						for (j = 0; j < n; j++) 
							v[j][k] = (-v[j][k]);
					}
					break;
				}
				if (its >= 30) {
					free((void*) rv1);
					fprintf(stderr, "No convergence after 30,000! iterations \n");
					return(0);
				}
	    
				/* shift from bottom 2 x 2 minor */
				x = (double)w[l];
				nm = k - 1;
				y = (double)w[nm];
				g = rv1[nm];
				h = rv1[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = PYTHAG(f, 1.0);
				f = ((x - z) * (x + z) + h * ((y / (f + mysign(g, f))) - h)) / x;
	          
				/* next QR transformation */
				c = s = 1.0;
				for (j = l; j <= nm; j++) 
				{
					i = j + 1;
					g = rv1[i];
					y = (double)w[i];
					h = s * g;
					g = c * g;
					z = PYTHAG(f, h);
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
					z = PYTHAG(f, h);
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
			}
		}
		free((void*) rv1);
		return(1);
	}

	void ring(double x[6][3], int m, double ori[3])
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


			dsvd(xx, m, 3, w,v);

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



	int dsvd2(double *a, int m, int n, double *w, double v[3][3])
	{
		int flag, i, its, j, jj, k, l, nm;
		double c, f, h, s, x, y, z;
		double anorm = 0.0, g = 0.0, scale = 0.0;
		double *rv1;
	  
		if (m < n) 
		{
			fprintf(stderr, "#rows must be > #cols \n");
			return(0);
		}
	  
		rv1 = (double *)malloc((unsigned int) n*sizeof(double));

	/* Householder reduction to bidiagonal form */
		for (i = 0; i < n; i++) 
		{
			/* left-hand reduction */
			l = i + 1;
			rv1[i] = scale * g;
			g = s = scale = 0.0;
			if (i < m) 
			{
				for (k = i; k < m; k++) 
					scale += fabs((double)a[k*n+i]);
				if (scale) 
				{
					for (k = i; k < m; k++) 
					{
						a[k*n+i] = (double)((double)a[k*n+i]/scale);
						s += ((double)a[k*n+i] * (double)a[k*n+i]);
					}
					f = (double)a[i*n+i];
	                
					g = -mysign(sqrt(s), f);
	                
					h = f * g - s;
					a[i*n+i] = (double)(f - g);
					if (i != n - 1) 
					{
						for (j = l; j < n; j++) 
						{
							for (s = 0.0, k = i; k < m; k++) 
								s += ((double)a[k*n+i] * (double)a[k*n+j]);
							f = s / h;
							for (k = i; k < m; k++) 
								a[k*n+j] += (double)(f * (double)a[k*n+i]);
						}
					}
					for (k = i; k < m; k++) 
						a[k*n+i] = (double)((double)a[k*n+i]*scale);
				}
			}
			w[i] = (double)(scale * g);
	    
			/* right-hand reduction */
			g = s = scale = 0.0;
			if (i < m && i != n - 1) 
			{
				for (k = l; k < n; k++) 
					scale += fabs((double)a[i*n+k]);
				if (scale) 
				{
					for (k = l; k < n; k++) 
					{
						a[i*n+k] = (double)((double)a[i*n+k]/scale);
						s += ((double)a[i*n+k] * (double)a[i*n+k]);
					}
					f = (double)a[i*n+l];
					g = -mysign(sqrt(s), f);
					h = f * g - s;
					a[i*n+l] = (double)(f - g);
					for (k = l; k < n; k++) 
						rv1[k] = (double)a[i*n+k] / h;
					if (i != m - 1) 
					{
						for (j = l; j < m; j++) 
						{
							for (s = 0.0, k = l; k < n; k++) 
								s += ((double)a[j*n+k] * (double)a[i*n+k]);
							for (k = l; k < n; k++) 
								a[j*n+k] += (double)(s * rv1[k]);
						}
					}
					for (k = l; k < n; k++) 
						a[i*n+k] = (double)((double)a[i*n+k]*scale);
				}
			}
			anorm = mymax(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
		}
	  
		/* accumulate the right-hand transformation */
		for (i = n - 1; i >= 0; i--) 
		{
			if (i < n - 1) 
			{
				if (g) 
				{
					for (j = l; j < n; j++)
						v[j][i] = (double)(((double)a[i*n+j] / (double)a[i*n+l]) / g);
						/* double division to avoid underflow */
					for (j = l; j < n; j++) 
					{
						for (s = 0.0, k = l; k < n; k++) 
							s += ((double)a[i*n+k] * (double)v[k][j]);
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
	  
		/* accumulate the left-hand transformation */
		for (i = n - 1; i >= 0; i--) 
		{
			l = i + 1;
			g = (double)w[i];
			if (i < n - 1) 
				for (j = l; j < n; j++) 
					a[i*n+j] = 0.0;
			if (g) 
			{
				g = 1.0 / g;
				if (i != n - 1) 
				{
					for (j = l; j < n; j++) 
					{
						for (s = 0.0, k = l; k < m; k++) 
							s += ((double)a[k*n+i] * (double)a[k*n+j]);
						f = (s / (double)a[i*n+i]) * g;
						for (k = i; k < m; k++) 
							a[k*n+j] += (double)(f * (double)a[k*n+i]);
					}
				}
				for (j = i; j < m; j++) 
					a[j*n+i] = (double)((double)a[j*n+i]*g);
			}
			else 
			{
				for (j = i; j < m; j++) 
					a[j*n+i] = 0.0;
			}
			++a[i*n+i];
		}

		/* diagonalize the bidiagonal form */
		for (k = n - 1; k >= 0; k--) 
		{                             /* loop over singular values */
			for (its = 0; its < 30; its++) 
			{                         /* loop over allowed iterations */
				flag = 1;
				for (l = k; l >= 0; l--) 
				{                     /* test for splitting */
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
							h = PYTHAG(f, g);
							w[i] = (double)h; 
							h = 1.0 / h;
							c = g * h;
							s = (- f * h);
							for (j = 0; j < m; j++) 
							{
								y = (double)a[j*n+nm];
								z = (double)a[j*n+i];
								a[j*n+nm] = (double)(y * c + z * s);
								a[j*n+i] = (double)(z * c - y * s);
							}
						}
					}
				}
				z = (double)w[k];
				if (l == k) 
				{                  /* convergence */
					if (z < 0.0) 
					{              /* make singular value nonnegative */
						w[k] = (double)(-z);
						for (j = 0; j < n; j++) 
							v[j][k] = (-v[j][k]);
					}
					break;
				}
				if (its >= 30) {
					free((void*) rv1);
					fprintf(stderr, "No convergence after 30,000! iterations \n");
					return(0);
				}
	    
				/* shift from bottom 2 x 2 minor */
				x = (double)w[l];
				nm = k - 1;
				y = (double)w[nm];
				g = rv1[nm];
				h = rv1[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = PYTHAG(f, 1.0);
				f = ((x - z) * (x + z) + h * ((y / (f + mysign(g, f))) - h)) / x;
	          
				/* next QR transformation */
				c = s = 1.0;
				for (j = l; j <= nm; j++) 
				{
					i = j + 1;
					g = rv1[i];
					y = (double)w[i];
					h = s * g;
					g = c * g;
					z = PYTHAG(f, h);
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
					z = PYTHAG(f, h);
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
						y = (double)a[jj*n+j];
						z = (double)a[jj*n+i];
						a[jj*n+j] = (double)(y * c + z * s);
						a[jj*n+i] = (double)(z * c - y * s);
					}
				}
				rv1[l] = 0.0;
				rv1[k] = f;
				w[k] = (double)x;
			}
		}
		free((void*) rv1);
		return(1);
	}




	void regression_plane(double* x, int m, double ori[3])
	{
			int i,j;
			int id;
			double w[3];
			double v[3][3];
			double d;
			double *xx;


			xx=new double [m*3];

			for(i=0;i<m;i++)
			for(j=0;j<3;j++)
					xx[i*3+j]=x[i*3+j];


			dsvd2(xx, m, 3, w,v);

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


	void rotation_around_axis(double point[3], double ori[3], double theta)
	{
		double x,y,z;
		double u,v,w;
		double t1_dot,t2_cos,t3_sin;

		x=point[0];
		y=point[1];
		z=point[2];

		u=ori[0];
		v=ori[1];
		w=ori[2];

		t1_dot=u*x+v*y+w*z;
		t2_cos=cos(theta);
		t3_sin=sin(theta);

		point[0]=u*t1_dot*(1-t2_cos)+x*t2_cos+(-w*y+v*z)*t3_sin;
		point[1]=v*t1_dot*(1-t2_cos)+y*t2_cos+(+w*x-u*z)*t3_sin;
		point[2]=w*t1_dot*(1-t2_cos)+z*t2_cos+(-v*x+u*y)*t3_sin;

		return;
	}

	void project(double ori[3], double p1[3], double p2[3])
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

	double effect(double x[6][3], int m, double ori[3], double p1[3])
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
			project(ori,p1,p2);
			for(i=0;i<m-1;i++)
			{
					for(j=0;j<3;j++)
					{
							t1[j]=x[i][j]-p1[j];
							t2[j]=x[i+1][j]-x[i][j];
							t3[j]=x[i+1][j]-p1[j];
					}
					leg1=veclength(t1);
					leg3=veclength(t3);
					d=1/(leg1*leg1*leg1)+1/(leg3*leg3*leg3);
					cross(t4,t2,t1);
					tt=dot(t4,ori);

					for(j=0;j<3;j++)
					{
							t1[j]=x[i][j]-p2[j];
							t3[j]=x[i+1][j]-p2[j];
					}
					leg1=veclength(t1);
					leg2=veclength(t2);
					leg3=veclength(t3);
					s=area(leg1,leg2,leg3);
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
					leg1=veclength(t1);
					leg3=veclength(t3);
					d=1/(leg1*leg1*leg1)+1/(leg3*leg3*leg3);
					cross(t4,t2,t1);
					tt=dot(t4,ori);

					for(j=0;j<3;j++)
					{
							t1[j]=x[i][j]-p2[j];
							t3[j]=x[0][j]-p2[j];
					}
					leg1=veclength(t1);
					leg2=veclength(t2);
					leg3=veclength(t3);
					s=area(leg1,leg2,leg3);
					if(tt<0)
							s=-s;
					ss+=s*d;
			}

			return ss;
	}



	double gaussrand()
	{
		static double V1, V2, S;
		static int phase = 0;
		double X;

		if(phase == 0) {
			do {
				double U1 = (double)rand() / RAND_MAX;
				double U2 = (double)rand() / RAND_MAX;

				V1 = 2 * U1 - 1;
				V2 = 2 * U2 - 1;
				S = V1 * V1 + V2 * V2;
				} while(S >= 1 || S == 0);

			X = V1 * sqrt(-2 * log(S) / S);
		} else
			X = V2 * sqrt(-2 * log(S) / S);

		phase = 1 - phase;

		return X;
	}



vector<int> cluster_pick2(int ndata,vector<double> hs, vector<double> ns, double scale)
{
	vector<int> back;
	int i;
	int b1,b2,b3,b4,fb1,fb2,fb3,fb4;
	double min_dis;
	double d1,d2,d3,d4;
	vector<int> touse;

	if(ndata==3)
	{
		min_dis=1000000.0;
		for(b1=0;b1<ndata;b1++)
		{
			for(b2=b1+1;b2<ndata;b2++)
			{
				d1=(hs.at(b1)-hs.at(b2))*scale;
				d2=ns.at(b1)-ns.at(b2);
				d1=d1*d1+d2*d2;
				if(d1<min_dis)
				{
					min_dis=d1;
					fb1=b1;
					fb2=b2;
				}
			}
		}
		back.resize(ndata,0);
		back.at(fb1)=1;
		back.at(fb2)=1;
	}


	else if(ndata==4)
	{
		min_dis=1000000.0;
		for(b1=0;b1<ndata;b1++)
		{
			for(b2=b1+1;b2<ndata;b2++)
			{
				d1=(hs.at(b1)-hs.at(b2))*scale;
				d2=ns.at(b1)-ns.at(b2);
				d1=d1*d1+d2*d2;


				touse.clear();
				for(i=0;i<ndata;i++)
				{
					if(i!=b1 && i!=b2)
						touse.push_back(i);
				}

				b3=touse.at(0);
				b4=touse.at(1);
				d3=(hs.at(b3)-hs.at(b4))*scale;
				d4=ns.at(b3)-ns.at(b4);
				d1=d3*d3+d4*d4+d1;

				if(d1<min_dis)
				{
					min_dis=d1;
					fb1=b1;
					fb2=b2;
					fb3=b3;
					fb4=b4;
				}
			}
		}
		back.resize(ndata,0);
		back.at(fb1)=1;
		back.at(fb2)=1;
	}

	return back;
}



// ///////////////////////////////////////////////////////////////////////////////////////////////////////////

}//end of ldw_math



namespace Sequence
{
	string code2name(char c)
	{
		string t;
		t="UNK";
		switch(c)
		{
		case 'A':
			t="ALA";
			break;
		case 'G':
			t="GLY";
			break;
		case 'V':
			t="VAL";
			break;
		case 'L':
			t="LEU";
			break;
		case 'M':
			t="MET";
			break;
		case 'I':
			t="ILE";
			break;
		case 'S':
			t="SER";
			break;
		case 'T':
			t="THR";
			break;
		case 'C':
			t="CYS";
			break;
		case 'P':
			t="PRO";
			break;
		case 'N':
			t="ASN";
			break;
		case 'Q':
			t="GLN";
			break;
		case 'F':
			t="PHE";
			break;
		case 'Y':
			t="TYR";
			break;
		case 'W':
			t="TRP";
			break;
		case 'K':
			t="LYS";
			break;
		case 'R':
			t="ARG";
			break;
		case 'H':
			t="HIS";
			break;
		case 'D':
			t="ASP";
			break;
		case 'E':
			t="GLU";
			break;
		case 'X':
			t="UNK";
			break;
		}

		return t;
	}


	int code2pos(char code)
	{ 
		int pos;
		switch (code)
		{
		case 'A':
			pos=0;
			break;
		case 'R':
			pos=1;
			break;
		case 'N':
			pos=2;
			break;
		case 'D':
			pos=3;
			break;
		case 'C':
			pos=4;
			break;
		case 'Q':
			pos=5;
			break;
		case 'E':
			pos=6;
			break;
		case 'G':
			pos=7;
			break;
		case 'H':
			pos=8;
			break;
		case 'I':
			pos=9;
			break;
		case 'L':
			pos=10;
			break;
		case 'K':
			pos=11;
			break;
		case 'M':
			pos=12;
			break;
		case 'F':
			pos=13;
			break;
		case 'P':
			pos=14;
			break;
		case 'S':
			pos=15;
			break;
		case 'T':
			pos=16;
			break;
		case 'W':
			pos=17;
			break;
		case 'Y':
			pos=18;
			break;
		case 'V':
			pos=19;
			break;
		default:
			pos=7;
			//cout<<"unknown AA code "<<code<<endl;
		}

		return pos;
	}


	char name2code(string in)
	{
		char c;
		if(in.compare("ALA")==0)
			c='A';
		else if(in.compare("GLY")==0)
			c='G';
		else if(in.compare("VAL")==0)
			c='V';	
		else if(in.compare("LEU")==0)
			c='L';
		else if(in.compare("MET")==0)
			c='M';
		else if(in.compare("ILE")==0)
			c='I';
		else if(in.compare("SER")==0)
			c='S';
		else if(in.compare("THR")==0)
			c='T';
		else if(in.compare("CYS")==0)
			c='C';
		else if(in.compare("CYX")==0)
			c='C';
		else if(in.compare("PRO")==0)
			c='P';
		else if(in.compare("ASN")==0)
			c='N';
		else if(in.compare("GLN")==0)
			c='Q';
		else if(in.compare("PHE")==0)
			c='F';
		else if(in.compare("TYR")==0)
			c='Y';
		else if(in.compare("TRP")==0)
			c='W';
		else if(in.compare("HIS")==0)
			c='H';
		else if(in.compare("HIE")==0)
			c='H';
		else if(in.compare("HID")==0)
			c='H';
		else if(in.compare("ARG")==0)
			c='R';
		else if(in.compare("LYS")==0)
			c='K';
		else if(in.compare("ASP")==0)
			c='D';
		else if(in.compare("GLU")==0)
			c='E';
		else if(in.compare("UNK")==0)
			c='X';
		else if(in.compare("MIS")==0)
			c='X';
		else
			c='B';

		return c;
	}

	void code2array(char code, int buffer[20])
	{
		int pos;

		for(pos=0;pos<20;pos++)
			buffer[pos]=0;	
		pos=code2pos(code);
		buffer[pos]=1;
		return;
	}



	void code2same(char code, int buffer[20])
	{
		int pos;
		int i;
		
		static int blosum[20][20]=
		{{4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0},
		{-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3},
		{-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3},
		{-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3},
		{0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},
		{-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2},
		{-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2},
		{0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3},
		{-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3},
		{-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3},
		{-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1},
		{-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2},
		{-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1},
		{-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1},
		{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2},
		{1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2},
		{0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0},
		{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3},
		{-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1},
		{0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4}};

		pos=code2pos(code);
		for(i=0;i<20;i++)
			buffer[i]=blosum[pos][i];

		return;
	}


	vector<double> expand(char code,vector<double> *in)
	{
		int i,pos,n;
		vector<double> out;

		pos=code2pos(code); //pos =0 ... 19
		n=in->size();

		out.resize(n*20);

		for(i=0;i<n;i++)
		{
			out.at(i+pos*n)=in->at(i);
		}
		return out;
	}


	vector<double> expand2(int type,vector<double> *in)
	{
		vector<double> out;
		int pos,n;
		int i;

		pos=type-1; //pos=0,1,2,....9
		n=in->size();

		out.resize(n*10);
		for(i=0;i<n;i++)
			out.at(i+pos*n)=in->at(i);

		return out;
	}




	vector<int> align(string c1,string c2)
	//c1: sequence of MD run. c2: sequence from bmrb.
	{
		int s1,s2;
		int i,j;
		int *score,*point;
		int score1,score2,score3;
		vector<int> out;

		s1=c1.length();
		s2=c2.length();

		out.resize(s1+1);

		score=new int [(s1+1)*(s2+1)];
		point=new int [(s1+1)*(s2+1)];

		for(i=0;i<=s1;i++)
		for(j=0;j<=s2;j++)
		{
			score[j*(s1+1)+i]=0;
			point[j*(s1+1)+i]=0;
		}


		for(i=1;i<=s1;i++)
		{
			score[i]=-i;
			point[i]=3;
		}
		for(j=1;j<=s2;j++)
		{
			score[j*(s1+1)]=-j;	
			point[j*(s1+1)]=2;
		}

		for(i=1;i<=s1;i++)
		{
			for(j=1;j<=s2;j++)
			{
				if(c1[i-1] == c2[j-1] || c1[i-1]=='X' || c2[j-1]=='X') // 'X' can match any AA
					score1=score[(j-1)*(s1+1)+i-1]+1;
				else
					score1=score[(j-1)*(s1+1)+i-1]-1;
				score2=score[(j-1)*(s1+1)+i]-1;
				score3=score[j*(s1+1)+i-1]-1;

				if(score1>score2)
				{
					if(score1>score3)
					{
						score[j*(s1+1)+i]=score1;
						point[j*(s1+1)+i]=1;
					}
					else
					{
						score[j*(s1+1)+i]=score3;
						point[j*(s1+1)+i]=3;
					}
				}
				else
				{
					if(score2>score3)
					{
						score[j*(s1+1)+i]=score2;
						point[j*(s1+1)+i]=2;
					}
					else
					{
						score[j*(s1+1)+i]=score3;
						point[j*(s1+1)+i]=3;
					}
				}
			}
		}


	/*	for(i=0;i<=s1;i++)
		{
			for(j=0;j<=s2;j++)
			{
				cout<<score[j*(s1+1)+i]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;

		for(i=0;i<=s1;i++)
		{
			for(j=0;j<=s2;j++)
			{
				cout<<point[j*(s1+1)+i]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;*/


		i=s1;
		j=s2;
		while(point[j*(s1+1)+i]!=0)
		{
			if(point[j*(s1+1)+i]==1)
			{
				//cout<<c1[i-1];
				//cout<<c2[j-1];
				if(c1[i-1]==c2[j-1] || c1[i-1]=='X' || c2[j-1]=='X' ) 
					out[i]=j;
				i--;
				j--;
				//cout<<endl;
			}
			else if(point[j*(s1+1)+i]==2)
			{
				//cout<<"-";
				//cout<<c2[j-1];
				j--;
				//cout<<endl;
			}
			else if(point[j*(s1+1)+i]==3)
			{
				//cout<<c1[i-1];
				//cout<<"-";
				i--;
				//cout<<endl;
			}
		}


		/*for(i=1;i<out.size();i++)
		{
			cout<<i<<" "<<out.at(i)<<endl;
		}*/

		delete [] score;
		delete [] point;
		
		return out;
	}


	vector<int> aligno(string c1,string c2, string &out1, string &out2, string &out3)
	//c1: sequence of MD run. c2: sequence from bmrb.
	{
		int s1,s2;
		int i,j;
		int *score,*point;
		int score1,score2,score3;
		vector<int> out;

		s1=c1.length();
		s2=c2.length();

		out1.clear();
		out2.clear();

		out.resize(s1+1);

		score=new int [(s1+1)*(s2+1)];
		point=new int [(s1+1)*(s2+1)];

		for(i=0;i<=s1;i++)
		for(j=0;j<=s2;j++)
		{
			score[j*(s1+1)+i]=0;
			point[j*(s1+1)+i]=0;
		}


		for(i=1;i<=s1;i++)
		{
			score[i]=-i;
			point[i]=3;
		}
		for(j=1;j<=s2;j++)
		{
			score[j*(s1+1)]=-j;	
			point[j*(s1+1)]=2;
		}

		for(i=1;i<=s1;i++)
		{
			for(j=1;j<=s2;j++)
			{
				if(c1[i-1] == c2[j-1] || c1[i-1]=='X' || c2[j-1]=='X') // 'X' can match any AA
					score1=score[(j-1)*(s1+1)+i-1]+1;
				else
					score1=score[(j-1)*(s1+1)+i-1]-1;
				score2=score[(j-1)*(s1+1)+i]-1;
				score3=score[j*(s1+1)+i-1]-1;

				if(score1>score2)
				{
					if(score1>score3)
					{
						score[j*(s1+1)+i]=score1;
						point[j*(s1+1)+i]=1;
					}
					else
					{
						score[j*(s1+1)+i]=score3;
						point[j*(s1+1)+i]=3;
					}
				}
				else
				{
					if(score2>score3)
					{
						score[j*(s1+1)+i]=score2;
						point[j*(s1+1)+i]=2;
					}
					else
					{
						score[j*(s1+1)+i]=score3;
						point[j*(s1+1)+i]=3;
					}
				}
			}
		}



		i=s1;
		j=s2;
		while(point[j*(s1+1)+i]!=0)
		{
			if(point[j*(s1+1)+i]==1)
			{
				out1.append(1,c1[i-1]);
				out2.append(1,c2[j-1]);
				if(c1[i-1]==c2[j-1] || c1[i-1]=='X' || c2[j-1]=='X' ) 
				{
					out[i]=j;
					out3.append(1,c1[i-1]);
				}
				else
					out3.append(1,'-');
				i--;
				j--;
			}
			else if(point[j*(s1+1)+i]==2)
			{
				out1.append(1,'-');
				out2.append(1,c2[j-1]);
				out3.append(1,' ');
				j--;
			}
			else if(point[j*(s1+1)+i]==3)
			{
				out1.append(1,c1[i-1]);
				out2.append(1,'-');
				out3.append(1,' ');
				i--;
			}
		}

		delete [] score;
		delete [] point;
		
		return out;
	}

}
