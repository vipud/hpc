
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>
#include <time.h>
#include <sstream>
using namespace std;

#include "config.h"
#include "supply.h"
#include "bmrb.h"
#include  "aa.h"
#include "pdb.h"
//#include "munkres.h"




vector<double> compare(vector< vector<double> > t,vector<double> shift)
{
	int i,j,n;
	double e,ee;
	bool b;
	vector<bool> bs;
	double sum;
	int count;
	vector<double> diff;

	for(i=0;i<t.at(0).size();i++)
	{
		b=0;
		for(j=0;j<t.size();j++)
		{
			if(t.at(j).at(i)<=0.001 || t.at(j).at(i)>990)
				b=1;
		}
		bs.push_back(b);
		
	}

	cout<<"RMSDs are";
	sum=0.0;
	count=0;
	for(i=1;i<t.size();i++)
	{
		ee=0;
		n=0;
		for(j=0;j<t.at(0).size();j++)
		{
			if(bs.at(j)==0)
			{
				e=t.at(0).at(j)-t.at(i).at(j)-shift.at(i-1);
				ee+=e*e;
				n++;
				sum+=t.at(i).at(j)-t.at(0).at(j);
				count++;
			}
		}
		ee/=n;
		ee=sqrt(ee);
		cout<<" "<<ee;
		if(count>0) diff.push_back(sum/count);
		else diff.push_back(0.0);
	}
	cout<<endl;

	return diff;
}

void output(char *name,vector<int> ids, vector<char> codes, vector<string> names, vector< vector<double> > t)
{
	int i,j;
	FILE * fp;
	fp=fopen(name,"wt");

	for(i=0;i<ids.size();i++)
	{
		fprintf(fp,"%8d%8c%8s",ids.at(i),codes.at(i),names.at(i).c_str());
		for(j=0;j<t.size();j++)
			fprintf(fp,"%8.3f",t.at(j).at(i));
		fprintf(fp,"\n");
	}
	fclose(fp);
	return;
}
		
void output2(char *name,vector<int> ids, vector<char> codes, vector<string> names, vector<int> types, vector< vector<double> > t)
{
	int i,j;
	FILE * fp;
	fp=fopen(name,"wt");

	for(i=0;i<ids.size();i++)
	{
		fprintf(fp,"%8d%8d%8c%8s",ids.at(i),types.at(i),codes.at(i),names.at(i).c_str());
		for(j=0;j<t.size();j++)
			fprintf(fp,"%10.3f",t.at(j).at(i));
		fprintf(fp,"\n");
	}
	fclose(fp);
	return;
}


double process(int &count,double x,double y,double scale)
{
	double s;

	s=(x-y)*scale;
	s=fabs(s);

	if(x>900 || y>900) 
	{
		s=0;
		count--;
	}
	else if(x<-900 && y<-900)
	{
		s=-5.0;
	}
	else if(x<-900 || y<-900)
	{
		s=100;
	}
	else
		s*=s;

	
	return s;
};

double adj(double x, double y, double s)
{
	y=(y-x)*s+x;
	return y;
}

			


/*void bmrb_matrix(double scale,vector<int> &correct_assignment,vector<int> &false_assignment,vector<struct bb_group> bb1,vector<struct bb_group> bb2,int selection, int nnoise, double std, double cutoff, double wc, double wh, double wn, double wha)
{
	int n,count;
	int total;
	int i,j;
	struct bb_group b1,b2;
	double s,ss;
	Munkres *m;
	bool bfalse;
	int c,f;

	ofstream fout("cs_matrix.dat");
	ofstream fassign("assign.dat");

	n=bb1.size();
	if(n!=bb2.size())
		return;

	for(i=0;i<bb2.size();i++)
	{
		bb2.at(i).exp_ca=adj(bb1.at(i).exp_ca,bb2.at(i).exp_ca,scale);
		bb2.at(i).exp_cb=adj(bb1.at(i).exp_cb,bb2.at(i).exp_cb,scale);
		bb2.at(i).exp_co=adj(bb1.at(i).exp_co,bb2.at(i).exp_co,scale);
		bb2.at(i).follow_exp_n=adj(bb1.at(i).follow_exp_n,bb2.at(i).follow_exp_n,scale);
		bb2.at(i).follow_exp_h=adj(bb1.at(i).follow_exp_h,bb2.at(i).follow_exp_h,scale);
	}

	Matrix<double> matrix(n, n);  //nrow, ncolumn
	Matrix<double> matrix2(n, n);  //nrow, ncolumn
	Matrix<double> outmatrix(n, n);  //nrow, ncolumn
	
	for(i=0;i<n;i++)
	{
		b1=bb1.at(i);
		for(j=0;j<n;j++)
		{
			b2=bb2.at(j);
			ss=0;


			if(selection==1) //ca cb ha
			{
				count=3;
				ss+=process(count,b1.exp_ca,b2.exp_ca,wc);
				ss+=process(count,b1.exp_cb,b2.exp_cb,wc);
				ss+=process(count,b1.exp_ha,b2.exp_ha,wha);
			}
			else if(selection==2) //co h n
			{
				count=3;
				ss+=process(count,b1.exp_ca,b2.exp_ca,wc);
				ss+=process(count,b1.follow_exp_h,b2.follow_exp_h,wh);
				ss+=process(count,b1.follow_exp_n,b2.follow_exp_n,wn);
			}
			else if(selection==3) //ca cb h n
			{
				count=4;
				ss+=process(count,b1.exp_ca,b2.exp_ca,wc);
				ss+=process(count,b1.exp_cb,b2.exp_cb,wc);
				ss+=process(count,b1.follow_exp_h,b2.follow_exp_h,wh);
				ss+=process(count,b1.follow_exp_n,b2.follow_exp_n,wn);
			}
			else if(selection==4) //ca cb h n ha
			{
				count=5;
				ss+=process(count,b1.exp_ca,b2.exp_ca,wc);
				ss+=process(count,b1.exp_cb,b2.exp_cb,wc);
				ss+=process(count,b1.follow_exp_h,b2.follow_exp_h,wh);
				ss+=process(count,b1.follow_exp_n,b2.follow_exp_n,wn);
				ss+=process(count,b1.exp_ha,b2.exp_ha,wha);
			}
			else  // all 
			{
				count=6;
				ss+=process(count,b1.exp_ca,b2.exp_ca,wc);
				ss+=process(count,b1.exp_cb,b2.exp_cb,wc);
				ss+=process(count,b1.exp_co,b2.exp_co,wc);
				ss+=process(count,b1.follow_exp_h,b2.follow_exp_h,wh);
				ss+=process(count,b1.follow_exp_n,b2.follow_exp_n,wn);
				ss+=process(count,b1.exp_ha,b2.exp_ha,wha);
			}

			if(ss<0)
				ss=0.0;

			if(count>0)
				ss=sqrt(ss/count);
			else
				ss=10.0;

			fout<<ss<<" ";
			matrix(i,j)=ss;
		}
		fout<<endl;
	}
	fout<<endl<<endl;
	

	
	for ( int row = 0 ; row < n ; row++ )
	for ( int col = 0 ; col < n ; col++ )
	{
		outmatrix(row,col)=0.0;
	}

	for(i=0;i<nnoise;i++)
	{
		if(i>0)
		{
			for ( int row = 0 ; row < n ; row++ )
			{
				for ( int col = 0 ; col < n ; col++ )
				{
					matrix2(row,col)=matrix(row,col)+ldw_math::gaussrand()*std;
				}
			}
		}
		else
		{
			for ( int row = 0 ; row < n ; row++ )
			{
				for ( int col = 0 ; col < n ; col++ )
				{
					matrix2(row,col)=matrix(row,col);
				}
			}
		}

		m=new Munkres;
		m->solve(matrix2);
		delete m;

		for ( int row = 0 ; row < n ; row++ )
		{
			for ( int col = 0 ; col < n ; col++ )
			{
				fout<<matrix2(row,col)+1<<" ";
			}
			fout<<endl;
		}
		fout<<endl<<endl;

		for ( int row = 0 ; row < n ; row++ )
		{
			for ( int col = 0 ; col < n ; col++ )
			{
				outmatrix(row,col)=outmatrix(row,col)+matrix2(row,col)+1;
			}			
		}
	}

	f=c=0;
	for ( int row = 0 ; row < n ; row++ )
	{
		bfalse=0;
		for ( int col = 0 ; col < n ; col++ )
		{
			fassign << outmatrix(row,col) << " ";
			if(col!=row && outmatrix(row,col)>=cutoff*nnoise)
			{
				bfalse=1;
				false_assignment.push_back(bb1.at(row).id);
				false_assignment.push_back(bb2.at(col).id);
			}
		}
		fassign << endl;
		if(bfalse)
			f++;
		else if(outmatrix(row,row)>=cutoff*nnoise)
		{
			c++;
			correct_assignment.push_back(bb1.at(row).id);
		}
	}
		
	fassign.close();
	fout.close();
	cout<<"correctly assign "<<c<<" while false positive is "<<f<<" and total is "<<n<<endl;
	return;
}
*/

double cal_diff(double &count,double x,double y,double scale)
{
	double s;

	if(fabs(x)>900 || fabs(y)>900) 
	{
		s=0.0;
		count-=scale;
	}
	else if(fabs(x)<=0.01 || fabs(y)<=0.01)
	{
		s=0.0;
		count-=scale;
	}
	else
	{
		s=(x-y)*scale;
		s*=s;
	}

	return s;
};
			

void compare_cs(vector<struct bb_group> bb1,vector<struct bb_group> bb2, vector<double> &errors)
{
	int i;
	struct bb_group b1,b2;
	double error;
	double count;
	
	for(i=0;i<(int)bb1.size() && i<(int)bb2.size() ;i++)
	{
		b1=bb1.at(i);
		b2=bb2.at(i);

		count=3.0;
		error=0.0;
		error+=cal_diff(count,b1.exp_ca,b2.exp_ca,1.0);
		error+=cal_diff(count,b1.exp_cb,b2.exp_cb,1.0);
		error+=cal_diff(count,b1.exp_co,b2.exp_co,1.0);
		//error+=cal_diff(count,b1.exp_h,b2.exp_h,2.0);
		//error+=cal_diff(count,b1.exp_n,b2.exp_n,0.4);
		
		if(count>=1.0)
			error=error/count;
		else
			error=0.0;
		errors.push_back(sqrt(error));
	}
	return;
}


	
void computer_rmsd(vector < vector<struct bb_group> > *bbs, string filename)
{
	int i,j;
	vector<struct bb_group> bb1,bb2;
	vector<double> errors;
	vector<int> ids;
	vector< vector<double> > all_errors;
	ofstream fout(filename.c_str());

	bb1=bbs->at(0);
	for(i=0;i<(int)bb1.size();i++)
		ids.push_back(bb1.at(i).id);
	for(i=1;i<(int)bbs->size();i++)
	{
		bb2.clear();
		errors.clear();
		bb2=bbs->at(i);
		compare_cs(bb1,bb2,errors);
		all_errors.push_back(errors);
	}

	for(i=0;i<(int)ids.size();i++)
	{
		fout<<ids.at(i);
		for(j=0;j<all_errors.size();j++)
			fout<<" "<<all_errors.at(j).at(i);
		fout<<endl;
	}
}


int main(int argc,char ** argv)
{
	int i,j;
	class CPdb pdb;
	class CBmrb bmrb,bmrb0;
	string filename;
	vector<struct bb_group> bb;
	vector<struct bb_group> bb2;
	vector < vector<struct bb_group> > bbs;
	vector<struct proton> proton;
	vector< vector<double> > cas,cbs,cos,hs,ns,hydrogens,has;
	vector<double> ca,cb,co,h,n,hydrogen,ha;
	vector<int> ids,ids2;
	vector<char> codes,codes2;
	vector<int> types;
	vector<string> names,names2;
	vector<string> args,args2;
	CCommandline cmdline;
	stringstream ss;
	string p;
	string bmrbname;
	string pre_name;
	bool bfirst;
	vector<int> false_assignment;
	vector<int> correct_assignment;
	string seq;



	vector<double> shift,shift1,shift2,shift3;
	

	bbs.clear();

	cout<<"NMR suit Version "<<NMR_VERSION_MAJOR<<"."<<NMR_VERSION_MINOR<<endl;

	
	args.push_back("-pdb");args2.push_back("pdb.pdb");
	args.push_back("--bmrb");args2.push_back("bmrb.dat bmrb_pre.dat");
	args.push_back("--sparta");args2.push_back("");
	args.push_back("-bb");args2.push_back("yes");
	args.push_back("-proton");args2.push_back("yes");
	args.push_back("-rmsd");args2.push_back("yes");


	cmdline.init(args,args2);
	cmdline.pharse(argc,argv);
	cmdline.print();


	FILE *fp=fopen(cmdline.query("-pdb").c_str(),"rt");
	if(fp==NULL)
	{
		bmrbname=cmdline.query("--bmrb");
		ss.clear();
		ss<<bmrbname;
		if(ss>>p)
		{
			bmrb0.clear();
			bmrb0.process(p.c_str());
			bmrb0.loaddetail("cs.dat");
			while(ss>>p)
			{
				bmrb.clear();
				bmrb.process(p.c_str());
				bmrb0.attach_bmrb(bmrb);
			}
		}

		bmrb0.print("detail.dat");

	}
	else
	{
		pdb.loadpdb(cmdline.query("-pdb").c_str());

		bfirst=1;

		bmrbname=cmdline.query("--bmrb");
		ss.clear();
		ss<<bmrbname;
		while(ss>>p)
		{
			pdb.loadpdb(cmdline.query("-pdb").c_str());
			bmrb.clear();
			bmrb.process(p.c_str());
			pdb.attach_bmrb(bmrb);
			bb.clear();
			pdb.getbb(&bb);
			bbs.push_back(bb);

			ca.clear();cb.clear();co.clear();h.clear();n.clear();ha.clear();


			for(j=0;j<bb.size();j++)
			{
				if(bfirst==1)
				{
					ids.push_back(bb.at(j).id);
					codes.push_back(bb.at(j).code);
					names.push_back(Sequence::code2name(bb.at(j).code));
				}
				ca.push_back(bb.at(j).exp_ca);
				cb.push_back(bb.at(j).exp_cb);
				co.push_back(bb.at(j).exp_co);
				 h.push_back(bb.at(j).exp_h);
				 n.push_back(bb.at(j).exp_n);
				ha.push_back(bb.at(j).exp_ha);
			}
			cas.push_back(ca);
			cbs.push_back(cb);
			cos.push_back(co);
			has.push_back(ha);
			hs.push_back(h);
			ns.push_back(n);

			proton.clear();
			pdb.allproton(&proton);
			hydrogen.clear();
			for(j=0;j<proton.size();j++)
			{
				if(bfirst==1)
				{
					ids2.push_back(proton.at(j).id);
					codes2.push_back(proton.at(j).code);
					names2.push_back(proton.at(j).name);
					types.push_back(proton.at(j).type);
				}
				hydrogen.push_back(proton.at(j).exp);
			}
			hydrogens.push_back(hydrogen);

			bfirst=0;
		}


		bmrbname=cmdline.query("--sparta");
		ss.clear();
		ss<<bmrbname;
		while(ss>>p)
		{
			pdb.loadpdb(cmdline.query("-pdb").c_str());
			bmrb.clear();
			bmrb.run_sparta(p,0);
			pdb.attach_bmrb(bmrb);
			pre_name=p.append("_bmrb.dat");
			pdb.print_prediction(pre_name);
		
			bb.clear();
			pdb.getbb(&bb);
			bbs.push_back(bb);

			ca.clear();cb.clear();co.clear();h.clear();n.clear();ha.clear();
			for(j=0;j<bb.size();j++)
			{
				if(bfirst==1)
				{
					ids.push_back(bb.at(j).id);
					codes.push_back(bb.at(j).code);
					names.push_back(Sequence::code2name(bb.at(j).code));
				}
				ca.push_back(bb.at(j).exp_ca);
				cb.push_back(bb.at(j).exp_cb);
				co.push_back(bb.at(j).exp_co);
				 h.push_back(bb.at(j).exp_h);
				 n.push_back(bb.at(j).exp_n);
				ha.push_back(bb.at(j).exp_ha);
			}
			cas.push_back(ca);
			cbs.push_back(cb);
			cos.push_back(co);
			has.push_back(ha);
			hs.push_back(h);
			ns.push_back(n);

			proton.clear();
			pdb.allproton(&proton);
			hydrogen.clear();
			for(j=0;j<proton.size();j++)
			{
				if(bfirst==1)
				{
					ids2.push_back(proton.at(j).id);
					codes2.push_back(proton.at(j).code);
					names2.push_back(proton.at(j).name);
					types.push_back(proton.at(j).type);
				}
				hydrogen.push_back(proton.at(j).exp);
			}
			hydrogens.push_back(hydrogen);

			bfirst=0;
		}

	

	
	if(cmdline.query("-bb").compare("yes")==0)
	{
		for(i=1;i<cas.size();i++)
			shift.push_back(0.0);

		cout<<"Ca ";shift1=compare(cas,shift);
		cout<<"Cb ";shift2=compare(cbs,shift);
		cout<<"Co ";shift3=compare(cos,shift);




		output("ca.out",ids,codes,names,cas);output("cb.out",ids,codes,names,cbs);output("co.out",ids,codes,names,cos);
		output("h.out",ids,codes,names,hs);output("n.out",ids,codes,names,ns);
		output("ha.out",ids,codes,names,has);

		shift.clear();
		for(i=1;i<cas.size();i++)
		{
			shift.push_back(-(shift1.at(i-1)+shift2.at(i-1))/2.0);
            cout<<"autoref of cacb is "<<-(shift1.at(i-1)+shift2.at(i-1))/2.0<<endl;
		}
        

		cout<<"Ca ";compare(cas,shift);
		cout<<"Cb ";compare(cbs,shift);

		shift.clear();
		for(i=1;i<cas.size();i++)
		{
			shift.push_back(-shift3.at(i-1));
            cout<<"autoref of co is "<<-shift3.at(i-1)<<endl;
		}

		cout<<"Co ";compare(cos,shift);

		shift.clear();
		for(i=1;i<cas.size();i++)
			shift.push_back(0.0);

		cout<<"ha ";compare(has,shift);
		cout<<"h  ";compare(hs,shift);
		cout<<"n  ";compare(ns,shift);
	}

	if(cmdline.query("-proton").compare("yes")==0)
	{
		shift.clear();
		for(i=0;i<hydrogens.size();i++)
			shift.push_back(0.0);
		output2("proton.out",ids2,codes2,names2,types,hydrogens);
		cout<<"Side chain proton ";
		compare(hydrogens,shift);
	}

	}



	return 0;
}



		


