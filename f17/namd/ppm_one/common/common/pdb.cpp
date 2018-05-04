#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <omp.h>


using namespace std;


#include "pdb.h"
#include "supply.h"
#include "debug.h"


int CDssp::loaddata(string name)
{
	char buffer[255];
	bool begin;
	struct dsspline aadssp;
	string line;

	sprintf(buffer,"dssp %s > dssp.out",name.c_str());
	system(buffer); 

	ifstream fin("dssp.out");
	if (!fin.is_open())
		return 1;

	
	begin=0;
	while(getline(fin,line))
	{
		if(line.find("RESIDUE AA")!=std::string::npos)
		{
			begin=1;
			continue;
		}
		if(begin==1)
		{
			aadssp.id=atoi(line.substr(5,5).c_str());
			aadssp.code=line[13];
			aadssp.ss=line[16];
			data.push_back(aadssp);
		}
	}
	return 0;
}

string CDssp::getseq()
{
	int i;
	string s;

	for(i=0;i<(int)data.size();i++)
		s.push_back(data.at(i).code);
	return s;
}

CDssp::CDssp()
{}
CDssp::~CDssp()
{}

//////////////////////////
// class CPdb           //
//////////////////////////

void CPdb::process_ambig(int flag)
{
	double st = omp_get_wtime();
	//meaning of flag:
	// 1. remove any ambig assignment
	// 2. combine all hb2, hb3 into mean value, if code is 1 or 2, otherwise remove ambig
	int i;
	if(flag==1)
	{
		for(i=0;i<(int)v.size();i++)
		{
			v.at(i)->remove_ambig(flag);
		}
	}
	
	//not really needed ??? if take mean in function get all proton
	if(flag==2)
	{
		for(i=0;i<(int)v.size();i++)
		{
			v.at(i)->methyl_ambig(flag);
			//v.at(i)->combine_hsamec(flag);
		}
	}

	if(flag==3)
	{
		//remove methyl ambig==2 case
	}
	cout << "pdb::process_ambig: " << omp_get_wtime()-st << " seconds" << endl;
}


void CPdb::clear()
{
	unsigned int i;
	for(i=0;i<(int)v.size();i++)
	{
		delete v.at(i);
	}
	v.clear();
	for(i=0;i<(int)ligand.size();i++)
	{
		delete ligand.at(i);
	}
	ligand.clear();
	chains.clear();
	chain_block.clear();
	chain_ligand.clear();
	blocks.clear();
	pdbseq.clear();
}

void CPdb::print_print(FILE *fp,vector<int> att1, vector<int> att2,int index,int type,double a,double b,double c,double w)
{
	int i,j;

	for(i=0;i<att1.size();i++)
	for(j=0;j<att2.size();j++)
	{	
		fprintf(fp,"%10d%10d%10d%10d%10d",att1.at(i),att2.at(j),1,index,type);
		fprintf(fp,"%10.4f%10.4f%10.4f%10.4f\n",a/10,b/10,c/10,w);  //unit is nm instead of A in Gromacs!
	}
	return;
}

void CPdb::print_gromacs(string filename,vector <struct diheline> *dihecons)
{
	FILE *fp;
	struct noeline tline;
	int i,j;

	fp=fopen(filename.c_str(),"wt");
	fprintf(fp,"[dihedral_restraints]\n");
	for(i=0;i<dihecons->size();i++)
	{
		for(j=0;j<4;j++) fprintf(fp,"%d ",dihecons->at(i).index[j].atoms.at(0).at(0));
		fprintf(fp," 1 1 ");
		fprintf(fp,"%8.3f%8.3f",(float)dihecons->at(i).middle,(float)dihecons->at(i).delta);
		fprintf(fp," 1 2\n");
	}
	fclose(fp);
	return;
}



void CPdb::print_gromacs(string filename,vector <struct noeline> *nmrcons)
{
	FILE *fp;
	struct noeline tline;
	int i,j,index;
	int n1,n2;
	double length1,length2;
	vector<int> att1,att2;
	string filename2;

	vector<int> nshort;
	vector<int> nlong;
	vector< vector<int> > whos;
	vector< vector<int> > pairs;
	vector< vector<int> > used;
	vector< vector<double> > weights;

	nshort.resize(v.size());
	nlong.resize(v.size());
	whos.resize(v.size());
	pairs.resize(v.size());
	weights.resize(v.size());
	used.resize(v.size());
	for(i=0;i<v.size();i++)
		used.at(i).resize(v.size(),0);


	for(i=0;i<pairs.size();i++)
	{
		pairs.at(i).resize(v.size());
		weights.at(i).resize(v.size());
	}

	for(i=0;i<nmrcons->size();i++)
	{
		if(fabs((float)(nmrcons->at(i).resid1-nmrcons->at(i).resid2))<=4)
		{
			nshort.at(nmrcons->at(i).resid1-1)++;
			nshort.at(nmrcons->at(i).resid2-1)++;
		}
		else
		{
			nlong.at(nmrcons->at(i).resid1-1)++;
			nlong.at(nmrcons->at(i).resid2-1)++;
			whos.at(nmrcons->at(i).resid1-1).push_back(nmrcons->at(i).resid2);
			whos.at(nmrcons->at(i).resid2-1).push_back(nmrcons->at(i).resid1);
			pairs.at(nmrcons->at(i).resid1-1).at(nmrcons->at(i).resid2-1)++;
			pairs.at(nmrcons->at(i).resid2-1).at(nmrcons->at(i).resid1-1)++;
		}

	}

	fp=fopen("cons_detail.dat","wt");
	for(i=0;i<v.size();i++)
	{
		fprintf(fp,"%d %d %d\n",i+1,nshort.at(i),nlong.at(i));
	}
	fclose(fp);

	fp=fopen("cons_whos.dat","wt");
	for(i=0;i<whos.size();i++)
	{
		fprintf(fp,"%d",i+1);
		for(j=0;j<whos.at(i).size();j++)
			fprintf(fp," %d",whos.at(i).at(j));
		fprintf(fp,"\n");
	}
	fclose(fp);


	int b1,b2,s1,s2,k1,k2;
	double c;
	for(i=0;i<v.size();i++)
	{
		for(j=0;j<v.size();j++)
		{
			weights.at(i).at(j)=1;
			if(pairs.at(i).at(j)>0)
				weights.at(i).at(j)=1.0/pairs.at(i).at(j);
			c=0;
			b1=-1; if(i==0) b1=0;
			b2=-1; if(j==0) b2=0;
			s1=1; if(i==v.size()-1) s1=0;
			s2=1; if(j==v.size()-1) s2=0;
			for(k1=i+b1;k1<=i+s1;k1++)
			for(k2=j+b2;k2<=j+s2;k2++)
				c+=pairs.at(k1).at(k2);
			c-=pairs.at(i).at(j);
			if(c>0) weights.at(i).at(j)/=sqrt(c);
		}
	}


	index=-1;
	fp=fopen(filename.c_str(),"wt");
	fprintf(fp,"[distance_restraints]\n");
	fprintf(fp,"; all long range noe distance restraints\n");
	for(i=0;i<nmrcons->size();i++)
	{
		if(fabs((float)(nmrcons->at(i).resid1-nmrcons->at(i).resid2))<=4)
			continue;
		tline=nmrcons->at(i);
		length1=0;
		length2=0;

		for(n1=0;n1<(int)nmrcons->at(i).index1.atoms.size();n1++)
		{
			length1=nmrcons->at(i).index1.length;
			for(n2=0;n2<(int)nmrcons->at(i).index2.atoms.size();n2++)
			{
				index++;
				length2=nmrcons->at(i).index1.length;
				att1=nmrcons->at(i).index1.atoms.at(n1);
				att2=nmrcons->at(i).index2.atoms.at(n2);
				print_print(fp,att1,att2,index,1,0,nmrcons->at(i).a+length1+length2,nmrcons->at(i).a+2+length1+length2,weights.at(tline.resid1-1).at(tline.resid2-1));
			}
		}
	}
	fclose(fp);

	fp=fopen("cns_distance.tbl","wt");
	for(i=0;i<nmrcons->size();i++)
	{
		if(nmrcons->at(i).bvio==1)
			continue;
		if(fabs((float)(nmrcons->at(i).resid1-nmrcons->at(i).resid2))<=4)
		{
			tline=nmrcons->at(i);
			fprintf(fp,"assign (resid %d and name %s ) (resid %d and name %s ) %f %f %f\n",
				tline.oldresid1,tline.oldatomname1.c_str(),tline.oldresid2,tline.oldatomname2.c_str(),
				(float)tline.b,(float)tline.c,(float)tline.a);
		}
		else
		{
			tline=nmrcons->at(i);
			if(used.at(tline.resid1-1).at(tline.resid2-1)==0)
			{
				fprintf(fp,"assign (resid %d and name %s ) (resid %d and name %s ) %f %f %f\n",
					tline.oldresid1,tline.oldatomname1.c_str(),tline.oldresid2,tline.oldatomname2.c_str(),
					(float)tline.b,(float)tline.c,(float)tline.a);
				used.at(tline.resid1-1).at(tline.resid2-1)++;
				used.at(tline.resid2-1).at(tline.resid1-1)++;
			}
		}
	}
	fclose(fp);


	index=-1;
	filename2=filename;
	filename2.insert(filename2.find("."),"2");
	fp=fopen(filename2.c_str(),"wt");
	fprintf(fp,"[distance_restraints]\n");
	fprintf(fp,"; vilated noe distance restraints\n");
	for(i=0;i<nmrcons->size();i++)
	{
		if(nmrcons->at(i).bvio==0)
			continue;
		tline=nmrcons->at(i);
		length1=0;
		length2=0;

		for(n1=0;n1<(int)nmrcons->at(i).index1.atoms.size();n1++)
		{
			length1=nmrcons->at(i).index1.length;
			for(n2=0;n2<(int)nmrcons->at(i).index2.atoms.size();n2++)
			{
				index++;
				length2=nmrcons->at(i).index1.length;
				att1=nmrcons->at(i).index1.atoms.at(n1);
				att2=nmrcons->at(i).index2.atoms.at(n2);
				print_print(fp,att1,att2,index,1,0,nmrcons->at(i).a+length1+length2,nmrcons->at(i).a+2+length1+length2,1.0);
			}
		}
	}
	fclose(fp);


	index=-1;
	filename2=filename;
	filename2.insert(filename2.find("."),"3");
	fp=fopen(filename2.c_str(),"wt");
	fprintf(fp,"[distance_restraints]\n");
	fprintf(fp,"; fullfiled noe distance restraints\n");
	for(i=0;i<nmrcons->size();i++)
	{
		if(nmrcons->at(i).bvio==1)
			continue;
		tline=nmrcons->at(i);
		length1=0;
		length2=0;

		for(n1=0;n1<(int)nmrcons->at(i).index1.atoms.size();n1++)
		{
			length1=nmrcons->at(i).index1.length;
			for(n2=0;n2<(int)nmrcons->at(i).index2.atoms.size();n2++)
			{
				index++;
				length2=nmrcons->at(i).index1.length;
				att1=nmrcons->at(i).index1.atoms.at(n1);
				att2=nmrcons->at(i).index2.atoms.at(n2);
				print_print(fp,att1,att2,index,1,0,nmrcons->at(i).a+length1+length2,nmrcons->at(i).a+2+length1+length2,1.0);
			}
		}
	}
	fclose(fp);


	index=-1;
	filename2=filename;
	filename2.insert(filename2.find("."),"4");
	fp=fopen(filename2.c_str(),"wt");
	fprintf(fp,"[distance_restraints]\n");
	fprintf(fp,"; all noe distance restraints\n");
	for(i=0;i<nmrcons->size();i++)
	{
		tline=nmrcons->at(i);
		length1=0;
		length2=0;

		for(n1=0;n1<(int)nmrcons->at(i).index1.atoms.size();n1++)
		{
			length1=nmrcons->at(i).index1.length;
			for(n2=0;n2<(int)nmrcons->at(i).index2.atoms.size();n2++)
			{
				index++;
				length2=nmrcons->at(i).index1.length;
				att1=nmrcons->at(i).index1.atoms.at(n1);
				att2=nmrcons->at(i).index2.atoms.at(n2);
				print_print(fp,att1,att2,index,1,0,nmrcons->at(i).a+length1+length2,nmrcons->at(i).a+2+length1+length2,1.0);
			}
		}
	}
	fclose(fp);

	return;
}


int CPdb::loadnoedata(string filename,string filename2,vector <struct noeline> *nmrcons)
{
	string nmrseq;
	int adj,adj2;
	vector<int> out;
	char c1,c2;
	ifstream fin(filename.c_str());
	ifstream fin2(filename2.c_str());
	istringstream iss;
	string line,p;
	string resname1,atomname1,resname2,atomname2,resname,atomname;
	int resid1,resid2,order,resid;
	double a,b,c;
	int i;
	struct noeline tline;
	bool bfirst;
	int index,index_old;
	int norder1,norder2;

	i=0;
	bfirst=1;
	while(getline(fin,line))
	{
		i++;
		iss.clear();
		iss.str(line);
		iss>>index;
		iss>>p>>order>>p>>p;
		iss>>resid;
		iss>>resname;
		iss>>atomname;
		iss>>p>>p>>p>>p>>p>>p;

		if(bfirst)
		{
			bfirst=0;
			resid1=resid2=-1;
			index_old=index;
			norder1=norder2=0;
		}

		if(index>index_old) //we now reach a new entry
		{
			tline.id=index_old;
			tline.resid1=resid1;
			tline.resid2=resid2;
			tline.resname1=resname1;
			tline.resname2=resname2;
			tline.atomname1=atomname1;
			tline.atomname2=atomname2;
			tline.multi=0;
			if(norder1>1 || norder2>1)
				tline.multi=1;
			nmrcons->push_back(tline);

			resid1=resid2=-1;
			index_old=index;
			norder1=norder2=0;
		}

		if(order==1)
		{
			norder1++;
			resid1=resid;
			resname1=resname;
			atomname1=atomname;
		}
		else if(order==2)
		{
			norder2++;
			resid2=resid;
			resname2=resname;
			atomname2=atomname;
		}
		else
		{
			cout<<"Unsupported order "<<order<<endl;
			break;
		}
	}

	//process the final entry.
	if(resid1>0 && resid2>0)
	{
		tline.id=index_old;
		tline.resid1=resid1;
		tline.resid2=resid2;
		tline.resname1=resname1;
		tline.resname2=resname2;
		tline.atomname1=atomname1;
		tline.atomname2=atomname2;
		nmrcons->push_back(tline);
	}


	i=0;
	while(getline(fin2,line))
	{
		iss.clear();
		iss.str(line);
		iss>>index;
		if(index!=nmrcons->at(i).id)
		{
			cout<<"readin noe data error. index is NOT as expected in second part"<<endl;
			break;
		}
		iss>>p>>p>>p>>p>>p>>p;
		iss>>b>>c>>a;
		iss>>p>>p;
		nmrcons->at(i).a=a;
		nmrcons->at(i).b=b;
		nmrcons->at(i).c=c;


		c1=Sequence::name2code(nmrcons->at(i).resname1);
		c2=Sequence::name2code(nmrcons->at(i).resname2);
		nmrcons->at(i).resid1--;
		nmrcons->at(i).resid2--;
		if(nmrcons->at(i).resid1+1>nmrseq.size())
			nmrseq.resize(nmrcons->at(i).resid1+1,'U');
		nmrseq.at(nmrcons->at(i).resid1)=c1;
		if(nmrcons->at(i).resid2+1>nmrseq.size())
			nmrseq.resize(nmrcons->at(i).resid2+1,'U');
		nmrseq.at(nmrcons->at(i).resid2)=c2;

		i++;
	}

	out=Sequence::align(pdbseq,nmrseq);
	adj=0;adj2=0;
	for(i=0;i<out.size();i++)
	{
		if(out.at(i)!=0)
		{
			adj+=(out.at(i)-i);
			adj2++;
		}
	}
	adj/=adj2;


	//remove entry with multiply atoms
	for(i=nmrcons->size()-1;i>=0;i--)
	{
		nmrcons->at(i).resid1++;
		nmrcons->at(i).resid2++;
		if(nmrcons->at(i).multi==1)
		{
			cout<<"Remove noe entry with umbigirious assignment. "<<nmrcons->at(i).id<<endl;
			nmrcons->erase(nmrcons->begin()+i);
		}
	}


	for(i=0;i<nmrcons->size();i++)
	{
		nmrcons->at(i).index1.length=0.0;
		nmrcons->at(i).index2.length=0.0;
		nmrcons->at(i).resid1-=adj;
		nmrcons->at(i).resid2-=adj;

		if(nmrcons->at(i).resid1>=1 && nmrcons->at(i).resid2>=1 && nmrcons->at(i).resid1<=(int)v.size() && nmrcons->at(i).resid2<=(int)v.size() )
		{
			nmrcons->at(i).index1=v.at(nmrcons->at(i).resid1-1)->query(nmrcons->at(i).atomname1);
			nmrcons->at(i).index2=v.at(nmrcons->at(i).resid2-1)->query(nmrcons->at(i).atomname2);
		}
		else
		{
			nmrcons->erase(nmrcons->begin()+i);
			i--;
		}
	}
	return adj;

}




int CPdb::loadnoedata(string filename,vector <struct noeline> *nmrcons)
{

	ifstream fin(filename.c_str());
	string line;
	vector<string> block1, block2;
	bool bstart,inloop,begin1,begin2;

	bstart=0;
	inloop=0;
	begin1=0;
	begin2=0;

	while(getline(fin,line))
	{
		if(line.size()<4)
			continue;

		if(begin2==1 &&  line.find("stop_")!=string::npos)
			bstart=begin2=0;
		if(begin1==1 &&  line.find("stop_")!=string::npos)
			begin1=0;


		if(begin1==1)
			block1.push_back(line);
		if(begin2==1)
			block2.push_back(line);

		if(line.find("loop_")!=string::npos)
			inloop=1;
		if(line.find("stop_")!=string::npos)
			inloop=0;

		if(inloop==0 && ( line.find("save_CNS/XPLOR_distance_constraints")!=string::npos || line.find("save_DYANA/DIANA_distance_constraints")!=string::npos) )
		{
			getline(fin,line);
			getline(fin,line);
			getline(fin,line);
			getline(fin,line);
			
			if(line.find("_Distance_constraint_list.Constraint_type")!=string::npos)
			{
				if( line.find("NOE")!=string::npos )
					bstart=1;
			}
			getline(fin,line);
			getline(fin,line);
		}

		if(bstart==1 && inloop==1 && line.find("_Dist_constraint.Distance_constraint_list_ID")!=string::npos)
			begin1=1;
		if(bstart==1 && inloop==1 && line.find("_Dist_constraint_value.Distance_constraint_list_ID")!=string::npos)
			begin2=1;
	}

	return actualload(&block1,&block2,nmrcons);
}
      



 
int CPdb::actualload(vector<string> *b1, vector<string> *b2, vector <struct noeline> *nmrcons)
{
	string nmrseq;
	int adj,adj2;
	vector<int> out;
	char c1,c2;
	istringstream iss;
	string line,p;
	string resname1,atomname1,resname2,atomname2,resname,atomname;
	string oldatomname,oldatomname1,oldatomname2;
	int resid1,resid2,order,resid,oldresid,oldresid1,oldresid2;
	double a,b,c;
	int i,ii,iii;
	struct noeline tline;
	bool bfirst;
	int index,index_old;
	int norder1,norder2;


	bfirst=1;
	for(i=0;i<b1->size();i++)
	{
		
		iss.clear();
		iss.str(b1->at(i));
		iss>>index;
		iss>>p>>order>>p>>p;
		iss>>resid;
		iss>>resname;
		iss>>atomname;
		iss>>p>>oldresid>>p>>oldatomname>>p>>p;

		if(bfirst)
		{
			bfirst=0;
			resid1=resid2=-1;
			index_old=index;
			norder1=norder2=0;
		}

		if(index>index_old) //we now reach a new entry
		{
			tline.id=index_old;
			tline.resid1=resid1;
			tline.resid2=resid2;
			tline.resname1=resname1;
			tline.resname2=resname2;
			tline.atomname1=atomname1;
			tline.atomname2=atomname2;
			tline.oldresid1=oldresid1;
			tline.oldresid2=oldresid2;
			tline.oldatomname1=oldatomname1;
			tline.oldatomname2=oldatomname2;
			tline.multi=0;
			if(norder1>1 || norder2>1)
				tline.multi=1;
			nmrcons->push_back(tline);

			resid1=resid2=-1;
			index_old=index;
			norder1=norder2=0;
		}

		if(order==1)
		{
			norder1++;
			resid1=resid;
			resname1=resname;
			atomname1=atomname;
			oldresid1=oldresid;
			oldatomname1=oldatomname;
		}
		else if(order==2)
		{
			norder2++;
			resid2=resid;
			resname2=resname;
			atomname2=atomname;
			oldresid2=oldresid;
			oldatomname2=oldatomname;
		}
		else
		{
			cout<<"Unsupported order "<<order<<endl;
			break;
		}
	}

	//process the final entry.
	if(resid1>0 && resid2>0)
	{
		tline.id=index_old;
		tline.resid1=resid1;
		tline.resid2=resid2;
		tline.resname1=resname1;
		tline.resname2=resname2;
		tline.atomname1=atomname1;
		tline.atomname2=atomname2;
		nmrcons->push_back(tline);
		tline.oldresid1=oldresid1;
		tline.oldresid2=oldresid2;
		tline.oldatomname1=oldatomname1;
		tline.oldatomname2=oldatomname2;
	}


	cout<<"fillin sequence of NOE"<<endl;
	for(i=0;i<b2->size();i++)
	{
		//cout<<"i is "<<i<<" and b2's size is "<<b2->size()<<" and b1's size is "<<b1->size()<<" and nmrcons's size is "<<nmrcons->size()<<endl;
		iss.clear();
		iss.str(b2->at(i));
		iss>>index;

		iss>>p;
		iss>>p;
		iss>>p;
		iss>>p;
		iss>>p;
		iss>>p;
		iss>>p; b=atof(p.c_str());
		iss>>p; c=atof(p.c_str());
		iss>>p; a=atof(p.c_str());
		iss>>p;
		iss>>p;
		iii=-1;
		for(ii=0;ii<nmrcons->size();ii++)
		{
			if(index==nmrcons->at(ii).id)
				iii=ii;
		}

		// not found!
		if(iii<0)
			break;

		nmrcons->at(iii).a=a;
		nmrcons->at(iii).b=b;
		nmrcons->at(iii).c=c;
	}

	for(i=0;i<nmrcons->size();i++)
	{
		c1=Sequence::name2code(nmrcons->at(i).resname1);
		c2=Sequence::name2code(nmrcons->at(i).resname2);

		if(nmrcons->at(i).resid1>nmrseq.size())
		{
			nmrseq.resize(nmrcons->at(i).resid1,'U');
			//cout<<"resize nmrseq to "<<nmrcons->at(i).resid1+1<<endl;
		}
		//cout<<nmrcons->at(i).resid1<<" is filled with "<<c1<<"size of nmrseq is "<<nmrseq.size()<<endl;
		//cout.flush();
		nmrseq.at(nmrcons->at(i).resid1-1)=c1;
		//cout<<"done1"<<endl;
		if(nmrcons->at(i).resid2>nmrseq.size())
		{
			nmrseq.resize(nmrcons->at(i).resid2,'U');
			//cout<<"resize nmrseq to "<<nmrcons->at(i).resid2+1<<endl;
		}
		//cout<<nmrcons->at(i).resid2<<" is filled with "<<c2<<"size of nmrseq is "<<nmrseq.size()<<endl;
		//cout.flush();
		nmrseq.at(nmrcons->at(i).resid2-1)=c2;
		//cout<<"done2"<<endl;
	}

	//cout<<"just before alighment"<<endl;
	out=Sequence::align(pdbseq,nmrseq);
	adj=0;adj2=0;
	for(i=out.size()/10;i<out.size()-out.size()/10;i++)
	{
		if(out.at(i)!=0)
		{
			adj+=(out.at(i)-i);
			adj2++;
		}
	}
	adj/=adj2;
	//cout<<"Alignment of NOE data is "<<adj<<endl;


	//remove entry with multiply atoms
	for(i=nmrcons->size()-1;i>=0;i--)
	{
		if(nmrcons->at(i).multi==1)
		{
			cout<<"Remove noe entry with umbigirious assignment. "<<nmrcons->at(i).id<<endl;
			nmrcons->erase(nmrcons->begin()+i);
		}
	}


	for(i=0;i<nmrcons->size();i++)
	{
        nmrcons->at(i).group=1;
		nmrcons->at(i).index1.length=0.0;
		nmrcons->at(i).index2.length=0.0;
		nmrcons->at(i).resid1-=adj;
		nmrcons->at(i).resid2-=adj;

		if(nmrcons->at(i).resid1>=1 && nmrcons->at(i).resid2>=1 && nmrcons->at(i).resid1<=(int)v.size() && nmrcons->at(i).resid2<=(int)v.size() )
		{
			nmrcons->at(i).index1=v.at(nmrcons->at(i).resid1-1)->query(nmrcons->at(i).atomname1);
			nmrcons->at(i).index2=v.at(nmrcons->at(i).resid2-1)->query(nmrcons->at(i).atomname2);
		}
		else
		{
			nmrcons->erase(nmrcons->begin()+i);
			i--;
		}
	}
	return adj;

}

void CPdb::outputnoe(string filename, vector <struct noeline> * nmrcons)
{
	int i;
	ofstream fout;

	fout.open(filename.c_str());

	for(i=0;i<nmrcons->size();i++)
	{
		fout<<nmrcons->at(i).group<<" ";
		fout<<nmrcons->at(i).resid1<<" ";
		fout<<nmrcons->at(i).resname1<<" ";
		fout<<nmrcons->at(i).atomname1<<" ";
		fout<<nmrcons->at(i).resid2<<" ";
		fout<<nmrcons->at(i).resname2<<" ";
		fout<<nmrcons->at(i).atomname2<<" ";
		fout<<nmrcons->at(i).b<<" ";
		fout<<nmrcons->at(i).c<<" ";
		fout<<nmrcons->at(i).a;
		fout<<endl;
	}
	fout.close();

	return;
}

void CPdb::inputnoe(string filename, vector <struct noeline> * nmrcons)
{
	int i;
	ifstream fin;
	struct noeline noe;
	string line,p;
	vector<string> ps;

	fin.open(filename.c_str());

	istringstream iss;

	i=0;
	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		ps.clear();
		
		while(iss>>p)
		{
			ps.push_back(p);
		}

		if(ps.size()>10)
		{
			noe.group=atoi(ps.at(0).c_str());
			noe.id=i++;
			noe.multi=0;
			noe.resid1=atoi(ps.at(2).c_str());
			noe.resname1=ps.at(3);
			noe.atomname1=ps.at(4);
			noe.resid2=atoi(ps.at(5).c_str());
			noe.resname2=ps.at(6);
			noe.atomname2=ps.at(7);
			noe.b=atof(ps.at(8).c_str());
			noe.c=atof(ps.at(9).c_str());
			noe.a=atof(ps.at(10).c_str());
			nmrcons->push_back(noe);
		}
	}
	fin.close();


	for(i=0;i<nmrcons->size();i++)
	{
		nmrcons->at(i).index1.length=0.0;
		nmrcons->at(i).index2.length=0.0;
		if(nmrcons->at(i).resid1>=1 && nmrcons->at(i).resid2>=1 && nmrcons->at(i).resid1<=(int)v.size() && nmrcons->at(i).resid2<=(int)v.size() )
		{
			nmrcons->at(i).index1=v.at(nmrcons->at(i).resid1-1)->query(nmrcons->at(i).atomname1);
			nmrcons->at(i).index2=v.at(nmrcons->at(i).resid2-1)->query(nmrcons->at(i).atomname2);
		}
	}

	return;
}




int CPdb::load_exactnoe(string filename,vector <struct noeline> *nmrcons)
{
	ifstream fin(filename.c_str());
	string line,p;
	istringstream iss;
	struct noeline tline;
	string nmrseq;
	char c1,c2;
	int i;

	i=0;
	while(getline(fin,line))
	{
		i++;
		tline.id=i;
		iss.clear();
		iss.str(line);
		iss>>tline.resid1>>tline.resname1>>tline.atomname1;
		iss>>tline.resid2>>tline.resname2>>tline.atomname2;
		iss>>tline.a>>tline.c;
		nmrcons->push_back(tline);
	}

	for(i=0;i<nmrcons->size();i++)
	{
		c1=Sequence::name2code(nmrcons->at(i).resname1);
		c2=Sequence::name2code(nmrcons->at(i).resname2);
		if(nmrcons->at(i).resid1>=nmrseq.size())
		{
			nmrseq.resize(nmrcons->at(i).resid1+1,'U');
		}
		nmrseq.at(nmrcons->at(i).resid1)=c1;
		if(nmrcons->at(i).resid2>=nmrseq.size())
		{
			nmrseq.resize(nmrcons->at(i).resid2+1,'U');
		}
		nmrseq.at(nmrcons->at(i).resid2)=c2;
	}

	int adj,adj2;
	string out1,out2,out3;
	vector<int> out=Sequence::aligno(pdbseq,nmrseq,out1,out2,out3);
	adj=0;adj2=0;
	for(i=out.size()/10;i<out.size()-out.size()/10;i++)
	{
		if(out.at(i)!=0)
		{
			adj+=(out.at(i)-i);
			adj2++;
		}
	}
	adj/=adj2;
	cout<<out1<<endl;
	cout<<out2<<endl;
	cout<<out3<<endl;



	for(i=0;i<nmrcons->size();i++)
	{
		nmrcons->at(i).index1.length=0.0;
		nmrcons->at(i).index2.length=0.0;
		nmrcons->at(i).resid1-=adj;nmrcons->at(i).resid1++;
		nmrcons->at(i).resid2-=adj;nmrcons->at(i).resid2++;

		if(nmrcons->at(i).resid1>=1 && nmrcons->at(i).resid2>=1 && nmrcons->at(i).resid1<=(int)v.size() && nmrcons->at(i).resid2<=(int)v.size() )
		{
			nmrcons->at(i).index1=v.at(nmrcons->at(i).resid1-1)->query(nmrcons->at(i).atomname1);
			nmrcons->at(i).index2=v.at(nmrcons->at(i).resid2-1)->query(nmrcons->at(i).atomname2);
		}
		else
		{
			nmrcons->erase(nmrcons->begin()+i);
			i--;
		}
	}

	return 0;
}





int CPdb::loaddihecons(string filename,vector<struct diheline> *dihes)		
{
	ifstream fin(filename.c_str());
	string line;
	vector<string> block;
	bool bstart,inloop,begin;


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
			dihecons_actualload(&block,dihes);
			block.clear();
		}

		if(begin==1)
			block.push_back(line);


		if(line.find("loop_")!=string::npos)
			inloop=1;
		if(line.find("stop_")!=string::npos)
			inloop=0;

		if(inloop==0 && line.find("save_CNS/XPLOR_dihedral")!=string::npos)
			bstart=1;

		if(bstart==1 && inloop==1 && line.find("_Torsion_angle_constraint.Torsion_angle_constraint_list_ID")!=string::npos)
			begin=1;
	}

	if(block.size()>0) 
	{
		dihecons_actualload(&block,dihes);
		block.clear();
	}
	return 0;
}


int CPdb::dihecons_actualload(vector<string>* block,vector<struct diheline> *dihes)
{
	istringstream iss;
	string p;
	int i,j,m;
	string nmrseq;
	vector<int> out;
	char c1,c2;
	int adj,adj2;
	bool bremove;

	dihes->resize(block->size());
	for(i=0;i<block->size();i++)
	{
		iss.clear();
		iss.str(block->at(i));
		iss>>dihes->at(i).id;
		iss>>p;
		for(j=0;j<4;j++)
		{
			iss>>p>>p;
			iss>>m;dihes->at(i).resid.push_back(m);
			iss>>p;dihes->at(i).resname.push_back(p);
			iss>>p;dihes->at(i).atomname.push_back(p);
		}
		iss>>dihes->at(i).upper>>dihes->at(i).lower;
		dihes->at(i).middle=(dihes->at(i).lower+dihes->at(i).upper)/2.0;
		dihes->at(i).delta=fabs(dihes->at(i).upper-dihes->at(i).lower)/2.0;

		c1=Sequence::name2code(dihes->at(i).resname[0]);
		c2=Sequence::name2code(dihes->at(i).resname[3]);
		for(j=0;j<4;j++)	dihes->at(i).resid[j]--;
		if(dihes->at(i).resid[0]+1>nmrseq.size())
			nmrseq.resize(dihes->at(i).resid[0]+1,'U');
		nmrseq.at(dihes->at(i).resid[0])=c1;
		if(dihes->at(i).resid[3]+1>nmrseq.size())
			nmrseq.resize(dihes->at(i).resid[3]+1,'U');
		nmrseq.at(dihes->at(i).resid[3])=c2;
	}

	out=Sequence::align(pdbseq,nmrseq);
	adj=0;adj2=0;
	for(i=0;i<out.size();i++)
	{
		if(out.at(i)!=0)
		{
			adj+=(out.at(i)-i);
			adj2++;
		}
	}
	adj/=adj2;

	for(i=dihes->size()-1;i>=0;i--)
	{
		bremove=0;
		for(j=0;j<4;j++)
		{
			dihes->at(i).resid[j]++;
			dihes->at(i).resid[j]-=adj;
			if(dihes->at(i).resid[j]<1 || dihes->at(i).resid[j]>v.size())
				bremove=1;
		}
		if(bremove)
		{
			cerr<<"Remove dihes cons ID "<<dihes->at(i).id<<" because it is out-of-bound"<<endl;
			dihes->erase(dihes->begin()+i);
		}
	}




	for(i=0;i<dihes->size();i++)
	{
		dihes->at(i).index.resize(4);
		for(j=0;j<4;j++)
			dihes->at(i).index[j]=v.at(dihes->at(i).resid[j]-1)->query(dihes->at(i).atomname[j]);
	}
	return 0;
}

int CPdb::buildpdb(string seq)
{	
	cout << "BUILDPDB" << endl;
	string residue;
	CAminoacid *t;
	int nres;
	int i;

	for(i=0;i<seq.size();i++)
	{
		residue=Sequence::code2name(seq.at(i));
		t=NULL;
		if(residue=="ALA")	t=new CAla;
		else if(residue=="ARG")	t=new CArg;
		else if(residue=="ASN")	t=new CAsn;
		else if(residue=="ASP")	t=new CAsp;
		else if(residue=="GLN")	t=new CGln;
		else if(residue=="GLU")	t=new CGlu;
		else if(residue=="GLY")	t=new CGly;
		else if(residue=="HIS")	t=new CHis;
		else if(residue=="HID")	t=new CHis;
		else if(residue=="HIE")	t=new CHis;
		else if(residue=="ILE")	t=new CIle;
		else if(residue=="LEU") t=new CLeu;
		else if(residue=="LYS")	t=new CLys;
		else if(residue=="LYP")	t=new CLys;
		else if(residue=="MET")	t=new CMet;
		else if(residue=="PHE")	t=new CPhe;
		else if(residue=="PRO")	t=new CPro;
		else if(residue=="SER")	t=new CSer;
		else if(residue=="THR")	t=new CThr;
		else if(residue=="TYR")	t=new CTyr;
		else if(residue=="VAL")	t=new CVal;
		else if(residue=="TRP")	t=new CTrp;
		else if(residue=="CYS")	t=new CCys;
		else if(residue=="CYN")	t=new CCys;
		else if(residue=="CYX")	t=new CCyx;
		else if(residue=="UNK")	t=new CUnk;
		else t=new CUnk;

		t->setresidue(i+1);

		v.push_back(t);
	}

	nres=v.size();
	natom=0;
	natom2=0;
	chains.push_back(v.size());

	return natom;
}

void CPdb::setup(CAminoacid *t, CLigand *tt, int iligand, string residue , int index_old,vector<string> block)
{
	return;
}

// Updated function for OpenACC
// Some extra arrays are created to avoid having to load complex data structures onto GPU
int CPdb::loadpdb(string filename)
{
	double st = omp_get_wtime();

	string line,part;
	ifstream fin(filename.c_str());
	pdbfilename=filename;
	bool first;
	int iligand;
	int atomindex;
	int index,index_old;
	string chain,chain_old;
	string residue,residue_old;
	struct proteinblock pdbblock;
	string add;
	char buffer[10];
	

	cout<<"load in "<<pdbfilename<<endl;

	clear(); //clear all data before load in new PDB information
	iligand=0;
	atomindex=0;
	first=1;
	index_old=1;
	chain_old=" ";
	while(getline(fin,line))
	{
		part=line.substr(0,6);
		if(part=="SEQADV")
		{

		}



		part=line.substr(0,6);
		if(part=="END" || part=="ENDMDL")
		{
			break;
		}

		part=line.substr(0,3);
		if(part=="TER")
		{
			//following is ligand or water in most case
			//or another chain! need to be taken care of this situation
			iligand=1;
			continue;
		}

		//neglect all entries that are not ATOM or HETATM
		part=line.substr(0,6);
		if(part!="ATOM  " && part!="HETATM")
			continue;


		//ignore water molecules.
		part=line.substr(17,3);
		if(part=="HOH" || part=="WAT")
			continue;

		
		index=atoi(line.substr(22,4).c_str());
		chain=line.substr(21,1);
		residue=line.substr(17,3);

		//Put in real atom index, this is important for coor processing !!
		atomindex++;
		sprintf(buffer,"%5d",atomindex);
		add=buffer;
		line.replace(6,5,add);

		if(first==1)
		{
			first=0;
			index_old=index;
			chain_old=chain;
			residue_old=residue;
			pdbblock.block.push_back(line);
			pdbblock.iligand=iligand;
			pdbblock.residue=residue;
			pdbblock.index=index;
		}
		else if(index==index_old && chain==chain_old )
		{
			pdbblock.block.push_back(line);
		}
		else if(chain==chain_old && index>index_old ) // new residue
		{
			blocks.push_back(pdbblock);
			pdbblock.block.clear();
			pdbblock.block.push_back(line);  
			pdbblock.iligand=iligand;
			pdbblock.index=index;
			pdbblock.residue=residue;
			index_old=index;
			residue_old=residue;
		}
		else if(chain==chain_old && index<index_old ) //most likely a new chain
		{
			blocks.push_back(pdbblock);
			chain_block.push_back(blocks.size());
			iligand=0;
			pdbblock.block.clear();
			pdbblock.block.push_back(line);	
			pdbblock.iligand=iligand;
			pdbblock.index=index;
			pdbblock.residue=residue;
			index_old=index;
			residue_old=residue;
			chain_old=chain;
		}
		else if(chain!=chain_old)
		{
			blocks.push_back(pdbblock);
			chain_block.push_back(blocks.size());
			iligand=0;
			pdbblock.block.clear();
			pdbblock.block.push_back(line);	
			pdbblock.iligand=iligand;
			pdbblock.index=index;
			pdbblock.residue=residue;
			index_old=index;
			residue_old=residue;
			chain_old=chain;
		}
	}
	
	if(pdbblock.block.size()>0)
		blocks.push_back(pdbblock);
	chain_block.push_back(blocks.size());

	//start part 2 here
	int start,begin,stop,stop2;
	CAminoacid *t,*t2;
	CLigand *tt;
	int n;
	int i,j,ii;
	natom=natom2=0;


	for(ii=0;ii<chain_block.size();ii++)
	{
		if(ii==0)
			start=0;
		else
			start=chain_block.at(ii-1);
		stop=chain_block.at(ii);
		n=0;
		stop2=stop;


		pdbblock=blocks.at(stop-1);
		if(pdbblock.iligand==0)
		{
			//even without TER, res has only one line is mostly a ligand,if it is the last entry
			if(pdbblock.block.size()==1 )
			{
					blocks.at(stop-1).iligand=1;
					stop2=stop-1;
			}
		}

		for(i=stop-2;i>=start;i--)
		{
			pdbblock=blocks.at(i);
			if(pdbblock.iligand==0 && blocks.at(i+1).iligand==1)
			{
				//even without TER, res has only one line is mostly a ligand, if following one is also a ligand
				if(pdbblock.block.size()==1 )
				{
					blocks.at(i).iligand=1;
					stop2=i;
				}
			}
		}


		for(i=start;i<stop;i++)
		{
			
			pdbblock=blocks.at(i);
			if(pdbblock.iligand==0)
			{
				residue=pdbblock.residue;
				if(residue=="ALA")	t=new CAla;
				else if(residue=="ARG")	t=new CArg;
				else if(residue=="ASN")	t=new CAsn;
				else if(residue=="ASP")	t=new CAsp;
				else if(residue=="GLN")	t=new CGln;
				else if(residue=="GLU")	t=new CGlu;
				else if(residue=="GLY")	t=new CGly;
				else if(residue=="HIS")	t=new CHis;	
				else if(residue=="HID")	t=new CHis;
				else if(residue=="HIE")	t=new CHis;
				else if(residue=="ILE")	t=new CIle;
				else if(residue=="LEU")	t=new CLeu;
				else if(residue=="LYS")	t=new CLys;
				else if(residue=="LYP")	t=new CLys;
				else if(residue=="MET")	t=new CMet;
				else if(residue=="PHE")	t=new CPhe;
				else if(residue=="PRO")	t=new CPro;
				else if(residue=="SER")	t=new CSer;
				else if(residue=="THR")	t=new CThr;
				else if(residue=="TYR")	t=new CTyr;
				else if(residue=="VAL")	t=new CVal;
				else if(residue=="TRP") t=new CTrp;
				else if(residue=="CYS") t=new CCys;
				else if(residue=="CYN") t=new CCys;
				else if(residue=="CYX") t=new CCyx;
				else if(iligand>=2) tt=new CLigand;
				else
				{
					t=new CUnk;
					#ifndef IGNORE_UNKNOWN
					cout<<"Warning! unrecognized residue name "<<residue<<" for residue "<<index_old<<endl;
					printblock(pdbblock.block);
					#endif
				}

				n++;t->setresidue(n);
				if(i==start) t->setnterminal();
				if(i==stop2-1) t->setcterminal();
				else if(blocks.at(i+1).iligand==1) t->setcterminal();
				t->process(pdbblock.block);
				natom+=pdbblock.block.size();
				
				if(i>start)
				{
					for(j=blocks.at(i-1).index+1;j<pdbblock.index;j++)
					{
						t2=new CMiss;
						n++;t2->setresidue(n);
						v.push_back(t2);
						nmiss++;
					}
				}
				v.push_back(t);
			}
			else //ligand
			{
				tt=new CLigand;
				n++;tt->setresidue(n);
				tt->process(pdbblock.block);
				natom2+=pdbblock.block.size();
				ligand.push_back(tt);
			}
		}//for(i=start;i<stop;i++)

		chains.push_back(v.size());
		chain_ligand.push_back(ligand.size());
	}



	//re-index residues
	for(j=0;j<(int)chains.size();j++)
	{
		if(j==0)
			begin=0;
		else
			begin=chains.at(j-1);
		stop=chains.at(j);

		for(i=begin;i<stop;i++)
		{
			v.at(i)->setresidue(i+1);
		}
	}


	for(j=0;j<(int)chains.size();j++)
	{
		if(j==0)
			begin=0;
		else
			begin=chains.at(j-1);
		stop=chains.at(j);

		//set up chain information
		for(i=begin;i<stop;i++)
			v.at(i)->chain=j;

		for(i=begin;i<stop-1;i++)
		{	
			v.at(i)->setfollowingn(v.at(i+1)->get("N").index);
		}
		v.at(stop-1)->setfollowingn(v.at(stop-1)->get("O").index);
		
		for(i=begin+1;i<stop;i++)
		{
			v.at(i)->setpreviousc(v.at(i-1)->get("C").index);
		}
		v.at(begin)->setpreviousc(v.at(begin)->get("H1").index);

		//residue index from the PDB file may be inconsistent
	}

	for(i=0;i<(int)v.size();i++)
	{
		pdbseq.push_back(v.at(i)->OneLetterName);
	}

////////// Begin updated code //////////////////////////////////////////
	v_oneletternames = new char[v.size()]; // Remove from complex data structure for easier GPU use
	code_pos = new int[v.size()]; // Remove from complex data structure for easier GPU use
	for(i=0;i<v.size();i++){
		v_oneletternames[i]=v.at(i)->OneLetterName;
		code_pos[i]=Sequence::code2pos(v_oneletternames[i]);
	}
	v_arr = v.data(); // Grab point to underlying array which are better suited for GPU
	v_size = v.size();
////////// End updated code ////////////////////////////////////////////

	cout << "loadpdb: " << omp_get_wtime() - st << " seconds" << endl;
	return (natom+natom2);

}


		
// Altered to use "END" or "ENDMDL" to support some PDB files
int CPdb::loadpdb_old(string filename)
{
	cout << "LOADPDB_OLD" << endl;
	bool bres;
	int i,j,n;
	int begin,stop;
	int atomindex;
	int iligand;
	vector<string> block;
	string line, part, residue,chain,chain_old,add;
	int index,index_old;
	bool first;
	char buffer[6];

	pdbfilename=filename;

	clear();
	
	ifstream fin(filename.c_str());
	CAminoacid *t;
	CLigand *tt;

	
	natom=natom2=0;
	first=1;
	atomindex=0;
	iligand=0;
	nmiss=0;


	cout<<"load in "<<pdbfilename<<endl;
	while(getline(fin,line))
	{
		part=line.substr(0,6);
		if(part=="END" || part=="ENDMDL")
		{
			break;
		}

		part=line.substr(0,3);
		if(part=="TER")
		{
			//following is ligand or water in most case
			//or another chain! need to be taken care of this situation
			iligand=1;
			continue;
		}

		part=line.substr(0,6);
		if(part!="ATOM  " && part!="HETATM")
			continue;

		part=line.substr(17,3);
		if(part=="HOH" || part=="WAT")
			continue;

		atomindex++;
		index=atoi(line.substr(22,4).c_str());
		chain=line.substr(21,1);


		if(first==1)
		{
			first=0;
			index_old=index;
			chain_old=chain;
			residue=line.substr(17,3);
			n=0;
		}
		bres=0;

		if(index==index_old  && chain==chain_old) //same residue
		{
			if(line.substr(17,3)==residue)
			{
				sprintf(buffer,"%5d",atomindex);
				add=buffer;
				line.replace(6,5,add);
				block.push_back(line);
			}
			else
			{
				#ifndef IGNORE_UNKNOWN
				cout<<"Warning: same residue index but different residue name!   ";
				cout<<line<<endl;
				#endif
			}
		}
		else if(chain != chain_old) //new chain, so that there is also new residue
		{
			chains.push_back((int)v.size()+1);
			if(chain.compare("A")!=0 ) iligand=0;
			bres=1;
			n=0;
		}
		else if(index!=index_old)//same chain, but new residue.
		{
			bres=1;
		}
		else //must be error, so we break;
		{
			cout<<"Error in processing PDB file! Processing aborted!"<<endl;
			break;
		}


		if(bres==1)  //
		{
			t=NULL;
			tt=NULL;
		if(residue=="ALA")	t=new CAla;
		else if(residue=="ARG")	t=new CArg;
		else if(residue=="ASN")	t=new CAsn;
		else if(residue=="ASP")	t=new CAsp;
		else if(residue=="GLN")	t=new CGln;
		else if(residue=="GLU")	t=new CGlu;
		else if(residue=="GLY")	t=new CGly;
		else if(residue=="HIS")	t=new CHis;	
		else if(residue=="HID")	t=new CHis;
		else if(residue=="HIE")	t=new CHis;
		else if(residue=="ILE")	t=new CIle;
		else if(residue=="LEU")	t=new CLeu;
		else if(residue=="LYS")	t=new CLys;
		else if(residue=="LYP")	t=new CLys;
		else if(residue=="MET")	t=new CMet;
		else if(residue=="PHE")	t=new CPhe;
		else if(residue=="PRO")	t=new CPro;
		else if(residue=="SER")	t=new CSer;
		else if(residue=="THR")	t=new CThr;
		else if(residue=="TYR")	t=new CTyr;
		else if(residue=="VAL")	t=new CVal;
		else if(residue=="TRP") t=new CTrp;
		else if(residue=="CYS") t=new CCys;
		else if(residue=="CYN") t=new CCys;
		else if(residue=="CYX") t=new CCyx;
		else if(iligand>=2) tt=new CLigand;
		else
		{
			t=new CUnk;
			#ifndef IGNORE_UNKNOWN
			cout<<"Warning! unrecognized resiude name "<<residue<<" fore residue "<<index_old<<endl;
			printblock(block);
			#endif
		}

			if(t!=NULL)
			{
				n++;t->setresidue(n);
				if(n==1) t->setnterminal();
				t->process(block);
				natom+=block.size();
				v.push_back(t);
				if(iligand>0)
					iligand++;
				//cout<<natom<<endl;
			}

			if(tt!=NULL)
			{
				tt->process(block);
				natom2+=block.size();
				ligand.push_back(tt);
			}
			
			//insert missing residues if necessary
			if(iligand<2)
			{
				for(i=index_old+1;i<index;i++)
				{
					t=new CMiss;
					n++;t->setresidue(n);
					v.push_back(t);
					nmiss++;
				}
			}

			

			//prepare for new residues
			index_old=index;
			chain_old=chain;
			block.clear();
			index=atoi(line.substr(22,4).c_str());
			chain=line.substr(21,1);

			sprintf(buffer,"%5d",atomindex);
			add=buffer;
			line.replace(6,5,add);

			block.push_back(line);
			residue=line.substr(17,3);
		}

	}

	

	if(block.size()>0)
	{
		if(iligand>=1)
			iligand++;
		
		t=NULL;
		tt=NULL;
		
		if(residue=="ALA")	t=new CAla;
		else if(residue=="ARG")	t=new CArg;
		else if(residue=="ASN")	t=new CAsn;
		else if(residue=="ASP")	t=new CAsp;
		else if(residue=="GLN")	t=new CGln;
		else if(residue=="GLU")	t=new CGlu;
		else if(residue=="GLY")	t=new CGly;
		else if(residue=="HIS")	t=new CHis;	
		else if(residue=="HID")	t=new CHis;
		else if(residue=="HIE")	t=new CHis;
		else if(residue=="ILE")	t=new CIle;
		else if(residue=="LEU")	t=new CLeu;
		else if(residue=="LYS")	t=new CLys;
		else if(residue=="LYP")	t=new CLys;
		else if(residue=="MET")	t=new CMet;
		else if(residue=="PHE")	t=new CPhe;
		else if(residue=="PRO")	t=new CPro;
		else if(residue=="SER")	t=new CSer;
		else if(residue=="THR")	t=new CThr;
		else if(residue=="TYR")	t=new CTyr;
		else if(residue=="VAL")	t=new CVal;
		else if(residue=="TRP") t=new CTrp;
		else if(residue=="CYS") t=new CCys;
		else if(residue=="CYN") t=new CCys;
		else if(residue=="CYX") t=new CCyx;
		else if(iligand>=2) tt=new CLigand;
		else
		{
			t=new CUnk;
			#ifndef IGNORE_UNKNOWN
			cout<<"Warning! unrecognized resiude name "<<residue<<" fore residue "<<index_old<<endl;
			printblock(block);
			#endif
		}

		if(t!=NULL)
		{
			n++;t->setresidue(n);
			t->process(block);
			natom+=block.size();
			v.push_back(t);
		}
		if(tt!=NULL)
		{
			tt->process(block);
			ligand.push_back(tt);
			natom2+=block.size();
		}
		chains.push_back((int)v.size());
	}

	//re-index residues
	for(j=0;j<(int)chains.size();j++)
	{
		if(j==0)
			begin=0;
		else
			begin=chains.at(j-1);
		stop=chains.at(j);

		for(i=begin;i<stop;i++)
		{
			v.at(i)->setresidue(i+1-begin);
		}
	}


	for(j=0;j<(int)chains.size();j++)
	{
		if(j==0)
			begin=0;
		else
			begin=chains.at(j-1);
		stop=chains.at(j);

		for(i=begin;i<stop-1;i++)
		{	
			v.at(i)->setfollowingn(v.at(i+1)->get("N").index);
		}
		v.at(stop-1)->setfollowingn(v.at(stop-1)->get("O").index);
		
		for(i=begin+1;i<stop;i++)
		{
			v.at(i)->setpreviousc(v.at(i-1)->get("C").index);
		}
		v.at(begin)->setpreviousc(v.at(begin)->get("H1").index);

		//residue index from the PDB file may be inconsistent
	}

	for(i=0;i<(int)v.size();i++)
	{
		pdbseq.push_back(v.at(i)->OneLetterName);
	}

	return (natom+natom2);
}


bool CPdb::bisunknownmissing()
{
	int i;
	bool flag1,flag2,flag3;

	flag1=0;
	flag2=0;
	flag3=0;

	for(i=0;i<(int)v.size();i++)
	{
		if(v.at(i)->OneLetterName!='X' && flag2==0)
		{
			flag1=1;
			continue;
		}
		if(v.at(i)->OneLetterName=='X' && flag1==1)
		{
			flag1=0;
			flag2=1;
			continue;
		}
		if(v.at(i)->OneLetterName!='X' && flag2==1)
		{
			flag3=1;
			flag1=0;
			flag2=0;
			break;
		}
	}

	return flag3;
}

		



void CPdb::printpdb(char * filename, int flag)
{
	int i,n;
	FILE *fp;

	fp=fopen(filename,"wt");
	n=1;
	for(i=0;i<v.size();i++)
	{
		if(flag==0) v.at(i)->printpdb(fp,atomname,x,y,z,b,n);
		else v.at(i)->printpdb(fp,atomname,x,y,z,b2,n);
	}
	fclose(fp);
	return;
}
	



void CPdb::printblock(vector<string> block)
{
	int i;
	cout<<"*********************************************************"<<endl;
	for(i=0;i<(int)block.size();i++)
		cout<<block.at(i)<<endl;
	cout<<"*********************************************************"<<endl;
	return;
}

vector<int> CPdb::getselect(string input)
{
	int n,n1,n2;
	vector<int> res;
	vector<int> sel;
	vector<string> atoms;
	size_t found_res,found;
	string part1,part2,part3;
	string p,p1,p2;


	stringstream ss;

	try
	{
		found_res=input.find(':');
		if(found_res==string::npos)
			throw "Can't find :";
		if(found_res!=0)
			throw "First charactor is not :";
		found_res=input.find('@');
		//cout<<"input is "<<input<<" found_res is "<<found_res<<endl;
		if(found_res==string::npos)
			throw "Can't find @";
		part1=input.substr(1,found_res-1);
		part2=input.substr(found_res+1,input.length()-found_res-1);
		//cout<<part1<<endl;
		//cout<<part2<<endl;

		//extract & part to part3
		found_res=input.find('&');
		if(found_res!=string::npos)
		{
			part3=part1.substr(found_res,part1.length()-found_res);
			part1=part1.substr(0,found_res-1);
		}


		//processing part1 now
		ss.str(part1);
		while(getline(ss,p,','))
		{
			found=p.find('-');
			if(found==string::npos)
			{
				if(p=="%")
					n=v.size();
				else
					n=atoi(p.c_str());
				if(n<=0)
					throw "Residue number <=0 ??";
				res.push_back(n);
			}
			else
			{
				p1=p.substr(0,found);
				p2=p.substr(found+1,p.length()-found);
				if(p1=="%")
					n1=v.size();
				else
					n1=atoi(p1.c_str());
				if(p2=="%")
					n2=v.size();
				else
					n2=atoi(p2.c_str());
				if(n1<=0 || n2<=0)
					throw "Residue number <=0 ??";
				for(n=n1;n<=n2;n++)
					res.push_back(n);

			}
		}
		//process part3 now
		if(part3.size()>0)
		{
			for(int i=res.size()-1;i>=0;i--)
			{
				n=res.at(i)-1;
				if(n<0)
					throw "Res index is less than 1";
				if(n>=v.size())
					throw "Res index is larger than nres";
				char c=v.at(n)->getdssp();
				if(part3.find(c)==string::npos)
					res.erase(res.begin()+i);
			}
		}

		//precess part2 now
		ss.clear();
		ss.str(part2);
		while(getline(ss,p,','))
		{
			atoms.push_back(p);
		}
	


		unsigned int i,j;
		for(i=0;i<(int)res.size();i++)
		{
			n=res.at(i)-1;
			if(n<0)
				throw "Res index is less than 1";
			if(n>=v.size())
				throw "Res index is larger than nres";

			for(j=0;j<atoms.size();j++)
			{
				
				if(atoms.at(j).compare("allheavy")==0)
				{
					v.at(n)->bbheavycoor(&sel);
					v.at(n)->heavycoor(&sel);
				}
				else if(atoms.at(j).compare("all")==0)
				{
					v.at(n)->allcoor(&sel);
				}
				else if(atoms.at(j).compare("scheavy")==0)
				{
					v.at(n)->heavycoor(&sel);
				}
				else if(atoms.at(j).compare("bbheavy")==0)
				{
					v.at(n)->bbheavycoor(&sel);
				}				
				else
				{
					n1=v.at(n)->get(atoms.at(j).c_str()).index;
					if(n1>0)
						sel.push_back(n1);
					else
						cout<<"Res "<<n+1<<v.at(n)->OneLetterName<<" doesn't have "<<atoms.at(j).c_str()<<" atom"<<endl;
				}
			}
		}
	}


	catch (const char * str)
	{
		cout<<"Input selection phasing exception raised because "<<str<<". Please check."<<endl;
		cout<<"select all Ca atoms of all residues instead!"<<endl;
		sel.clear();getca(&sel);
	}
	return sel;
}


vector<int> CPdb::getselect(string input,vector<int> res)
{
	int n,n1;
	vector<int> sel;
	vector<string> atoms;
	string part1,part2;
	string p,p1,p2;


	stringstream ss;


	ss.clear();
	ss.str(input);
	while(getline(ss,p,','))
	{
		atoms.push_back(p);
	}
	

	unsigned int i,j;
	for(i=0;i<(int)res.size();i++)
	{
		n=res.at(i)-1;

		for(j=0;j<atoms.size();j++)
		{
			if(atoms.at(j).compare("allheavy")==0)
			{
				v.at(n)->bbheavycoor(&sel);
				v.at(n)->heavycoor(&sel);
			}
			else if(atoms.at(j).compare("all")==0)
			{
				v.at(n)->allcoor(&sel);
			}
			else if(atoms.at(j).compare("scheavy")==0)
			{
				v.at(n)->heavycoor(&sel);
			}
			else if(atoms.at(j).compare("bbheavy")==0)
			{
				v.at(n)->bbheavycoor(&sel);
			}					
			else
			{
				n1=v.at(n)->get(atoms.at(j).c_str()).index;
				if(n1>0)
					sel.push_back(n1);
				else
					cout<<"Res "<<n+1<<v.at(n)->OneLetterName<<" doesn't have "<<atoms.at(j).c_str()<<" atom"<<endl;
			}
		}
	}
	return sel;
}



vector<int> CPdb::getselect(vector<int> &num, string input)
{
    int n,n1;
    unsigned int j;
    vector<int> sel;
    vector<string> atoms;
    string part1,part2;
    string p,p1,p2;
    
    
    stringstream ss;
    
    
    ss.clear();
    ss.str(input);
    while(getline(ss,p,','))
    {
        atoms.push_back(p);
    }
    
    
    
    for(n=0;n<v.size();n++)
    {
        
        for(j=0;j<atoms.size();j++)
        {
            if(atoms.at(j).compare("allheavy")==0)
            {
                v.at(n)->bbheavycoor(&sel);
                v.at(n)->heavycoor(&sel);
            }
            else if(atoms.at(j).compare("all")==0)
            {
                v.at(n)->allcoor(&sel);
            }
            else if(atoms.at(j).compare("scheavy")==0)
            {
                v.at(n)->heavycoor(&sel);
            }
            else if(atoms.at(j).compare("bbheavy")==0)
            {
                v.at(n)->bbheavycoor(&sel);
            }					
            else
            {
                n1=v.at(n)->get(atoms.at(j).c_str()).index;
                if(n1>0)
                    sel.push_back(n1);
                else
                    cout<<"Res "<<n+1<<v.at(n)->OneLetterName<<" doesn't have "<<atoms.at(j).c_str()<<" atom"<<endl;
            }
        }
        num.push_back(sel.size());
    }
    return sel;
}



void CPdb::getca(vector<int> *index)
{
	unsigned int i;
	for(i=0;i<(int)v.size();i++)
	{
		//if(i>40 && i<87) continue;
		index->push_back(v.at(i)->getca());
	}

	return;
}

vector<int> CPdb::getselectca(vector<int> res)
{
	unsigned int i;
	int n;
	vector<int> index;

	for(i=0;i<(int)res.size();i++)
	{
		n=res.at(i);
		if(n>=1 && n<=v.size()) index.push_back(v.at(n-1)->getca());
	}
	return index;
}




void CPdb::getdihe(vector<dihe_group> *t, vector<int> *n)
{
	double st = omp_get_wtime();
	unsigned int i;

	for(i=0;i<(int)v.size();i++)
	{
		v[i]->bbdihe(t);	
		v[i]->dihe(t);
		n->push_back(t->size());
	}
	cout << "pdb::getdihe: " << omp_get_wtime()-st << " seconds" << endl;
}


void CPdb::getbbdihe(vector<dihe_group> *t)
{
	unsigned int i;

	for(i=0;i<(int)v.size();i++)
	{
		v[i]->bbdihe(t);	
	}
}


void CPdb::getbbdihe(vector<dihe_group> *t,int *ibb)
{
	unsigned int i;

	for(i=0;i<(int)v.size();i++)
	{
		if( i==0 || i==(int)v.size()-1)
		{
			ibb[i]=0;
		}
		else
		{
			v[i]->bbdihe(t);	
			ibb[i]=1;
		}
	}
}

void CPdb::getbbdihe_nopro(vector<dihe_group> *t,int *ibb)
{
	unsigned int i;

	for(i=0;i<(int)v.size();i++)
	{
		if(v.at(i)->OneLetterName=='P' || i==0 || i==(int)v.size()-1)
		{
			ibb[i]=0;
		}
		else
		{
			v[i]->bbdihe(t);	
			ibb[i]=1;
		}
	}
}



void CPdb::getdihe(vector<dihe_group> *t)
{
	unsigned int i;


	for(i=0;i<(int)v.size();i++)
	{
		v[i]->bbdihe(t);	
		v[i]->dihe(t);
	}
}

void CPdb::getbb(vector<bb_group> *t)
{
	double st = omp_get_wtime();
	int i,j;
	int begin,stop;

	for(j=0;j<chains.size();j++)
	{
		if(j==0)
			begin=0;
		else
			begin=chains.at(j-1);
		stop=chains.at(j);

		if(stop-begin<3)
			cout<<"Chain "<<j<<" only has 1 or 2 residues! I don't know how to do."<<endl;
		//cout << j << ": before" << endl;
		v[begin+0]->bb(t);
		//cout << j << ": middle" << endl;
		v[begin+1]->follow_bb(t);
		//cout << j << ": after" << endl;

		for(i=begin+1;i<stop-1;i++)
		{
			v[i]->bb(t);
			v[i-1]->previous_bb(t);
			v[i+1]->follow_bb(t);
		}
		v[stop-1]->bb(t);
	}
	cout << "pdb::getbb: " << omp_get_wtime()-st << " seconds" << endl;
	return;
}

void CPdb::getbb_assign(vector<bb_group> *t)
{
    int i,j;
    int begin,stop;
    
    for(j=0;j<chains.size();j++)
    {
        if(j==0)
            begin=0;
        else
            begin=chains.at(j-1);
        stop=chains.at(j);
        
        if(stop-begin<3)
            cout<<"Chain "<<j<<" only has 1 or 2 residues! I don't know how to do."<<endl;
        
        v[begin+0]->bb(t);
        v[begin+1]->follow_bb_assign(t);
        
        for(i=begin+1;i<stop-1;i++)
        {
            v[i]->bb(t);
            v[i-1]->previous_bb(t);
            v[i+1]->follow_bb_assign(t);
        }
        v[stop-1]->bb(t);
    }
    return;
}

void CPdb::bbhbond(vector<bbhbond_group> *t)
{	
	double st = omp_get_wtime();
	unsigned int i;
	for(i=0;i<(int)v.size();i++)
	{
		v[i]->bbhbond(t);
	}
	cout << "pdb::bbhbond: " << omp_get_wtime()-st << " seconds" << endl;
	return;
}





void CPdb::schbond(vector<bbhbond_group> *t)
{
	double st = omp_get_wtime();
	unsigned int i;
	for(i=0;i<v.size();i++)
	{
		v[i]->schbond(t);
	}
	cout << "pdb::schbond: " << omp_get_wtime()-st << " seconds" << endl;
	return;
}

void CPdb::caha(vector<index_three> *t)
{
	int i;
	for(i=0;i<(int)v.size();i++)
	{
		if(v[i]->OneLetterName != 'G')
			v[i]->caha(t);
	}
	return;
}

void CPdb::bbco(vector<co_group> *t)
{
	int i;
	for(i=0;i<(int)v.size()-1;i++)
	{
		v[i]->bbco(t);
	}
	return;
}


void CPdb::bbnh(vector<nh_group> *t)
{
	double st = omp_get_wtime();
	int i;
	for(i=0+1;i<(int)v.size();i++)
	{
		if(v[i]->OneLetterName != 'P')
			v[i]->bbnh(t);
	}
	cout << "pdb::bbnh: " << omp_get_wtime()-st << " seconds" << endl;
	return;
}


void CPdb::ired(vector<struct ired> *t)
{
	int i;
	for(i=0;i<(int)v.size();i++)
	{
		v[i]->ired(t,i);
	}
	return;
}

void CPdb::clearred(void)
{
	int i;
	for(i=0;i<v.size();i++)
		v.at(i)->clearred();
	return;
}


void CPdb::loadred(vector<struct ired> *t)
{
	int i;

	clearred();
	for(i=0;i<t->size();i++)
	{
		if(t->at(i).pos>=0 && t->at(i).pos<v.size())
			v.at(t->at(i).pos)->loadred(t->at(i));
	}

	return;
}
		
void CPdb::loadred(string filename)
{
	ifstream fin(filename.c_str());
	string p,line;
	string name1,name2;
	istringstream iss;
	vector<struct ired> red;
	struct ired t;
	int i,j,base;
	string seq1,seq2;
	vector<int> out;
	double dtotal;
	int n,total;

	clearred();

	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		iss>>t.id>>t.code>>t.s2.name1>>t.s2.name2>>t.s2.exp;
		red.push_back(t);
	}
	fin.close();


	if(red.size()>0) base=red.at(0).id-1;
	for(i=0;i<(int)red.size();i++)
		red.at(i).id-=base;

	base=0;
	seq2.clear();
	for(i=0;i<(int)red.size();i++)
	{
		if(red.at(i).id>base)
		{
			base=red.at(i).id;
			if(seq2.size()<base)
				seq2.resize(base,'X');
			seq2.at(base-1)=red.at(i).code;
		}
	}

	for(i=0;i<(int)v.size();i++)
	{
		seq1.push_back(v.at(i)->OneLetterName);
	}
	out=Sequence::align(seq1,seq2);

	dtotal=0.0;
	n=0;
	for(i=2;i<(int)out.size()-1;i++)
	{
		j=out.at(i);
		if(j!=0 && out.at(i)==out.at(i-1)+1 && out.at(i)==out.at(i+1)-1)
		{
			dtotal+=j-i;
			n++;
		}
	}
	dtotal/=n;
	if(dtotal>=0) total=(int)(dtotal+0.5);
	else total=(int)(dtotal-0.5);



	for(i=0;i<(int)red.size();i++)
	{
		red.at(i).pos=red.at(i).id-total-1;
	}

	loadred(&red);

};
	


void CPdb::getring(vector<ring_group> *t)
{	
	double st = omp_get_wtime();
	int i;

	for(i=0;i<(int)v.size();i++)
	{
		v[i]->ring(t);
	}

	for(i=t->size()-1;i>=0;i--)
	{
		if(t->at(i).x1==3 || t->at(i).x1==4)
		{
			if(t->at(i).x2<0 || t->at(i).x3<0 || t->at(i).x4<0 || t->at(i).x5<0 || t->at(i).x6<0 )
				t->erase(t->begin()+i);
		}
		else
		{
			if(t->at(i).x2<0 || t->at(i).x3<0 || t->at(i).x4<0 || t->at(i).x5<0 || t->at(i).x6<0 || t->at(i).x7<0)
				t->erase(t->begin()+i);
		}
	}
	cout << "pdb::getring: " << omp_get_wtime()-st << " seconds" << endl;
	return;
}

void CPdb::getred(int in,vector<struct ired> *t)
{
	in--;

	if(in<0 || in>v.size()-1)
		;
	else
		v.at(in)->ired(t,in);

	return;
}

char CPdb::code(int in)
{
	char c;
	in=in-1;
	if(in<0 || in >(int)v.size()-1)
		c='X';
	else
		c=v_oneletternames[in];
		//c=v[in]->OneLetterName;

	return c;
}

int CPdb::chain(int in)
{
	int r;
	in=in-1;
	if(in<0 || in>(int)v.size()-1)
		r=-1;
	else
		r=v[in]->chain;
	return r;
}

void CPdb::name(int in,char *name)
{	in=in-1;
	if(in<0 || in >(int)v.size()-1)
		strcpy(name,"err");
	else
		strcpy(name,v[in]->ThreeLetterName);
}

void CPdb::proton(vector<struct proton> *sel, int flag)
{
	int i;
	sel->clear();
	for(i=0;i<(int)v.size();i++)
	{
		v.at(i)->proton2(sel);
	}

	if(flag==1)
	{
		for(i=0;i<sel->size();i++)
		{
			bool bmiss=0;
			for(int j=0;j<sel->at(i).nh;j++)
			{
				if(sel->at(i).hpos[j]<0)
					bmiss=1;
			}
			if(bmiss==1)
			{
				#ifndef IGNORE_UNKNOWN
				cerr<<"Residue "<<sel->at(i).id<<" "<<sel->at(i).code<<" contain missing protons "<<sel->at(i).name<<" , removed"<<endl;
				#endif
				sel->erase(sel->begin()+i);
				i--;
			}
		}
	}
}


void CPdb::proton(vector<struct proton> *sel)
{
	double st = omp_get_wtime();
	int i;
	sel->clear();
	for(i=0;i<(int)v.size();i++)
	{
		v.at(i)->proton2(sel);
	}
	for(i=0;i<sel->size();i++)
	{
		bool bmiss=0;
		for(int j=0;j<sel->at(i).nh;j++)
		{
			if(sel->at(i).hpos[j]<0)
				bmiss=1;
		}
		if(bmiss==1)
		{
			#ifndef IGNORE_UNKNOWN
			cerr<<"Residue "<<sel->at(i).id<<" "<<sel->at(i).code<<" contain missing protons "<<sel->at(i).name<<" , removed"<<endl;
			#endif
			sel->erase(sel->begin()+i);
			i--;
		}
	}
	cout << "pdb::proton: " << omp_get_wtime()-st << " seconds" << endl;
}


// New function for OpenACC
// Sequential function, but is a significant performance upgrade from original
// Created new function as it functions very differently from original but produces same results
void CPdb::proton_acc(vector<struct proton> * sel)
{
	double st = omp_get_wtime();
	vector<struct proton> tmp;
	int i,j;
	bool bmiss;
	sel->clear();
	for(i=0;i<(int)v.size();i++)
	{
		v.at(i)->proton2(&tmp);
	}
	sel->reserve(tmp.size());
	for(i=0;i<tmp.size();i++)
	{
		bmiss=0;
		for(j=0;j<tmp.at(i).nh;j++)
		{
			if(tmp.at(i).hpos[j]<0){
				bmiss=1;
				break;
			}
		}
		if(bmiss!=1)
		{
			sel->push_back(tmp.at(i));
		} else {
			#ifndef IGNORE_UNKNOWN
			cerr<<"Residue "<<tmp.at(i).id<<" "<<tmp.at(i).code<<" contain missing protons "<<tmp.at(i).name<<" , removed"<<endl;
			#endif
		}
	}
	cout << "pdb::proton_nofilter: " << omp_get_wtime()-st << " seconds" << endl;
}


void CPdb::allproton3(vector<struct proton> *sel)
{
	double st = omp_get_wtime();
	int i;
	sel->clear();
	for(i=0;i<(int)(int)v.size();i++)
	{
		v.at(i)->proton3(sel);
	}
	for(i=0;i<sel->size();i++)
	{
		bool bmiss=0;
		for(int j=0;j<sel->at(i).nh;j++)
		{
			if(sel->at(i).hpos[j]<0)
				bmiss=1;
		}
		if(bmiss==1)
		{
			#ifndef IGNORE_UNKNOWN
			cerr<<"Residue "<<sel->at(i).id<<" "<<sel->at(i).code<<" contain missing protons "<<sel->at(i).name<<" , removed"<<endl;
			#endif
			sel->erase(sel->begin()+i);
			i--;
		}
	}
	cout << "pdb::allproton3: " << omp_get_wtime()-st << " seconds" << endl;
}


// New function for OpenACC
// Sequential function, but is a significant performance upgrade from original
// Created new function as it functions very differently from original but produces same results
void CPdb::allproton3_acc(vector<struct proton> *sel)
{
	double st = omp_get_wtime();
	int i,j;
	vector<struct proton> tmp;
	sel->clear();
	for(i=0;i<(int)(int)v.size();i++)
	{
		v.at(i)->proton3(&tmp);
	}
	sel->reserve(tmp.size());
	for(i=0;i<tmp.size();i++)
	{
		bool bmiss=0;
		for(j=0;j<tmp.at(i).nh;j++)
		{
			if(tmp.at(i).hpos[j]<0)
				bmiss=1;
		}
		if(bmiss!=1)
		{
			sel->push_back(tmp.at(i));
			
		} else {
			#ifndef IGNORE_UNKNOWN
			cerr<<"Residue "<<tmp.at(i).id<<" "<<tmp.at(i).code<<" contain missing protons "<<tmp.at(i).name<<" , removed"<<endl;
			#endif
		}
	}
	cout << "pdb::allproton3_acc: " << omp_get_wtime()-st << " seconds" << endl;
}


void CPdb::allproton(vector<struct proton> *sel)
{
	double st = omp_get_wtime();
	int i;
	sel->clear();
	for(i=0;i<(int)(int)v.size();i++)
	{
		v.at(i)->proton(sel);
	}
	for(i=0;i<sel->size();i++)
	{
		bool bmiss=0;
		for(int j=0;j<sel->at(i).nh;j++)
		{
			if(sel->at(i).hpos[j]<0)
				bmiss=1;
		}
		if(bmiss==1)
		{
			//cerr<<"Residue "<<sel->at(i).id<<" "<<sel->at(i).code<<" contain missing protons "<<sel->at(i).name<<" , removed"<<endl;
			sel->erase(sel->begin()+i);
			i--;
		}
	}
	cout << "pdb::allproton: " << omp_get_wtime()-st << " seconds" << endl;
}


// New function for OpenACC
// Sequential function, but is a significant performance upgrade from original
// Created new function as it functions very differently from original but produces same results
void CPdb::allproton_acc(vector<struct proton> *sel)
{
	double st = omp_get_wtime();
	int i,j;
	vector<struct proton> tmp;
	sel->clear();
	for(i=0;i<(int)(int)v.size();i++)
	{
		v.at(i)->proton(&tmp);
	}
	sel->reserve(tmp.size());
	for(i=0;i<tmp.size();i++)
	{
		bool bmiss=0;
		for(j=0;j<tmp.at(i).nh;j++)
		{
			if(tmp.at(i).hpos[j]<0)
				bmiss=1;
		}
		if(bmiss!=1)
		{
			sel->push_back(tmp.at(i));
		} else {
			#ifndef IGNORE_UNKNOWN
			cerr<<"Residue "<<tmp.at(i).id<<" "<<tmp.at(i).code<<" contain missing protons "<<tmp.at(i).name<<" , removed"<<endl;
			#endif
		}
	}
	cout << "pdb::allproton: " << omp_get_wtime()-st << " seconds" << endl;
}

void CPdb::ani(vector<struct ani_group> *anistropy)
{
	double st = omp_get_wtime();
	int i,j;
	int begin,stop;

	for(j=0;j<chains.size();j++)
	{
		if(j==0)
			begin=0;
		else
			begin=chains.at(j-1);
		stop=chains.at(j);

	
		for(i=begin;i<stop-1;i++)
		{
			v[i]->bbani(anistropy);
			v[i]->ani(anistropy);
		}
		if(stop>1) v[stop-1]->ani(anistropy);

	}


	for(i=anistropy->size()-1;i>=0;i--)
	{
		if(anistropy->at(i).pos[0]<0 || anistropy->at(i).pos[1]<0 || anistropy->at(i).pos[2]<0 )
			anistropy->erase(anistropy->begin()+i);
	}
	cout << "pdb::ani: " << omp_get_wtime()-st << " seconds" << endl;
	return;
}


// New function for OpenACC
// Sequential function, but is a significant performance upgrade from original
// Created new function as it functions very differently from original but produces same results
void CPdb::ani_acc(vector<struct ani_group> *anistropy)
{
	double st = omp_get_wtime();
	int i,j;
	vector<struct ani_group> tmp;
	
	int begin,stop;

	for(j=0;j<chains.size();j++)
	{
		if(j==0)
			begin=0;
		else
			begin=chains.at(j-1);
		stop=chains.at(j);

	
		for(i=begin;i<stop-1;i++)
		{
			v[i]->bbani(&tmp);
			v[i]->ani(&tmp);
		}
		if(stop>1) v[stop-1]->ani(&tmp);

	}
	anistropy->reserve(tmp.size());
	for(i=0; i<tmp.size(); i++)
	{
		if(tmp.at(i).pos[0]>=0 && tmp.at(i).pos[1]>=0 && tmp.at(i).pos[2]>=0)
			anistropy->push_back(tmp.at(i));
	}
	cout << "pdb::ani_acc: " << omp_get_wtime()-st << " seconds" << endl;
	return;
}

void CPdb::heavycoor()
{
	int i;
	heavy.clear();
	boundary.clear();
	allcoor.clear();
	for(i=0;i<(int)v.size();i++)
	{
		v.at(i)->bbheavycoor(&heavy);
		v.at(i)->heavycoor(&heavy);
		boundary.push_back(heavy.size());
		v.at(i)->allcoor(&allcoor);
	}

	for(i=0;i<(int)ligand.size();i++)
	{
		ligand.at(i)->heavycoor(&heavy);
		boundary.push_back(heavy.size());
	}

	return;
}



void CPdb::output(string filename)
{
	int i,j;
	int begin,stop;
	FILE *fp;

	fp=fopen(filename.c_str(),"wt");

	for(i=0;i<(int)v.size();i++)
	{
		if(i==0)
			begin=0;
		else
			begin=boundary.at(i-1);
		stop=boundary.at(i);
		for(j=begin;j<stop;j++)
		{
			if(heavy.at(j)>0)
			{
				fprintf(fp,"%d %s %s",i+1,v.at(i)->ThreeLetterName,atomname.at(heavy.at(j)-1).c_str());
				//fprintf(fp," %f %f %f",x.at(heavy.at(j)-1),y.at(heavy.at(j)-1),z.at(heavy.at(j)-1));
				fprintf(fp," %f %f\n",b.at(heavy.at(j)-1),b2.at(heavy.at(j)-1));
			}
		}
	}
	fclose(fp);
	return;
}


void CPdb::caoutput(string filename)
{
	int i,j;
	FILE *fp;

	fp=fopen(filename.c_str(),"wt");

	for(i=0;i<(int)v.size();i++)
	{
		if(i==0)
			j=1;
		else
			j=boundary.at(i-1)+1;
		if(heavy.at(j)>0)
		{
			fprintf(fp,"%d %s",i+1,v.at(i)->ThreeLetterName);
			fprintf(fp," %f %f\n",b.at(heavy.at(j)-1),b2.at(heavy.at(j)-1));
		}
	}
	fclose(fp);
	return;
}


// New function for OpenACC
// Runs sequentially but fits better with the structure of the OpenACC code
void CPdb::attach_bbprediction(int id,double ca, double cb, double co, double n, double h, double ha)
{
	v.at(id-1)->attach_bbprediction(ca, cb, co, n, h, ha);
}


void CPdb::attach_bbprediction(int id,double pre[])
{
	v.at(id-1)->attach_bbprediction(pre);
}

void CPdb::attach_protonprediction(int id,string name,double cs)
{
	v.at(id-1)->attach_protonprediction(name,cs);
}

void CPdb::print_prediction()
{
	string test="bmrb_pre.dat";
	print_prediction(test);
}


void CPdb::print_debug(string name){
	int i,j;
	FILE *fp=fopen(name.c_str(),"wt");
	char toprint0[]="loop_\n    _Residue_seq_code\n     _Residue_label\n";
	char toprint[]="      loop_\n      _Atom_shift_assign_ID\n      _Residue_seq_code\n      _Residue_label\n      _Atom_name\n      _Atom_type\n      _Chem_shift_value\n      _Chem_shift_value_error\n      _Chem_shift_ambiguity_code\n";


	fprintf(fp,toprint0);
	for(i=0;i<(int)v.size();i++)
	{
		fprintf(fp," %d %s",i+1,v.at(i)->ThreeLetterName);
		if(i%3==0) fprintf(fp,"\n");
	}
	fprintf(fp,"\n      stop_\n");
	
	
	fprintf(fp,toprint);
	j=1;
	for(i=0;i<(int)v.size();i++)
		v.at(i)->print_prediction(&j,fp);
	fprintf(fp,"      stop_\n");
	fclose(fp);
}

void CPdb::print_prediction(string name)
{
	int i,j;
	FILE *fp=fopen(name.c_str(),"wt");
	char toprint0[]="loop_\n    _Residue_seq_code\n     _Residue_label\n";
	char toprint[]="      loop_\n      _Atom_shift_assign_ID\n      _Residue_seq_code\n      _Residue_label\n      _Atom_name\n      _Atom_type\n      _Chem_shift_value\n      _Chem_shift_value_error\n      _Chem_shift_ambiguity_code\n";


	fprintf(fp,toprint0);
	for(i=0;i<(int)v.size();i++)
	{
		fprintf(fp," %d %s",i+1,v.at(i)->ThreeLetterName);
		if(i%3==0) fprintf(fp,"\n");
	}
	fprintf(fp,"\n      stop_\n");
	
	
	fprintf(fp,toprint);
	j=1;
	for(i=0;i<(int)v.size();i++)
		v.at(i)->print_prediction(&j,fp);
	fprintf(fp,"      stop_\n");
	fclose(fp);

	fp=fopen("bb_predict.dat","wt");
	fprintf(fp,"%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","index","residue","pre_ca","exp_ca",
		"pre_cb","exp_cb","pre_c","exp_c","pre_h","exp_h","pre_n","exp_n","pre_ha","exp_ha");
	for(i=0;i<(int)v.size();i++)
		v.at(i)->print_bbprediction(fp);
	fclose(fp);

	fp=fopen("proton_predict.dat","wt");
	fprintf(fp,"%8s %8s %8s %8s %8s\n","#  index","aa type","type","cs","exp");
	for(i=0;i<(int)v.size();i++)
		v.at(i)->print_protonprediction(fp);
	fclose(fp);

}

double CPdb::test_bmbr(class CBmrb bmrb)
{
	int i;
	vector<int> index;
	string seq1,seq2;
	vector<int> out;
	int total;



	index.clear();
	seq2=bmrb.getseq(index);
	for(i=0;i<(int)v.size();i++)
	{
		seq1.push_back(v.at(i)->OneLetterName);
	}
	out=Sequence::align(seq1,seq2);

	total=0;
	for(i=1;i<out.size();i++)
	{
		if(out.at(i)>0 && seq1.at(i-1)!='X')
			total++;
	}

	return ((double)total)/seq2.length();
}


int CPdb::attach_bmrb(class CBmrb bmrb)
{
	double st = omp_get_wtime();
	int i,j;
	vector<int> index;
	string seq1,seq2;
	vector<int> out;
	int total,n,m;
	double dtotal;
	char c1,c2;

	for(i=0;i<v.size();i++)
		v.at(i)->clearexp();


	index.clear();
	seq2=bmrb.getseq();
	for(i=0;i<(int)v.size();i++)
	{
		seq1.push_back(v.at(i)->OneLetterName);
	}
	out=Sequence::align(seq1,seq2);

	dtotal=0.0;
	n=0;
	for(i=2;i<(int)out.size()-1;i++)
	{
		j=out.at(i);
		if(j!=0 && out.at(i)==out.at(i-1)+1 && out.at(i)==out.at(i+1)-1)
		{
			dtotal+=j-i;
			n++;
		}
	}

	dtotal/=n;


	if(dtotal>=0) total=(int)(dtotal+0.5);
	else total=(int)(dtotal-0.5);



	n=0;
	m=0;
	for(i=1;i<(int)out.size();i++)
	{
		for(j=0;j<bmrb.getsize();j++)
		{
			if(bmrb.getdata(j).res-total==i)
			{
				c1=v.at(i-1)->OneLetterName;
				c2=Sequence::name2code(bmrb.getdata(j).name);
				if(c1==c2)
				{
					v.at(i-1)->loadexp(bmrb.getdata(j));
					n++;
				}
				else
				{
					m++;
					v.at(i-1)->set_mismatch();
				}
			}
		}
	}
	cout << "pdb::attach_bmrb: " << omp_get_wtime()-st << " seconds" << endl;
	return m;
}

void CPdb::clear_cs()
{
	int i;
	for(i=0;i<(int)v.size();i++)
	{
		v.at(i)->bexploaded=0;
		v.at(i)->exploaded=0;
		v.at(i)->clearexp();
	}
	return;
}



void CPdb::attach_rmsf(vector<double> t)
{
	int i;
	for(i=0;i<v.size();i++)
	{
		v.at(i)->attach_rmsf(t);
	}
	return;
}

void CPdb::print_rmsf(string filename)
{
	FILE *fp;
	int i;

	fp=fopen(filename.c_str(),"rt");
	for(i=0;i<v.size();i++)
	{
		v.at(i)->print_rmsf(fp);
	}
	fclose(fp);
	return;
}




//flag=1: bmrb rc
//flag==2: wang' rc
//flag==3: wang's rc with neighboring correction
void CPdb::attach_coil(int flag=1)
{
	int i;

	if(flag==3 && v.size()>=3 )
	{
		v.at(0)->set_coil_wc('X',v.at(1)->OneLetterName);
		for(i=0+1;i<v.size()-1;i++)
			v.at(i)->set_coil_wc(v.at(i-1)->OneLetterName,v.at(i+1)->OneLetterName);
		v.at(i)->set_coil_wc(v.at(i-1)->OneLetterName,'X');
	}
	else
	{
		for(i=0;i<(int)v.size();i++)
		{
			if(v.at(i)->bexploaded==0)
				v.at(i)->set_coil(flag);
		}
	}


	return;
}

void CPdb::attach_mean()
{
	int i;
	for(i=0;i<(int)v.size();i++)
	{
		if(v.at(i)->bexploaded==0)
			v.at(i)->set_mean();
	}
	return;
}

void CPdb::attach_dssp(string filename)
{
	int i,j;
	string seq1,seq2;
	vector<int> out;


	dssp.loaddata(filename);
	for(i=0;i<(int)v.size();i++)
	{
		seq1.push_back(v.at(i)->OneLetterName);
	}
	seq2=dssp.getseq();
	out=Sequence::align(seq1,seq2);
	for(i=1;i<(int)out.size();i++)
	{
		j=out.at(i);
		if(j!=0)
			v.at(i-1)->setdssp(dssp.data.at(j-1).ss);
	}
	return;
}


void CPdb::attach_dssp()
{
	int i,j;
	string seq1,seq2;
	vector<int> out;


	dssp.loaddata(pdbfilename);
	for(i=0;i<(int)v.size();i++)
	{
		seq1.push_back(v.at(i)->OneLetterName);
	}
	seq2=dssp.getseq();
	out=Sequence::align(seq1,seq2);
	for(i=1;i<(int)out.size();i++)
	{
		j=out.at(i);
		if(j!=0)
			v.at(i-1)->setdssp(dssp.data.at(j-1).ss);
	}
	return;
}

char CPdb::getss(int i)
{
	return v.at(i-1)->getdssp();
}

struct noeatoms CPdb::query(int resid,string name)
{
	return v.at(resid-1)->query(name);
}


CPdb::CPdb()
{nmiss=0;};

CPdb::~CPdb()
{
	clear();
};


//CPdbcut
CPdb2::CPdb2()
{};

CPdb2::~CPdb2()
{ clear(); };

void CPdb2::clear()
{
	int i;

	for(i=0;i<pdbs.size();i++)
	{
		delete pdbs.at(i);
	}

	pdbs.clear();
	return;
}

double CPdb2::loadpdb(string filename)
{
	int i,j;
	ifstream fin(filename.c_str());
	string name;
	string line;
	string part1,part2,part3,part4;
	struct modify mod;
	vector<struct modify> mods;
	char chain;
	char old_chain;
	vector<string> block;
	vector< vector<string> > blocks;
	class CPdb *pdb;
	int ixray;
	char buffer[255];
	double resolution;

	blocks.clear();
	clear();

	name=filename.substr(0,filename.find(".pdb"));


	old_chain=0;
	mods.clear();
	ixray=1;
	resolution=0.0;
	while(getline(fin,line))
	{
		//"SEQADV 1NBF GLZ C  3";
		//"76  UNP  P62988    G";
		//"LY    76 MODIFIED RE";
		//"SIDUE";

		part1=line.substr(0,6);
		if(part1.compare("EXPDTA")==0)
		{
			part2=line.substr(10,5);
			if(part2.compare("X-RAY")==0)
				ixray=1;
			if(part2.compare("SOLUT")==0)
				ixray=0;
			if(part2.compare("SOLID")==0)
				ixray=0;
		}

		part1=line.substr(0,22);
		if(part1.compare("REMARK   2 RESOLUTION.")==0)
		{
			part2=line.substr(22,8);
			resolution=atof(part2.c_str());
		}


		if(ixray==0)
		{
			cout<<"neglect nmr structure "<<filename.c_str()<<endl;
			break;
		}
		if(resolution>2.0)
		{
			cout<<"neglect low resolution ("<<resolution<<") xray structure "<<filename.c_str()<<endl;
			break;
		}

		part1=line.substr(0,6);
		part2=line.substr(49,20);
		if(part1.compare("SEQADV")==0 && part2.find("MODIFIED")!=string::npos)
		{
			mod.from=line.substr(12,3);
			mod.to=line.substr(39,3);
			mod.chain=line.at(16);
			mod.res=atoi(line.substr(18,4).c_str());
			mods.push_back(mod);
		}

		part1=line.substr(0,4);
		part2=line.substr(0,6);
		part3=line.substr(0,3);
		if(part1!="ATOM" && part2!="HETATM" && part3 !="TER")
			continue;
		
		chain=line.at(21);

		if(chain==old_chain)
		{
			block.push_back(line);
		}
		else if(chain>old_chain || (old_chain=='Z' && chain=='A'))
		{
			if(block.size()>0)
				blocks.push_back(block);
			old_chain=chain;
			block.clear();
			block.push_back(line);
		}
		else
		{
			//cout<<"old chain id is "<<old_chain<<" and new chain id is "<<chain<<endl;
			break;
		}
	}

	if(block.size()>0)
		blocks.push_back(block);

	ofstream  fout;
	for(i=0;i<blocks.size();i++)
	{
		block=blocks.at(i);
		sprintf(buffer,"%s-%d.pdb",name.c_str(),i);
		fout.open(buffer,ofstream::out);
		line=block.at(0);
		fout<<"REMART   1  PDB seperated from "<<filename.c_str()<<" chain ID is "<<line.at(21)<<endl;
		for(j=0;j<block.size();j++)
		{
			line=block.at(j);
			line.at(21)=' ';
			fout<<line<<endl;
		}
		fout.close();
		pdb=new CPdb;
		pdb->loadpdb(buffer);
		pdbs.push_back(pdb);
	}

	if(ixray==0)
	{ 
		resolution=100.0;
	}

	return resolution;
};


