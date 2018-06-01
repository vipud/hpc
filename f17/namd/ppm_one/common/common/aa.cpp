#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>
#include <time.h>
using namespace std;


#include "aa.h"
#include "debug.h"


void CAminoacid::methyl_ambig(int flag)
{
	return;
}

void CAminoacid::printpdb(FILE *fp,vector<string> atomname,vector<double> x,vector<double> y,vector<double> z,vector<double> b, int &n)
{
	int i;
	struct Atom t;

	for(i=0;i<(int)atoms.size();i++)
	{
		t=atoms.at(i);
		if(t.index>0)
		{
			t.index--;
			fprintf(fp,"%6s","ATOM  ");
			fprintf(fp,"%5d",n); n++;
			fprintf(fp," ");
			fprintf(fp,"%4s",atomname.at(t.index).c_str());
			fprintf(fp," ");
			fprintf(fp,"%3s",ThreeLetterName);
			fprintf(fp," "); //reserved space
			fprintf(fp," "); //chain ID
			fprintf(fp,"%4d",residue);
			fprintf(fp,"    ");
			fprintf(fp,"%8.3f%8.3f%8.3f",x.at(t.index),y.at(t.index),z.at(t.index));
			fprintf(fp,"%6.2f",1.00); //occu
			fprintf(fp,"%6.2f",b.at(t.index));
			fprintf(fp,"\n");
		}
	}
	return;
}

int CAminoacid::get_proton(struct proton *t)
{
	int i;
	int n;
	int b;
	string cname;

	t->code=OneLetterName;
	t->id=residue;
	t->multy=0;

	b=1;
	n=0;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).cs_name==t->name)
		{
			t->hpos[n]=atoms.at(i).index;
			if(atoms.at(i).index<0)
			{
				b=0;
			}
			t->exp=atoms.at(i).cs_exp;
			t->type=atoms.at(i).proton_type;
			cname=atoms.at(i).carbon_name;
			n++;

		}
	}
	t->nh=n;

	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).name==cname)
		{
			t->cpos=atoms.at(i).index;
			t->exp_c=atoms.at(i).cs_exp;
			break;
		}
	}

	return b;
}

int CAminoacid::get_proton3(struct proton *t)
{
	int i;
	int n;
	int b;
	string cname;

	t->code=OneLetterName;
	t->id=residue;
	t->exp=0.0;
	t->multy=0;


	b=1;
	n=0;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).cs_name==t->name)
		{
			t->hpos[n]=atoms.at(i).index;
			if(atoms.at(i).index<0)
			{
				b=0;
			}
			t->exp=atoms.at(i).cs_exp;
			t->type=atoms.at(i).proton_type;
			cname=atoms.at(i).carbon_name;
			n++;

		}
	}
	t->nh=n;

	
	
	if(n==1) //not a methyl group or NH2 group or aromatic ring group
	{
		t->exp=0.0;
		b=1;
		n=0;
		for(i=0;i<(int)atoms.size();i++)
		{
			if(atoms.at(i).carbon_name==t->cname )
			{
				t->hpos[n]=atoms.at(i).index;
				if(atoms.at(i).index<0)
				{
					b=0;
				
				}
				t->exp+=atoms.at(i).cs_exp;
				if(n==0)
				{
					t->name=atoms.at(i).name;
					t->exp1=atoms.at(i).cs_exp;
					t->type=atoms.at(i).proton_type;
					cname=atoms.at(i).carbon_name;
				}
				else if(n==1)
				{
					t->multy=1;
					t->name2=atoms.at(i).name;
					t->exp2=atoms.at(i).cs_exp;
					t->type2=atoms.at(i).proton_type;
				}
				n++;
			}
		}
		t->nh=n;
		t->exp/=n;
	}






	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).name==cname)
		{
			t->cpos=atoms.at(i).index;
			t->exp_c=atoms.at(i).cs_exp;
			break;
		}
	}

	return b;
}

void CAminoacid::remove_ambig(int flag)
{
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).ambig>1 && atoms.at(i).ambig!=3 )
		{
			if(atoms.at(i).cs_exp>0.0)
				atoms.at(i).cs_exp=999.0;
		}
	}
	return;
};

//
void CAminoacid::combine_hsamec(int flag)
{
	int i,j;
	string oldname,old_basename;
	double cs;

	oldname=atoms.at(0).carbon_name;
	old_basename=atoms.at(0).base_name;
	j=0;

	for(i=1;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).proton==1 && atoms.at(i).bb==0 
			&&  atoms.at(j).proton==1 && atoms.at(j).bb==0 
			&&  atoms.at(i).carbon_name==oldname && atoms.at(i).base_name!=old_basename)
		{
			//mix i and j here
			cs=(atoms.at(i).cs_exp+atoms.at(j).cs_exp)/2.0;
			atoms.at(i).cs_exp=cs;
			atoms.at(j).cs_exp=cs;
		}
		oldname=atoms.at(i).carbon_name;
		old_basename=atoms.at(i).base_name;
		j=i;
	}

	return;
};


void CAminoacid::heavycoor(vector<int> *t)
{
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).bb==0 && atoms.at(i).proton==0)
			t->push_back(atoms.at(i).index);
	}
	return;
}

void CAminoacid::bbheavycoor(vector<int> *t)
{
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).bb==1 && atoms.at(i).proton==0 && atoms.at(i).name.compare("OXT")!=0)
			t->push_back(atoms.at(i).index);
	}
	return;
}


void CAminoacid::allcoor(vector<int> *t)
{
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).index>0)
			t->push_back(atoms.at(i).index);
	}
	return;
}

struct noeatoms CAminoacid::query(string name)
{
	struct noeatoms t;
	vector<int> tt;
	int index;

	if(name=="HN")
		name="H";

	tt.clear();
	index=get(name.c_str()).index;
	if(index<=0 && name=="H")
		index=get("H1").index; 
	if(index>0)
	{
		tt.push_back(index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	return t;
}


// New function for OpenACC
// Runs sequentially, but works better with how the OpenACC code is laid out
void CAminoacid::attach_bbprediction(double pre_ca, double pre_cb, double pre_c, double pre_n, double pre_h, double pre_ha)
{
	int i;

	if(OneLetterName=='G')
		pre_cb=999.0;
	else if(OneLetterName=='C')
		pre_ca=pre_cb=pre_c=pre_h=pre_n=999.0;
	else if(OneLetterName=='P')
		pre_h=pre_n=999.0;
	else if(OneLetterName=='U')
		pre_ca=pre_cb=pre_c=pre_h=pre_n=999.0;


	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).name=="CA")
			atoms.at(i).cs_pre=pre_ca;
		if(atoms.at(i).name=="CB")
			atoms.at(i).cs_pre=pre_cb;
		if(atoms.at(i).name=="C")
			atoms.at(i).cs_pre=pre_c;
		if(atoms.at(i).name=="H")
			atoms.at(i).cs_pre=pre_h;
		if(atoms.at(i).name=="N")
			atoms.at(i).cs_pre=pre_n;
		if(atoms.at(i).name=="HA" || atoms.at(i).name=="HA2" || atoms.at(i).name=="HA3")
			atoms.at(i).cs_pre=pre_ha;
	}

}



void CAminoacid::attach_bbprediction(double pre[5])
{
	int i;

	pre_ca=pre[0];
	pre_cb=pre[1];
	pre_c=pre[2];
	pre_h=pre[3];
	pre_n=pre[4];
	pre_ha=pre[5];

	if(OneLetterName=='G')
		pre_cb=999.0;
	else if(OneLetterName=='C')
		pre_ca=pre_cb=pre_c=pre_h=pre_n=999.0;
	else if(OneLetterName=='P')
		pre_h=pre_n=999.0;
	else if(OneLetterName=='U')
		pre_ca=pre_cb=pre_c=pre_h=pre_n=999.0;


	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).name=="CA")
			atoms.at(i).cs_pre=pre_ca;
		if(atoms.at(i).name=="CB")
			atoms.at(i).cs_pre=pre_cb;
		if(atoms.at(i).name=="C")
			atoms.at(i).cs_pre=pre_c;
		if(atoms.at(i).name=="H")
			atoms.at(i).cs_pre=pre_h;
		if(atoms.at(i).name=="N")
			atoms.at(i).cs_pre=pre_n;
		if(atoms.at(i).name=="HA" || atoms.at(i).name=="HA2" || atoms.at(i).name=="HA3")
			atoms.at(i).cs_pre=pre_ha;
	}

}

void  CAminoacid::attach_protonprediction(string name,double cs)
{	
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).cs_name==name)
		{
			atoms.at(i).cs_pre=cs;
		}
	}
	return;
}



void CAminoacid::print_prediction(int *index,FILE *fbmrb)
{
	int i;
	string oldname;

	oldname="xxx";
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).cs_pre<450 && atoms.at(i).cs_name!=oldname)
		{
			fprintf(fbmrb,"%8d %8d %8s %8s %8s %8.3f . 1\n",*index,original_residue,ThreeLetterName,atoms.at(i).cs_name.c_str(),
				atoms.at(i).cs_name.c_str(),atoms.at(i).cs_pre);
			(*index)++;
			oldname=atoms.at(i).cs_name;
		}
	}
	return;
}

void CAminoacid::print_bbprediction(FILE *fbmrb)
{
	fprintf(fbmrb,"%8d%8s",original_residue,ThreeLetterName);
	fprintf(fbmrb,"%8.3f%8.3f",get("CA").cs_pre,get("CA").cs_exp);
	fprintf(fbmrb,"%8.3f%8.3f",get("CB").cs_pre,get("CB").cs_exp);
	fprintf(fbmrb,"%8.3f%8.3f",get("C").cs_pre,get("C").cs_exp);
	fprintf(fbmrb,"%8.3f%8.3f",get("H").cs_pre,get("H").cs_exp);
	fprintf(fbmrb,"%8.3f%8.3f",get("N").cs_pre,get("N").cs_exp);
	if(OneLetterName=='G')
		fprintf(fbmrb,"%8.3f%8.3f",get("HA2").cs_pre,get("HA2").cs_exp);
	else
		fprintf(fbmrb,"%8.3f%8.3f",get("HA").cs_pre,get("HA").cs_exp);
	fprintf(fbmrb,"\n");
}

void CAminoacid::print_protonprediction(FILE *fbmrb)
{
	int i;
	string oldname;

	oldname="xxx";
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).cs_pre<990 && atoms.at(i).proton==1 && atoms.at(i).bb==0 && atoms.at(i).cs_name!=oldname)
		{
			oldname=atoms.at(i).cs_name;
			fprintf(fbmrb,"%8d %8c %8s %8.3f %8.3f\n",original_residue,OneLetterName,atoms.at(i).cs_name.c_str(),
				atoms.at(i).cs_pre,atoms.at(i).cs_exp);
		}
	}
	return;
}


void CAminoacid::attach_rmsf(vector<double> t)
{
	int i,j;

	for(i=0;i<atoms.size();i++)
	{
		j=atoms.at(i).index;
		if(j>0)
			atoms.at(i).rmsf=t.at(j-1);
	}
	return;
}

void CAminoacid::print_rmsf(FILE *fp)
{
	int i;

	for(i=0;i<atoms.size();i++)
	{
		if(atoms.at(i).rmsf>0)
			fprintf(fp,"%10d %9c %9s %10.4f\n",residue,OneLetterName,atoms.at(i).name.c_str(),atoms.at(i).rmsf);
	}
	return;
}


int CAminoacid::wang_correct_index(char c)
{
	int index;
	switch(c)
		{
			case 'A': index=0; break;
			case 'C': index=1; break;
			case 'D': index=2; break;
			case 'E': index=3; break;
			case 'F': index=4; break;
			case 'G': index=5; break;
			case 'H': index=6; break;
			case 'I': index=7; break;
			case 'K': index=8; break;
			case 'L': index=9; break;
			case 'M': index=10; break;
			case 'N': index=11; break;
			case 'P': index=12; break;
			case 'Q': index=13; break;
			case 'R': index=14; break;
			case 'S': index=15; break;
			case 'T': index=16; break;
			case 'V': index=17; break;
			case 'W': index=18; break;
			case 'Y': index=19; break;
			default: index=20;
		}

	return index;
}

void CAminoacid::set_coil_wc(char pre, char fol)
{
	int index,index1,index2;
	int i;
	double t[6];

	set_coil(2);
	index1=wang_correct_index(pre);
	index2=wang_correct_index(fol);
	index=wang_rc_index(OneLetterName);

	for(i=0;i<6;i++)
	{
		t[i]=wang_rc[index][i]+wang_c1[index1][i]+wang_c2[index2][i];
	}
	set_wang(t,index);
	return;
}

int CAminoacid::wang_rc_index(char c)
{
	int index;
	switch(c)
		{
			case 'I': index=0; break;
			case 'V': index=1; break;
			case 'D': index=2; break;
			case 'N': index=3; break;
			case 'F': index=4; break;
			case 'H': index=5; break;
			case 'W': index=6; break;
			case 'Y': index=7; break;
			case 'K': index=8; break;
			case 'L': index=9; break;
			case 'M': index=10; break;
			case 'Q': index=11; break;
			case 'R': index=12; break;
			case 'E': index=13; break;
			case 'T': index=14; break;
			case 'C': index=15; break;
			case 'S': index=16; break;
			case 'A': index=17; break;
			case 'G': index=18; break;
			case 'P': index=19; break;
			default: index=20;
		}

	return index;
}

void CAminoacid::set_wang(double *t,int index)
{
	int i;
	//unknown or missing residue. 
	if(index==20)
		return;


	get_address("C")->cs_wang=t[1];
	get_address("CA")->cs_wang=t[2];

	if(index!=19)
	{
		get_address("H")->cs_wang=t[4];
		get_address("N")->cs_wang=t[0];
	}

	if(index!=18) 
	{
		get_address("HA")->cs_wang=t[5];
		get_address("CB")->cs_wang=t[3];
	}
	else
	{
		get_address("HA1")->cs_wang=t[5];
		get_address("HA2")->cs_wang=t[5];
	}

	for(i=0;i<(int)atoms.size();i++)
	{
		atoms.at(i).cs_exp=atoms.at(i).cs_wang;
	}

	return;
}

void CAminoacid::set_coil(int flag=1)
{
	int i;
	if(flag==1)
	{
		for(i=0;i<(int)atoms.size();i++)
		{
			atoms.at(i).cs_exp=atoms.at(i).cs_coil;
		}
	}
	else if(flag==2)
	{
		int index=wang_rc_index(OneLetterName);
		set_wang(wang_rc[index],index);	
	}

	return;
}

void CAminoacid::set_mean()
{
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		atoms.at(i).cs_exp=atoms.at(i).cs_mean;
	}
	return;
}

void CAminoacid::clearexp()
{
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).cs_exp>0.0)
			atoms.at(i).cs_exp=999.0;
	}
	bexploaded=0;
	exploaded=0;
	return;
}



void CAminoacid::loadexp(struct CBmrbdata data)
{
	int i,j;
	bool bused;
	string name;
	string dataname;

	if(data.name.compare(ThreeLetterName)!=0)
	{
		cout<<"Error. at pos "<<residue<<" PDB is "<<ThreeLetterName<<" while BMRB is "<<data.name<<endl;
		bexploaded=0;
		exploaded=-1;
		for(i=0;i<(int)atoms.size();i++)
		{
			if(atoms.at(i).cs_exp>0.0)
				atoms.at(i).cs_exp=9999.0;
		}
		return;
	}

	for(i=0;i<(int)data.block.size();i++)
	{
		data.block.at(i).used=1;
		bused=0;
		dataname=data.block.at(i).type;

		for(j=0;j<(int)atoms.size();j++)
		{
			name=atoms.at(j).cs_name;
			if(dataname.at(0)=='H' && dataname.size()>1)
			{
				if(data.block.at(i).type.find(name)==0 && name!="H" )
				{
					bused=1;
					if(data.block.at(i).ambig<=3)
					{
						atoms.at(j).cs_exp=data.block.at(i).cs;
						atoms.at(j).ambig=data.block.at(i).ambig;
						data.block.at(i).used=1;
					}
				}
			}

			else if(dataname.at(0)=='C' && dataname.size()>2 && name.size()>1)
			{
				if(data.block.at(i).type.find(name)==0 )
				{
					bused=1;
					if(data.block.at(i).ambig<=3)
					{
						atoms.at(j).cs_exp=data.block.at(i).cs;
						atoms.at(j).ambig=data.block.at(i).ambig;
						data.block.at(i).used=1;
					}
				}
			}

			else if(dataname.at(0)=='C' && dataname.size()>1)
			{
				if(data.block.at(i).type.compare(name)==0 )
				{
					bused=1;
					if(data.block.at(i).ambig<=3)
					{
						atoms.at(j).cs_exp=data.block.at(i).cs;
						atoms.at(j).ambig=data.block.at(i).ambig;
						data.block.at(i).used=1;
					}
				}
			}

			else
			{
				if(data.block.at(i).type.compare(name)==0 )
				{
					bused=1;
					if(data.block.at(i).ambig<=3)
					{
						atoms.at(j).cs_exp=data.block.at(i).cs;
						atoms.at(j).ambig=data.block.at(i).ambig;
						data.block.at(i).used=1;
					}
				}
			}
		}
		if(bused==0)
			cout<<"In "<<residue<<" "<<ThreeLetterName<<", unmatched bmrb data point for "<<data.block.at(i).type.c_str()<<" cs is "<<data.block.at(i).cs<<endl;
	}
	bexploaded=1;
	exploaded=1;
	return;
}

void CAminoacid::set_mismatch(void)
{
	int i;

	bexploaded=0;
	exploaded=-1;

	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).cs_exp>0.0)
			atoms.at(i).cs_exp=9999.0;
	}

	return;
}


void CAminoacid::dihe(vector<dihe_group> * dihe_index)
{
	cout<<"run virtual fuction of base class, sth is wrong!"<<endl;
	return;
}

void CAminoacid::proton(vector<struct proton> *sel)
{
	struct proton t;
	int i;
	string oldname;

	oldname="XXX";
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).bb==0 && atoms.at(i).proton==1 && atoms.at(i).cs_name!=oldname)
		{
			oldname=atoms.at(i).cs_name;
			t.name=atoms.at(i).cs_name;
			if(get_proton(&t))
				sel->push_back(t);
		}
	}
	return;
}


//combine hb2 and hb3 into one.
void CAminoacid::proton3(vector<struct proton> *sel)
{
	struct proton t;
	int i;
	string oldname,cname;

	oldname="XXX";
	cname="YYYY";
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).bb==0 && atoms.at(i).proton==1 && atoms.at(i).cs_name!=oldname && atoms.at(i).carbon_name!=cname)
		{
			oldname=atoms.at(i).cs_name;
			cname=atoms.at(i).carbon_name;
			t.name=atoms.at(i).cs_name;
			t.cname=atoms.at(i).carbon_name;

			t.exp=t.exp1=t.exp2=0.0;
			t.cname2="";
			t.name2="";
			t.cpos=-1;
			t.nh=0;
			t.multy=0;
			t.type=t.type2=0;

			if(get_proton3(&t))
				sel->push_back(t);
		}
	}
	return;
}

void CAminoacid::proton2(vector<struct proton> *)
{
	return;
}

void CAminoacid::ring(vector<ring_group> *ring)
{
	return;
}

void CAminoacid::previous_bb(vector<bb_group> *bb)
{
	if(exploaded==-1)
		bb->at(bb->size()-1).previous_mut=1;
}


void CAminoacid::follow_bb(vector<bb_group> *bb)
{
	if(exploaded==-1)
		bb->at(bb->size()-1).follow_mut=1;

	bb->at(bb->size()-1).follow_hpos=get("H").index;
	bb->at(bb->size()-1).follow_npos=get("N").index;
	bb->at(bb->size()-1).follow_exp_n=get("N").cs_exp;
	bb->at(bb->size()-1).follow_exp_h=get("H").cs_exp;
	bb->at(bb->size()-1).follow_exp_ca=get("CA").cs_exp;
	
	if(OneLetterName=='P')
	{
		bb->at(bb->size()-1).follow_exp_n=999.0;
		bb->at(bb->size()-1).follow_exp_h=999.0;
	}

}

void CAminoacid::follow_bb_assign(vector<bb_group> *bb)
{
    if(exploaded==-1)
        bb->at(bb->size()-1).follow_mut=1;
    
    bb->at(bb->size()-1).follow_hpos=get("H").index;
    bb->at(bb->size()-1).follow_npos=get("N").index;
    bb->at(bb->size()-1).follow_exp_n=get("N").cs_exp;
    bb->at(bb->size()-1).follow_exp_h=get("H").cs_exp;
    bb->at(bb->size()-1).follow_exp_ca=get("CA").cs_exp;
    
    if(OneLetterName=='P')   //Proline won't have peaks in 3D experiments
    {
        bb->at(bb->size()-1).follow_exp_n=999.0;
        bb->at(bb->size()-1).follow_exp_h=999.0;
        
    
         bb->at(bb->size()-1).follow_exp_ca=999.0;
         bb->at(bb->size()-1).exp_ca=999.0;
         bb->at(bb->size()-1).exp_cb=999.0;
         bb->at(bb->size()-1).exp_co=999.0;
    
    }
    
}


void CAminoacid::bb(vector<bb_group> *bb)
{
	struct bb_group t;

	t.chain=chain;
	t.id=residue;
	t.id0=original_residue;
	t.code=OneLetterName;
	t.ss=ss;
	
	t.exploaded=1;
	t.previous_mut=0;
	t.follow_mut=0;

	if(exploaded!=1)
		t.exploaded=0;
	
	t.hpos=get("H").index;
	t.npos=get("N").index;
	t.capos=get("CA").index;
	t.cbpos=get("CB").index;
	t.copos=get("C").index;
	t.opos=get("O").index;
	

	t.exp_ca=get("CA").cs_exp;
	t.exp_cb=get("CB").cs_exp;
	t.exp_co=get("C").cs_exp;
	t.exp_n=get("N").cs_exp;
	t.exp_h=get("H").cs_exp;
	

	if(OneLetterName=='G')
	{
		t.exp_ha=get("HA2").cs_exp;
		t.hapos=get("HA2").index;
		t.exp_ha2=get("HA3").cs_exp;
		t.hapos2=get("HA3").index;
	}
	else
	{
		t.exp_ha=get("HA").cs_exp;
		t.hapos=get("HA").index;
	}


	if(t.cbpos==-2)
		t.exp_cb=-999.0; //GLY
	if(t.hpos==-2)
		t.exp_h=-999.0; //PRO

	t.pre_ca=get("CA").cs_pre;
	t.pre_cb=get("CB").cs_pre;
	t.pre_c=get("C").cs_pre;
	t.pre_n=get("N").cs_pre;
	t.pre_h=get("H").cs_pre;
	

	if(OneLetterName=='G')
	{
		t.pre_ha=get("HA2").cs_pre;
	}
	else
	{
		t.pre_ha=get("HA").cs_pre;
	}

	bb->push_back(t);
	return;
}

void CAminoacid::bbco(vector<struct co_group> *co)
{
	struct co_group t;

	t.id=residue;
	t.code=OneLetterName;
	t.cpos=get("C").index;
	t.opos=get("O").index;
	t.exp_c=get("C").cs_exp;
	t.exp_o=get("O").cs_exp;

	co->push_back(t);
	
	return;
}


void CAminoacid::ired(vector<struct ired> *t, int pos)
{
	int i;
	string s1,s2;
	struct ired red;

	red.id=residue;
	red.id0=original_residue;
	red.chain=chain;
	red.pos=pos;
	red.code=OneLetterName;

	for(i=0;i<(int)order_parameters.size();i++)
	{
		if(order_parameters.at(i).name1=="C0")
			red.index1=previousc;
		else
			red.index1=get(order_parameters.at(i).name1.c_str()).index;
		red.index2=get(order_parameters.at(i).name2.c_str()).index;

		red.s2=order_parameters.at(i);

		if(red.index1>=0 && red.index2>=0)
			t->push_back(red);
	}



	return;
}


void CAminoacid::clearred()
{
	int i;
	for(i=0;i<(int)order_parameters.size();i++)
	{	
		order_parameters.at(i).pre=0.0;
		order_parameters.at(i).exp=0.0;	
	}
	return;
}



void CAminoacid::loadred(struct ired red)
{
	int i;
	struct S2 t;


	{
		for(i=0;i<(int)order_parameters.size();i++)
		{
			t=order_parameters.at(i);
			if( (t.name1==red.s2.name1 && t.name2==red.s2.name2) || (t.name1==red.s2.name2 && t.name2==red.s2.name1) )
			{
				order_parameters.at(i).pre=red.s2.pre;
				order_parameters.at(i).exp=red.s2.exp;
			}
		}
	}

	return;
}



void CAminoacid::bbnh(vector<struct nh_group> *nh)
{
	struct nh_group t;

	t.id=residue;
	t.code=OneLetterName;
	t.hpos=get("H").index;
	t.npos=get("N").index;
	t.exp_h=get("H").cs_exp;
	t.exp_n=get("N").cs_exp;

	if(OneLetterName!='P')
		nh->push_back(t);
	
	return;
}

void CAminoacid::caha(vector<struct index_three> *caha)
{
	struct index_three t;
	t.x1=residue;
	t.x2=get("CA").index;
	t.x3=get("HA").index;
	if(OneLetterName!='G')
		caha->push_back(t);
	return;
}

void CAminoacid::bbhbond(vector<bbhbond_group> *bb)
{
	struct bbhbond_group t;

	t.id=residue;
	t.code=OneLetterName;
	t.npos=get("N").index;
	t.hpos=get("H").index;
	t.cpos=get("C").index;
	t.opos=get("O").index;
	t.type=1;

	bb->push_back(t);
	return;
}

void CAminoacid::schbond(vector<bbhbond_group> *)
{
	//remaining AA don't have sc donor or acceptor group.
	return;
}


void CAminoacid::bbani(vector<struct ani_group> *anistropy)
{
	struct ani_group t;
	t.type=1;
	t.id=residue;
	t.code=OneLetterName;
	t.pos[0]=followingn;
	t.pos[1]=get("C").index;
	t.pos[2]=get("O").index;

	anistropy->push_back(t);
}

void CAminoacid::ani(vector<struct ani_group> *anistropy)
{
	return;
}

void CAminoacid::sccoor(vector<int> *t)
{
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).bb==0 && atoms.at(i).proton==0)
			t->push_back(atoms.at(i).index);
	}
	return;
}


void CAminoacid::setcterminal(void)
{
	struct Atom t;
	
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="OXT";t.cs_name="OXT";t.base_name="OXT";atoms.push_back(t);

	return;
}


void CAminoacid::setnterminal(void)
{
	struct Atom t;
	int i;

	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).name=="H")
			atoms.erase(atoms.begin()+i);
	}

	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;
	t.bb=1;
	t.proton=1;
	t.ambig=0;
	t.name="H1";t.cs_name="H";t.base_name="H";atoms.push_back(t);
	t.name="H2";t.cs_name="H";t.base_name="H";atoms.push_back(t);
	t.name="H3";t.cs_name="H";t.base_name="H";atoms.push_back(t);
}


void CAminoacid::process(vector<string> block)
{
	int i,j,k;
	int index;
	char c;
	bool bmatch;
	string t,t2;
	string atomname,atomname2;
	string part;
	string name;
	vector<string> subblock;

	original_residue=atoi(block.at(0).substr(22,4).c_str());


	vector<int> toremove;
	toremove.resize(block.size(),0);
	for(i=0;i<(int)block.size();i++)
	{
		for(j=i+1;j<(int)block.size();j++)
		{
			atomname=block.at(i).substr(11,5);
			atomname2=block.at(j).substr(11,5);
			if(atomname2==atomname)
			{
				toremove.at(j)=1;
			}
		}
	}

	for(i=(int)toremove.size()-1;i>0;i--)
	{
		if(toremove.at(i)==1)
			block.erase(block.begin()+i);
	}
	toremove.clear();



	
	//replace D to H
	for(i=0;i<(int)block.size();i++)
	{
		t=block[i];
		atomname=t.substr(11,3);
		if(atomname.compare("  D")==0)
		{
			t2=t.substr(0,11);
			t2=t2.append("  H");
			t2=t2.append(t.substr(14));
			//cout<<"Atomname is: "<<atomname<<endl;
			//cout<<"from: "<<t<<endl;
			//cout<<"to:   "<<t2<<endl;
			block.at(i)=t2;
		}
	}

	for(i=0;i<(int)block.size();i++)
	{
		t=block[i];
		atomname=t.substr(11,2);
		if(atomname.compare(" D")==0)
		{
			t2=t.substr(0,11);
			t2=t2.append(" H");
			t2=t2.append(t.substr(13));
			//cout<<"Atomname is: "<<atomname<<endl;
			//cout<<"from: "<<t<<endl;
			//cout<<"to:   "<<t2<<endl;
			block.at(i)=t2;
		}
	}





	for(i=0;i<(int)block.size();i++)
	{
		t=block[i];
		atomname=t.substr(11,5);  //atom name here
		while((j=atomname.find(" ")) != string::npos)
		{
			atomname.replace(j, 1, "");
		}	
		part=t.substr(6,5); //atom index
		index=atoi(part.c_str());

		if(atomname=="OC1")  //C-terminal atoms 
			atomname="O";
		if(atomname=="OC2")
			atomname="OXT";

		if(atomname=="HN")
			atomname="H";


		c=atomname.at(0);
		if(c>'0' && c<'9')  //if name started with number, move that number to the end
		{
			atomname.erase(0,1);
			atomname.append(1,c);
		}

		
		if(OneLetterName=='I')   //gromacs use CD in ILE, should be CD1
		{
			if( atomname=="CD") atomname="CD1"; 	
			if( atomname=="HD1")atomname="HD11"; 	
			if( atomname=="HD2")atomname="HD12"; 	
			if( atomname=="HD3")atomname="HD13"; 	
		}

		c=atomname.at(0);
		if(c=='H' && atomname.compare("H")!=0 && atomname.compare("H1")!=0 && atomname.compare("H2")!=0 && atomname.compare("H3")!=0 )  //started with H, but it is not HN,h1,h2,or h3
		{
			//process protons
			bmatch=0;
			for(j=0;j<(int)atoms.size()  && bmatch==0 ;j++)
			{
				name=atoms.at(j).base_name;
				k=atoms.at(j).index;
				if(name.compare("H")!=0 && atomname.find(name)==0 && k<0)
				{
					atoms.at(j).index=index;
					bmatch=1;
				}
			}
			if(bmatch==0){  //cannot match atomname, print out error message.
				#ifndef IGNORE_UNKNOWN
				cout<<"Unknown atom name "<<atomname.c_str()<<" in residue "<<residue<<" "<<ThreeLetterName<<endl;
				#endif
			}
		}
		else  //heavy atoms or HN atom. try to match name directly.
		{
			bmatch=0;
			for(j=0;j<(int)atoms.size();j++)
			{
				if(atomname==atoms.at(j).name)
				{
					atoms.at(j).index=index;
					bmatch=1;
				}
			}
			if(bmatch==0) {  //cannot match atomname, print out error message.
			#ifndef IGNORE_UNKNOWN
			cout<<"Unknown atom name "<<atomname.c_str()<<" in residue "<<residue<<" "<<ThreeLetterName<<endl;
			#endif
			}
		}
	}
	return;
}


void CAminoacid::bbdihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;

	t.id=residue;
	t.code=OneLetterName;
	t.type=1;
	t.x1=previousc;
	t.x2=get("N").index;
	t.x3=get("CA").index;
	t.x4=get("C").index;
	t.bgood=0;

	if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0)
		t.bgood=1;
	dihe_index->push_back(t);

	t.id=residue;
	t.code=OneLetterName;
	t.type=2;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("C").index;
	t.x4=followingn;
	t.bgood=0;
	if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0)
		t.bgood=1;
	dihe_index->push_back(t);

	
	return;
}

struct Atom  CAminoacid::get(const char* name)
{
	int i;
	struct Atom t;

	t.index=-2;
	t.cs_exp=-999.0;
	t.cs_pre=-999.0;
	t.name="NO_EXIST";

	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).name==name)
		{
			t=atoms.at(i);
			break;
		}
	}

	if(i==atoms.size())  //didn't find it, so ...
	{
		t.index=-2;
		t.cs_exp=-999.0;
		t.cs_pre=-999.0;
		t.name="NO_EXIST";
	}
	return t;
}

struct Atom*  CAminoacid::get_address(const char* name)
{
	int i;
	struct Atom *t;


	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).name==name)
		{
			t=&(atoms.at(i));
			break;
		}
	}

	if(i==atoms.size())  //didn't find it, so ...
	{
		t=&atom_nouse;
	}
	return t;
}




CAminoacid::CAminoacid() 
{
	OneLetterName='B';
	bexploaded=0;
	exploaded=0;
	ss=0;
	
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=1;
	t.proton=0;
	t.ambig=0;
	t.proton=0;t.name="N";t.cs_name="N";t.base_name="N";atoms.push_back(t);
	t.proton=1;t.name="H";t.cs_name="H";t.base_name="H";atoms.push_back(t);
	t.proton=0;t.name="CA";t.cs_name="CA";t.base_name="CA";atoms.push_back(t);
	t.proton=1;t.name="HA";t.cs_name="HA";t.base_name="HA";t.carbon_name="CA";atoms.push_back(t);
	t.proton=0;t.name="C";t.cs_name="C";t.base_name="C";atoms.push_back(t);
	t.proton=0;t.name="O";t.cs_name="O";t.base_name="O";atoms.push_back(t);

	struct S2 s2;
	s2.name1="N";s2.name2="H";s2.exp=999.0;s2.pre=999.0;order_parameters.push_back(s2);
	s2.name1="CA";s2.name2="HA";s2.exp=999.0;s2.pre=999.0;order_parameters.push_back(s2);
	s2.name1="N";s2.name2="CA";s2.exp=999.0;s2.pre=999.0;order_parameters.push_back(s2);
	s2.name1="CA";s2.name2="C";s2.exp=999.0;s2.pre=999.0;order_parameters.push_back(s2);
	s2.name1="C0";s2.name2="N";s2.exp=999.0;s2.pre=999.0;order_parameters.push_back(s2);

	for(int i=0;i<6;i++)
	{
		coil_pre[i]=coil_fol[i]=0.0;
	}

}

int CAminoacid::wishart_index()
{
	int index;

	switch(OneLetterName)
	{
		case 'A': index=0;break;
		case 'R': index=1;break;
		case 'N': index=2;break;
		case 'D': index=3;break;
		case 'Q': index=4;break;
		case 'E': index=5;break;
		case 'G': index=6;break;
		case 'H': index=7;break;
		case 'I': index=8;break;
		case 'L': index=9;break;
		case 'K': index=10;break;
		case 'M': index=11;break;
		case 'F': index=12;break;
		case 'P': index=13;break;
		case 'S': index=14;break;
		case 'T': index=15;break;
		case 'W': index=16;break;
		case 'Y': index=17;break;
		case 'V': index=18;break;
		case 'C': index=20;break;
		default: index=21;
	}
	if(strcmp(ThreeLetterName,"XYX")==0) index=19;

	return index;
}



vector<double> CAminoacid::get_wishart()
{
	vector <double> t;
	for(int i=0;i<12;i++)
		t.push_back(wishart[wishart_index()][i]);
	return t;
}


CAminoacid::~CAminoacid() {}



//ALA

struct noeatoms CAla::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("MB")==0 || name.compare("QB")==0)
	{
		tt.push_back(get("HB1").index);
		tt.push_back(get("HB2").index);
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}


void CAla::dihe(vector<dihe_group> * dihe_index)
{
	return;
}


void CAla::proton2(vector<struct proton> *sel)
{
	struct proton t;

	t.type=1;
	t.name="HB";
	get_proton(&t);
	sel->push_back(t);

	return;
}


CAla::CAla() 
{
	//CAminoacid::CAminoacid();
	OneLetterName='A';	strcpy(ThreeLetterName,"ALA");

	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.proton=1;t.carbon_name="CB";
	t.proton_type=1;t.name="HB1";t.cs_name="HB";t.base_name="HB";atoms.push_back(t);
	t.proton_type=1;t.name="HB2";t.cs_name="HB";t.base_name="HB";atoms.push_back(t);
	t.proton_type=1;t.name="HB3";t.cs_name="HB";t.base_name="HB";atoms.push_back(t);

	get_address("H")->cs_coil=8.24;
	get_address("N")->cs_coil=124.2;
	get_address("CA")->cs_coil=52.1;
	get_address("HA")->cs_coil=4.39;
	get_address("CB")->cs_coil=19.3;
	get_address("C")->cs_coil=176.8;

	get_address("H")->cs_mean=8.19;
	get_address("N")->cs_mean=123.22;
	get_address("CA")->cs_mean=53.19;
	get_address("HA")->cs_mean=4.25;
	get_address("CB")->cs_mean=18.98;
	get_address("C")->cs_mean=177.79;

	get_address("N")->cs_wang=123.82;
	get_address("C")->cs_wang=177.28;
	get_address("CA")->cs_wang=52.46;
	get_address("CB")->cs_wang=18.98;
	get_address("H")->cs_wang=8.09;
	get_address("HA")->cs_wang=4.31;

	coil_pre[0]=-2.21;
	coil_pre[1]= 0.14;
	coil_pre[2]=-0.01;
	coil_pre[3]=-0.04;
	coil_pre[4]=-0.07;
	coil_pre[5]=-0.01;

	coil_fol[0]=-0.11;
	coil_fol[1]= 0.05;
	coil_fol[2]= 0.07;
	coil_fol[3]=-0.09;
	coil_fol[4]=-0.01;
	coil_fol[5]=-0.03;
}
CAla::~CAla() {}


//ARG
void CArg::schbond(vector<bbhbond_group> *grp)
{
	bbhbond_group t;
	t.id=residue;
	t.code=OneLetterName;
	t.cpos=-1;
	t.opos=-1;
	t.npos=get("NE").index;
	t.hpos=get("HE").index;
	t.type=13; 

	grp->push_back(t);
	return;
}

struct noeatoms CArg::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QG")==0)
	{
		tt.clear();
		tt.push_back(get("HG2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HG3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QD")==0)
	{
		tt.clear();
		tt.push_back(get("HD2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HD3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QH1")==0)
	{
		tt.clear();
		tt.push_back(get("HH11").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HH12").index);
		t.atoms.push_back(tt);
	}
	else if(name.compare("QH2")==0)
	{
		tt.clear();
		tt.push_back(get("HH21").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HH22").index);
		t.atoms.push_back(tt);
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}


void CArg::ani(vector<struct ani_group> *anistropy)
{
	struct ani_group t;
	CAminoacid::ani(anistropy);

	t.type=4;
	t.id=residue;
	t.code=OneLetterName;
	t.pos[0]=get("NH1").index;
	t.pos[1]=get("CZ").index;
	t.pos[2]=get("NH2").index;
	anistropy->push_back(t);
	
	return ;
}



void CArg::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;

	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;;
	t.x3=get("CB").index;;
	t.x4=get("CG").index;;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("CD").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);
	
	t.type=5;
	t.x1=get("CB").index;
	t.x2=get("CG").index;
	t.x3=get("CD").index;
	t.x4=get("NE").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);
	
	t.type=6;
	t.x1=get("CG").index;
	t.x2=get("CD").index;
	t.x3=get("NE").index;
	t.x4=get("CZ").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);
	
	t.type=7;
	t.x1=get("CD").index;
	t.x2=get("NE").index;
	t.x3=get("CZ").index;
	t.x4=get("NH2").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);
	
	return;
}


CArg::CArg() 
{
	//CAminoacid::CAminoacid();
	OneLetterName='R';strcpy(ThreeLetterName,"ARG");

	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD";t.cs_name="CD";t.base_name="CB";atoms.push_back(t);
	t.name="NE";t.cs_name="NE";t.base_name="CB";atoms.push_back(t);
	t.name="CZ";t.cs_name="CZ";t.base_name="CB";atoms.push_back(t);
	t.name="NH1";t.cs_name="NH1";t.base_name="CB";atoms.push_back(t);
	t.name="NH2";t.cs_name="NH2";t.base_name="CB";atoms.push_back(t);
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=11;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=12;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="CG";
	t.proton_type=13;t.name="HG2";t.cs_name="HG2";t.base_name="HG";atoms.push_back(t);
	t.proton_type=14;t.name="HG3";t.cs_name="HG3";t.base_name="HG";atoms.push_back(t);
	t.carbon_name="CD";
	t.proton_type=15;t.name="HD2";t.cs_name="HD2";t.base_name="HD";atoms.push_back(t);
	t.proton_type=16;t.name="HD3";t.cs_name="HD3";t.base_name="HD";atoms.push_back(t);
	t.carbon_name="NE";
	t.proton_type=17;t.name="HE";t.cs_name="HE";t.base_name="HE";atoms.push_back(t);
	t.carbon_name="NH1";
	t.proton_type=18;t.name="HH11";t.cs_name="HH1";t.base_name="HH1";atoms.push_back(t);
	t.proton_type=18;t.name="HH12";t.cs_name="HH1";t.base_name="HH1";atoms.push_back(t);
	t.carbon_name="NH2";
	t.proton_type=19;t.name="HH21";t.cs_name="HH2";t.base_name="HH2";atoms.push_back(t);
	t.proton_type=19;t.name="HH22";t.cs_name="HH2";t.base_name="HH2";atoms.push_back(t);
	
	get_address("H")->cs_coil=8.24;
	get_address("N")->cs_coil=121.7;
	get_address("CA")->cs_coil=55.9;
	get_address("HA")->cs_coil=4.47;
	get_address("CB")->cs_coil=31.0;
	get_address("C")->cs_coil=175.3;

	get_address("H")->cs_mean=8.23;
	get_address("N")->cs_mean=120.76;
	get_address("CA")->cs_mean=56.82;
	get_address("HA")->cs_mean=4.29;
	get_address("CB")->cs_mean=30.65;
	get_address("C")->cs_mean=176.46;

	get_address("H")->cs_wang=8.21;
	get_address("N")->cs_wang=120.75;
	get_address("CA")->cs_wang=56.18;
	get_address("HA")->cs_wang=4.26;
	get_address("CB")->cs_wang=30.36;
	get_address("C")->cs_wang=176.01;

	coil_pre[0]=-0.45;
	coil_pre[1]=-0.07;
	coil_pre[2]= 0.00;
	coil_pre[3]= 0.01;
	coil_pre[4]=-0.03;
	coil_pre[5]= 0.01;

	coil_fol[0]=-0.09;
	coil_fol[1]= 0.16;
	coil_fol[2]= 0.19;
	coil_fol[3]=-0.12;
	coil_fol[4]=-0.01;
	coil_fol[5]=-0.02;
}
CArg::~CArg() {}


// ASN
void CAsn::schbond(vector<bbhbond_group> *grp)
{
	bbhbond_group t;
	t.id=residue;
	t.code=OneLetterName;
	t.npos=-1;
	t.hpos=-1;
	t.cpos=get("CG").index;
	t.opos=get("OD1").index;
	t.type=23; 

	grp->push_back(t);
	return;
}


struct noeatoms CAsn::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QD")==0)
	{
		tt.clear();
		tt.push_back(get("HD21").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HD22").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}



void CAsn::ani(vector<struct ani_group> *anistropy)
{
	struct ani_group t;
	CAminoacid::ani(anistropy);

	t.type=2;
	t.id=residue;
	t.code=OneLetterName;
	t.pos[0]=get("OD1").index;
	t.pos[1]=get("CG").index;
	t.pos[2]=get("ND2").index;
	anistropy->push_back(t);
	
	return ;
}



void CAsn::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("ND2").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}


CAsn::CAsn() 
{
	//CAminoacid::CAminoacid();
	OneLetterName='N';strcpy(ThreeLetterName,"ASN");
	
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="OD1";t.cs_name="OD1";t.base_name="CB";atoms.push_back(t);
	t.name="ND2";t.cs_name="ND2";t.base_name="CB";atoms.push_back(t);
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=20;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=21;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="ND2";
	t.proton_type=22;t.name="HD21";t.cs_name="HD21";t.base_name="HD2";atoms.push_back(t);
	t.proton_type=23;t.name="HD22";t.cs_name="HD22";t.base_name="HD2";atoms.push_back(t);

	get_address("H")->cs_coil=8.42;
	get_address("N")->cs_coil=119.4;
	get_address("CA")->cs_coil=53.0;
	get_address("HA")->cs_coil=4.75;
	get_address("CB")->cs_coil=38.9;
	get_address("C")->cs_coil=174.9;

	get_address("H")->cs_mean=8.33;
	get_address("N")->cs_mean=118.89;
	get_address("CA")->cs_mean=53.57;
	get_address("HA")->cs_mean=4.66;
	get_address("CB")->cs_mean=38.69;
	get_address("C")->cs_mean=175.30;


	get_address("H")->cs_wang=8.35;
	get_address("N")->cs_wang=118.50;
	get_address("CA")->cs_wang=53.00;
	get_address("HA")->cs_wang=4.66;
	get_address("CB")->cs_wang=38.43;
	get_address("C")->cs_wang=174.84;

	coil_pre[0]=-0.76;
	coil_pre[1]=-0.19;
	coil_pre[2]= 0.18;
	coil_pre[3]=-0.20;
	coil_pre[4]= 0.01;
	coil_pre[5]=-0.03;


	coil_fol[0]=-0.03;
	coil_fol[1]=-0.23;
	coil_fol[2]=+0.24;
	coil_fol[3]=-0.06;
	coil_fol[4]= 0.04;
	coil_fol[5]=-0.04;
}
	
CAsn::~CAsn() {}




//ASP
void CAsp::schbond(vector<bbhbond_group> *grp)
{
	bbhbond_group t;
	t.id=residue;
	t.code=OneLetterName;
	t.npos=-1;
	t.hpos=-1;
	t.cpos=get("CG").index;
	t.opos=get("OD1").index;
	t.type=22; 

	grp->push_back(t);
	return;
}

	

	


struct noeatoms CAsp::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}


void CAsp::ani(vector<struct ani_group> *anistropy)
{
	struct ani_group t;
	CAminoacid::ani(anistropy);


	t.type=3;
	t.id=residue;
	t.code=OneLetterName;
	t.pos[0]=get("OD1").index;
	t.pos[1]=get("CG").index;
	t.pos[2]=get("OD2").index;
	anistropy->push_back(t);
	
	return ;
}

void CAsp::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("OD1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}

CAsp::CAsp() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='D';strcpy(ThreeLetterName,"ASP");

	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;

	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="OD1";t.cs_name="OD1";t.base_name="CB";atoms.push_back(t);
	t.name="OD2";t.cs_name="OD2";t.base_name="CB";atoms.push_back(t);
	
	t.proton=1;t.carbon_name="CB";
	t.proton_type=24;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=25;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
		
	get_address("H")->cs_coil=8.31;
	get_address("N")->cs_coil=121.5;
	get_address("CA")->cs_coil=53.8;
	get_address("HA")->cs_coil=4.71;
	get_address("CB")->cs_coil=41.2;
	get_address("C")->cs_coil=175.7;

	get_address("H")->cs_mean=8.30;
	get_address("N")->cs_mean=120.63;
	get_address("CA")->cs_mean=54.71;
	get_address("HA")->cs_mean=4.59;
	get_address("CB")->cs_mean=40.87;
	get_address("C")->cs_mean=176.43;

	get_address("H")->cs_wang=8.31;
	get_address("N")->cs_wang=120.37;
	get_address("CA")->cs_wang=54.00;
	get_address("HA")->cs_wang=4.62;
	get_address("CB")->cs_wang=40.78;
	get_address("C")->cs_wang=176.00;

	coil_pre[0]=-0.43;
	coil_pre[1]= 0.07;
	coil_pre[2]= 0.20;
	coil_pre[3]=-0.07;
	coil_pre[4]=-0.01;
	coil_pre[5]=-0.03;

	coil_fol[0]= 0.23;
	coil_fol[1]=-0.11;
	coil_fol[2]= 0.28;
	coil_fol[3]= 0.11;
	coil_fol[4]= 0.04;
	coil_fol[5]=-0.03;
}
CAsp::~CAsp() {}


//CYS
struct noeatoms CCys::query(string name)
{	
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}




void CCys::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("SG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}

CCys::CCys() 
{
	//CAminoacid::CAminoacid();
	OneLetterName='C';strcpy(ThreeLetterName,"CYS");
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="SG";t.cs_name="SG";t.base_name="CB";atoms.push_back(t);
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=26;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=27;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="SG";
	t.proton_type=28;t.name="HG";t.cs_name="HG";t.base_name="HG";atoms.push_back(t);


	get_address("H")->cs_coil=8.20;
	get_address("N")->cs_coil=120.1;
	get_address("CA")->cs_coil=58.2;
	get_address("HA")->cs_coil=4.96;
	get_address("CB")->cs_coil=29.4;
	get_address("C")->cs_coil=174.7;

	get_address("H")->cs_mean=8.38;
	get_address("N")->cs_mean=120.09;
	get_address("CA")->cs_mean=58.31;
	get_address("HA")->cs_mean=4.64;
	get_address("CB")->cs_mean=32.66;
	get_address("C")->cs_mean=174.98;

	get_address("H")->cs_wang=8.10;
	get_address("N")->cs_wang=118.10;
	get_address("CA")->cs_wang=58.24;
	get_address("HA")->cs_wang=4.59;
	get_address("CB")->cs_wang=29.54;
	get_address("C")->cs_wang=175.11;


	coil_pre[0]= 1.36;
	coil_pre[1]=-0.00;
	coil_pre[2]= 0.44;
	coil_pre[3]=-0.18;
	coil_pre[4]= 0.03;
	coil_pre[5]=-0.00;

	coil_fol[0]=-1.17;
	coil_fol[1]= 0.10;
	coil_fol[2]= 0.17;
	coil_fol[3]= 0.21;
	coil_fol[4]= 0.01;
	coil_fol[5]= 0.03;

}
CCys::~CCys() {}




struct noeatoms CCyx::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}



void CCyx::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("SG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}


CCyx::CCyx() 
{
	//CAminoacid::CAminoacid();
	OneLetterName='C';strcpy(ThreeLetterName,"CYX");
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="SG";t.cs_name="SG";t.base_name="CB";atoms.push_back(t);
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=29;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=30;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="SG";
	t.proton_type=31;t.name="HG";t.cs_name="HG";t.base_name="HG";atoms.push_back(t);


	get_address("H")->cs_coil=8.20;
	get_address("N")->cs_coil=120.1;
	get_address("CA")->cs_coil=58.2;
	get_address("HA")->cs_coil=4.96;
	get_address("CB")->cs_coil=29.4;
	get_address("C")->cs_coil=174.7;

	get_address("H")->cs_mean=8.38;
	get_address("N")->cs_mean=120.09;
	get_address("CA")->cs_mean=58.31;
	get_address("HA")->cs_mean=4.64;
	get_address("CB")->cs_mean=32.66;
	get_address("C")->cs_mean=174.98;

	get_address("H")->cs_wang=8.10;
	get_address("N")->cs_wang=118.10;
	get_address("CA")->cs_wang=58.24;
	get_address("HA")->cs_wang=4.59;
	get_address("CB")->cs_wang=29.54;
	get_address("C")->cs_wang=175.11;

	coil_pre[0]= 1.36;
	coil_pre[1]=-0.00;
	coil_pre[2]= 0.44;
	coil_pre[3]=-0.18;
	coil_pre[4]= 0.03;
	coil_pre[5]=-0.00;

	coil_fol[0]=-1.17;
	coil_fol[1]= 0.10;
	coil_fol[2]= 0.17;
	coil_fol[3]= 0.21;
	coil_fol[4]= 0.01;
	coil_fol[5]= 0.03;
}
CCyx::~CCyx() {}


//GLN
void CGln::schbond(vector<bbhbond_group> *grp)
{
	bbhbond_group t;
	t.id=residue;
	t.code=OneLetterName;
	t.npos=-1;
	t.hpos=-1;
	t.cpos=get("CD").index;
	t.opos=get("OE1").index;
	t.type=23; 

	grp->push_back(t);
	return;
}


struct noeatoms CGln::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QG")==0)
	{
		tt.clear();
		tt.push_back(get("HG2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HG3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QE")==0)
	{
		tt.clear();
		tt.push_back(get("HE21").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HE22").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}

void CGln::ani(vector<struct ani_group> *anistropy)
{
	struct ani_group t;
	CAminoacid::ani(anistropy);

	t.type=2;
	t.id=residue;
	t.code=OneLetterName;
	t.pos[0]=get("OE1").index;
	t.pos[1]=get("CD").index;
	t.pos[2]=get("NE2").index;
	anistropy->push_back(t);
	
	return ;
}

void CGln::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("CD").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=5;
	t.x1=get("CB").index;
	t.x2=get("CG").index;
	t.x3=get("CD").index;
	t.x4=get("NE2").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}


CGln::CGln() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='Q';strcpy(ThreeLetterName,"GLN");
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD";t.cs_name="CD";t.base_name="CB";atoms.push_back(t);
	t.name="OE1";t.cs_name="OE1";t.base_name="CB";atoms.push_back(t);
	t.name="NE2";t.cs_name="NE2";t.base_name="CB";atoms.push_back(t);

	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=32;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=33;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="CG";
	t.proton_type=34;t.name="HG2";t.cs_name="HG2";t.base_name="HG";atoms.push_back(t);
	t.proton_type=35;t.name="HG3";t.cs_name="HG3";t.base_name="HG";atoms.push_back(t);
	t.carbon_name="NE2";
	t.proton_type=36;t.name="HE21";t.cs_name="HE21";t.base_name="HE2";atoms.push_back(t);
	t.proton_type=37;t.name="HE22";t.cs_name="HE22";t.base_name="HE2";atoms.push_back(t);

	get_address("H")->cs_coil=8.21;
	get_address("N")->cs_coil=120.2;
	get_address("CA")->cs_coil=55.5;
	get_address("HA")->cs_coil=4.43;
	get_address("CB")->cs_coil=29.4;
	get_address("C")->cs_coil=175.7;

	get_address("H")->cs_mean=8.21;
	get_address("N")->cs_mean=119.85;
	get_address("CA")->cs_mean=56.63;
	get_address("HA")->cs_mean=4.26;
	get_address("CB")->cs_mean=29.16;
	get_address("C")->cs_mean=176.37;

	get_address("H")->cs_wang=8.20;
	get_address("N")->cs_wang=119.82;
	get_address("CA")->cs_wang=55.89;
	get_address("HA")->cs_wang=4.29;
	get_address("CB")->cs_wang=29.01;
	get_address("C")->cs_wang=175.75;

	coil_pre[0]=-0.09;
	coil_pre[1]= 0.07;
	coil_pre[2]= 0.13;
	coil_pre[3]= 0.06;
	coil_pre[4]= 0.01;
	coil_pre[5]=-0.00;

	coil_fol[0]=-0.31;
	coil_fol[1]= 0.10;
	coil_fol[2]= 0.37;
	coil_fol[3]=-0.13;
	coil_fol[4]= 0.01;
	coil_fol[5]=-0.05;
}
CGln::~CGln() {}


//GLU
void CGlu::schbond(vector<bbhbond_group> *grp)
{
	bbhbond_group t;
	t.id=residue;
	t.code=OneLetterName;
	t.npos=-1;
	t.hpos=-1;
	t.cpos=get("CD").index;
	t.opos=get("OE1").index;
	t.type=22; 

	grp->push_back(t);
	return;
}


struct noeatoms CGlu::query(string name)
{
struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QG")==0)
	{
		tt.clear();
		tt.push_back(get("HG2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HG3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}

void CGlu::ani(vector<struct ani_group> *anistropy)
{
	struct ani_group t;
	CAminoacid::ani(anistropy);


	t.type=3;
	t.id=residue;
	t.code=OneLetterName;
	t.pos[0]=get("OE1").index;
	t.pos[1]=get("CD").index;
	t.pos[2]=get("OE2").index;
	anistropy->push_back(t);
	
	return ;
}

void CGlu::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("CD").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=5;
	t.x1=get("CB").index;
	t.x2=get("CG").index;
	t.x3=get("CD").index;
	t.x4=get("OE1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}


CGlu::CGlu() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='E';strcpy(ThreeLetterName,"GLU");
	struct Atom t;
	
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD";t.cs_name="CD";t.base_name="CB";atoms.push_back(t);
	t.name="OE1";t.cs_name="OE1";t.base_name="CB";atoms.push_back(t);
	t.name="OE2";t.cs_name="OE2";t.base_name="CB";atoms.push_back(t);

	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=38;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=39;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="CG";
	t.proton_type=40;t.name="HG2";t.cs_name="HG2";t.base_name="HG";atoms.push_back(t);
	t.proton_type=41;t.name="HG3";t.cs_name="HG3";t.base_name="HG";atoms.push_back(t);
	
	get_address("H")->cs_coil=8.36;
	get_address("N")->cs_coil=121.4;
	get_address("CA")->cs_coil=56.3;
	get_address("HA")->cs_coil=4.39;
	get_address("CB")->cs_coil=30.3;
	get_address("C")->cs_coil=175.9;
	
	get_address("H")->cs_mean=8.33;
	get_address("N")->cs_mean=120.66;
	get_address("CA")->cs_mean=57.37;
	get_address("HA")->cs_mean=4.25;
	get_address("CB")->cs_mean=29.98;
	get_address("C")->cs_mean=176.93;


	get_address("H")->cs_wang=8.36;
	get_address("N")->cs_wang=120.62;
	get_address("CA")->cs_wang=56.66;
	get_address("HA")->cs_wang=4.28;
	get_address("CB")->cs_wang=29.87;
	get_address("C")->cs_wang=176.32;

	coil_pre[0]=-0.36;
	coil_pre[1]= 0.04;
	coil_pre[2]= 0.01;
	coil_pre[3]= 0.04;
	coil_pre[4]=-0.01;
	coil_pre[5]=-0.01;

	coil_fol[0]= 0.26;
	coil_fol[1]= 0.14;
	coil_fol[2]= 0.25;
	coil_fol[3]= 0.06;
	coil_fol[4]= 0.07;
	coil_fol[5]=-0.04;
}
CGlu::~CGlu() {}

//GLY

struct noeatoms CGly::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QA")==0)
	{
		tt.clear();
		tt.push_back(get("HA2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HA3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}



void CGly::dihe(vector<dihe_group> * dihe_index)
{
	return;
}

CGly::CGly() 
{
    //CAminoacid::CAminoacid();
	bexploaded=0;
	atoms.clear();
	order_parameters.clear();
	OneLetterName='G';strcpy(ThreeLetterName,"GLY");
	
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=1;
	t.proton=0;
	t.ambig=0;
	t.proton=0;t.name="N";t.cs_name="N";t.base_name="N";atoms.push_back(t);
	t.proton=1;t.name="H";t.cs_name="H";t.base_name="H";atoms.push_back(t);
	t.proton=0;t.name="CA";t.cs_name="CA";t.base_name="CA";atoms.push_back(t);
	t.proton=0;t.name="C";t.cs_name="C";t.base_name="C";atoms.push_back(t);
	t.proton=0;t.name="O";t.cs_name="O";t.base_name="O";atoms.push_back(t);
	
	t.proton=1;
	t.carbon_name="CA";
	t.proton_type=42;t.name="HA2";t.cs_name="HA2";t.base_name="HA";atoms.push_back(t);
	t.proton_type=43;t.name="HA3";t.cs_name="HA3";t.base_name="HA";atoms.push_back(t);


	struct S2 s2;
	s2.name1="N";s2.name2="H";s2.exp=999.0;s2.pre=999.0;order_parameters.push_back(s2);
	s2.name1="N";s2.name2="CA";s2.exp=999.0;s2.pre=999.0;order_parameters.push_back(s2);
	s2.name1="CA";s2.name2="C";s2.exp=999.0;s2.pre=999.0;order_parameters.push_back(s2);
	s2.name1="C0";s2.name2="N";s2.exp=999.0;s2.pre=999.0;order_parameters.push_back(s2);

	get_address("H")->cs_coil=8.31;
	get_address("N")->cs_coil=109.8;
	get_address("CA")->cs_coil=45.2;
	get_address("HA1")->cs_coil=4.12;
	get_address("HA2")->cs_coil=4.12;
	get_address("C")->cs_coil=173.9;

	get_address("H")->cs_mean=8.33;
	get_address("N")->cs_mean=109.60;
	get_address("CA")->cs_mean=45.36;
	get_address("HA1")->cs_mean=3.97;
	get_address("HA2")->cs_mean=3.90;
	get_address("C")->cs_mean=173.38;

	get_address("H")->cs_wang=8.37;
	get_address("N")->cs_wang=109.48;
	get_address("CA")->cs_wang=45.28;
	get_address("HA1")->cs_wang=3.97;
	get_address("HA2")->cs_wang=3.97;
	get_address("C")->cs_wang=174.01;

	coil_pre[0]=-0.43;
	coil_pre[1]=-0.09;
	coil_pre[2]=-0.26;
	coil_pre[3]= 0.17;
	coil_pre[4]=-0.10;
	coil_pre[5]= 0.04;


	coil_fol[0]= 0.13;
	coil_fol[1]= 0.47;
	coil_fol[2]= 0.12;
	coil_fol[3]=-0.07;
	coil_fol[4]= 0.06;
	coil_fol[5]=-0.03;
}
CGly::~CGly() {}



//HIS
struct noeatoms CHis::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}

void CHis::ring(vector<ring_group> *ring)
{
	ring_group t;
	t.x1=3; //type
	t.x2=get("CG").index;
	t.x3=get("ND1").index;
	t.x4=get("CE1").index;
	t.x5=get("NE2").index;
	t.x6=get("CD2").index;
	ring->push_back(t);
	return ;
}

void CHis::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("ND1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}

CHis::CHis() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='H';strcpy(ThreeLetterName,"HIS");

	struct Atom t;
	
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD2";t.cs_name="CD2";t.base_name="CB";atoms.push_back(t);
	t.name="ND1";t.cs_name="ND1";t.base_name="CB";atoms.push_back(t);
	t.name="CE1";t.cs_name="CE1";t.base_name="CB";atoms.push_back(t);
	t.name="NE2";t.cs_name="NE2";t.base_name="CB";atoms.push_back(t);
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=44;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=45;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="ND1";
	t.proton_type=46;t.name="HD1";t.cs_name="HD1";t.base_name="HD1";atoms.push_back(t);
	t.carbon_name="CD2";
	t.proton_type=47;t.name="HD2";t.cs_name="HD2";t.base_name="HD2";atoms.push_back(t);
	t.carbon_name="CE1";
	t.proton_type=48;t.name="HE1";t.cs_name="HE1";t.base_name="HE1";atoms.push_back(t);
	t.carbon_name="NE2";
	t.proton_type=49;t.name="HE2";t.cs_name="HE2";t.base_name="HE2";atoms.push_back(t);
		
	get_address("H")->cs_coil=8.29;
	get_address("N")->cs_coil=119.7;
	get_address("CA")->cs_coil=55.3;
	get_address("HA")->cs_coil=4.73;
	get_address("CB")->cs_coil=30.1;
	get_address("C")->cs_coil=174.4;
	
	get_address("H")->cs_mean=8.24;
	get_address("N")->cs_mean=119.67;
	get_address("CA")->cs_mean=56.55;
	get_address("HA")->cs_mean=4.60;
	get_address("CB")->cs_mean=30.22;
	get_address("C")->cs_mean=175.28;

	get_address("H")->cs_wang=8.18;
	get_address("N")->cs_wang=118.92;
	get_address("CA")->cs_wang=55.74;
	get_address("HA")->cs_wang=4.60;
	get_address("CB")->cs_wang=29.50;
	get_address("C")->cs_wang=174.78;

	coil_pre[0]=-0.05;
	coil_pre[1]=-0.13;
	coil_pre[2]= 0.16;
	coil_pre[3]= 0.07;
	coil_pre[4]=-0.01;
	coil_pre[5]=-0.05;

	coil_fol[0]=-0.09;
	coil_fol[1]=-0.05;
	coil_fol[2]= 0.22;
	coil_fol[3]=-0.24;
	coil_fol[4]=-0.01;
	coil_fol[5]=-0.09;

}
CHis::~CHis() {}



//ILE

struct noeatoms CIle::query(string name)
{
struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QG")==0)
	{
		tt.clear();
		tt.push_back(get("HG12").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HG13").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("MG")==0 || name.compare("QG2")==0)
	{
		tt.push_back(get("HG21").index);
		tt.push_back(get("HG22").index);
		tt.push_back(get("HG23").index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	else if(name.compare("MD")==0|| name.compare("QD1")==0)
	{
		tt.push_back(get("HD11").index);
		tt.push_back(get("HD12").index);
		tt.push_back(get("HD13").index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}

void CIle::proton2(vector<struct proton> *sel)
{
	struct proton t;

	t.type=9;
	t.name="HG2";
	get_proton(&t);
	sel->push_back(t);

	t.type=10;
	t.name="HD1";
	get_proton(&t);
	sel->push_back(t);

	return;
}


void CIle::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG1").index;
	t.x4=get("CD1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);
	
	return;
}



CIle::CIle() 
{
	//CAminoacid::CAminoacid();
	OneLetterName='I';strcpy(ThreeLetterName,"ILE");
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG2";t.cs_name="CG2";t.base_name="CB";atoms.push_back(t);
	t.name="CG1";t.cs_name="CG1";t.base_name="CB";atoms.push_back(t);
	t.name="CD1";t.cs_name="CD1";t.base_name="CB";atoms.push_back(t);
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=50;t.name="HB";t.cs_name="HB";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="CG2";
	t.proton_type=9;t.name="HG21";t.cs_name="HG2";t.base_name="HG2";atoms.push_back(t);
	t.proton_type=9;t.name="HG22";t.cs_name="HG2";t.base_name="HG2";atoms.push_back(t);
	t.proton_type=9;t.name="HG23";t.cs_name="HG2";t.base_name="HG2";atoms.push_back(t);
	
	t.carbon_name="CG1";
	t.proton_type=51;t.name="HG12";t.cs_name="HG12";t.base_name="HG1";atoms.push_back(t);
	t.proton_type=52;t.name="HG13";t.cs_name="HG13";t.base_name="HG1";atoms.push_back(t);

	t.carbon_name="CD1";
	t.proton_type=10;t.name="HD11";t.cs_name="HD1";t.base_name="HD1";atoms.push_back(t);
	t.proton_type=10;t.name="HD12";t.cs_name="HD1";t.base_name="HD1";atoms.push_back(t);
	t.proton_type=10;t.name="HD13";t.cs_name="HD1";t.base_name="HD1";atoms.push_back(t);
		
	get_address("H")->cs_coil=8.30;
	get_address("N")->cs_coil=122.0;
	get_address("CA")->cs_coil=60.4;
	get_address("HA")->cs_coil=4.31;
	get_address("CB")->cs_coil=38.7;
	get_address("C")->cs_coil=174.9;
		
	get_address("H")->cs_mean=8.22;
	get_address("N")->cs_mean=121.79;
	get_address("CA")->cs_mean=55.70;
	get_address("HA")->cs_mean=4.30;
	get_address("CB")->cs_mean=42.27;
	get_address("C")->cs_mean=177.06;

	get_address("H")->cs_wang=7.94;
	get_address("N")->cs_wang=120.58;
	get_address("CA")->cs_wang=60.79;
	get_address("HA")->cs_wang=4.18;
	get_address("CB")->cs_wang=38.43;
	get_address("C")->cs_wang=175.52;


	coil_pre[0]= 2.92;
	coil_pre[1]= 0.11;
	coil_pre[2]=-0.15;
	coil_pre[3]= 0.26;
	coil_pre[4]= 0.12;
	coil_pre[5]= 0.04;

	coil_fol[0]=-0.20;
	coil_fol[1]=-0.09;
	coil_fol[2]= 0.03;
	coil_fol[3]= 0.28;
	coil_fol[4]=-0.01;
	coil_fol[5]= 0.03;
}
CIle::~CIle() {}




//LEU

void CLeu::methyl_ambig(int flag)
{
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).base_name=="HD1" && atoms.at(i).ambig==2)
		{
			if(atoms.at(i).cs_exp>0.0)
				atoms.at(i).cs_exp=999.0;
		}
		if(atoms.at(i).base_name=="HD2" && atoms.at(i).ambig==2)
		{
			if(atoms.at(i).cs_exp>0.0)
				atoms.at(i).cs_exp=999.0;
		}
	}
	return;
}

struct noeatoms CLeu::query(string name)
{
struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("MD2")==0 || name.compare("QD2")==0)
	{
		tt.push_back(get("HD21").index);
		tt.push_back(get("HD22").index);
		tt.push_back(get("HD23").index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	else if(name.compare("MD1")==0|| name.compare("QD1")==0)
	{
		tt.push_back(get("HD11").index);
		tt.push_back(get("HD12").index);
		tt.push_back(get("HD13").index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	else if(name.compare("QD")==0)
	{
		tt.clear();
		tt.push_back(get("HD21").index);
		tt.push_back(get("HD22").index);
		tt.push_back(get("HD23").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HD11").index);
		tt.push_back(get("HD12").index);
		tt.push_back(get("HD13").index);
		t.atoms.push_back(tt);
		t.length=5.0;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}

void CLeu::proton2(vector<struct proton> *sel)
{
	struct proton t;


	t.type=7;
	t.name="HD1";
	get_proton(&t);
	sel->push_back(t);

	t.type=8;
	t.name="HD2";
	get_proton(&t);
	sel->push_back(t);

	return;
}


void CLeu::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;	
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("CD1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}

CLeu::CLeu() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='L';strcpy(ThreeLetterName,"LEU");
	
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;
	t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD1";t.cs_name="CD1";t.base_name="CB";atoms.push_back(t);
	t.name="CD2";t.cs_name="CD2";t.base_name="CB";atoms.push_back(t);
	
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=53;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=54;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	
	t.carbon_name="CG";
	t.proton_type=55;t.name="HG";t.cs_name="HG";t.base_name="HG";atoms.push_back(t);
	
	t.carbon_name="CD1";
	t.proton_type=7;t.name="HD11";t.cs_name="HD1";t.base_name="HD1";atoms.push_back(t);
	t.proton_type=7;t.name="HD12";t.cs_name="HD1";t.base_name="HD1";atoms.push_back(t);
	t.proton_type=7;t.name="HD13";t.cs_name="HD1";t.base_name="HD1";atoms.push_back(t);

	t.carbon_name="CD2";
	t.proton_type=8;t.name="HD21";t.cs_name="HD2";t.base_name="HD2";atoms.push_back(t);
	t.proton_type=8;t.name="HD22";t.cs_name="HD2";t.base_name="HD2";atoms.push_back(t);
	t.proton_type=8;t.name="HD23";t.cs_name="HD2";t.base_name="HD2";atoms.push_back(t);
		
	get_address("H")->cs_coil=8.18;
	get_address("N")->cs_coil=122.3;
	get_address("CA")->cs_coil=54.5;
	get_address("HA")->cs_coil=4.47;
	get_address("CB")->cs_coil=42.5;
	get_address("C")->cs_coil=176.4;

		
	get_address("H")->cs_mean=8.22;
	get_address("N")->cs_mean=121.79;
	get_address("CA")->cs_mean=55.70;
	get_address("HA")->cs_mean=4.30;
	get_address("CB")->cs_mean=42.27;
	get_address("C")->cs_mean=177.06;

	get_address("H")->cs_wang=8.06;
	get_address("N")->cs_wang=121.57;
	get_address("CA")->cs_wang=54.77;
	get_address("HA")->cs_wang=4.36;
	get_address("CB")->cs_wang=42.14;
	get_address("C")->cs_wang=176.70;


	coil_pre[0]=-0.76;
	coil_pre[1]= 0.13;
	coil_pre[2]=-0.07;
	coil_pre[3]=-0.03;
	coil_pre[4]=-0.06;
	coil_pre[5]= 0.00;

	coil_fol[0]=-0.49;
	coil_fol[1]= 0.06;
	coil_fol[2]= 0.10;
	coil_fol[3]=-0.10;
	coil_fol[4]=-0.07;
	coil_fol[5]=-0.02;
}
CLeu::~CLeu() {}




//LYS

struct noeatoms CLys::query(string name)
{
struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QG")==0)
	{
		tt.clear();
		tt.push_back(get("HG2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HG3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QD")==0)
	{
		tt.clear();
		tt.push_back(get("HD2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HD3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QE")==0)
	{
		tt.clear();
		tt.push_back(get("HE2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HE3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}

void CLys::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("CD").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=5;
	t.x1=get("CB").index;
	t.x2=get("CG").index;
	t.x3=get("CD").index;
	t.x4=get("CE").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=6;
	t.x1=get("CG").index;
	t.x2=get("CD").index;
	t.x3=get("CE").index;
	t.x4=get("NZ").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}


CLys::CLys() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='K';strcpy(ThreeLetterName,"LYS");
	
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD";t.cs_name="CD";t.base_name="CB";atoms.push_back(t);
	t.name="CE";t.cs_name="CE";t.base_name="CB";atoms.push_back(t);
	t.name="NZ";t.cs_name="NZ";t.base_name="CB";atoms.push_back(t);

	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=56;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=57;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="CG";
	t.proton_type=58;t.name="HG2";t.cs_name="HG2";t.base_name="HG";atoms.push_back(t);
	t.proton_type=59;t.name="HG3";t.cs_name="HG3";t.base_name="HG";atoms.push_back(t);
	t.carbon_name="CD";
	t.proton_type=60;t.name="HD2";t.cs_name="HD2";t.base_name="HD";atoms.push_back(t);
	t.proton_type=61;t.name="HD3";t.cs_name="HD3";t.base_name="HD";atoms.push_back(t);
	t.carbon_name="CE";
	t.proton_type=62;t.name="HE2";t.cs_name="HE2";t.base_name="HE";atoms.push_back(t);
	t.proton_type=63;t.name="HE3";t.cs_name="HE3";t.base_name="HE";atoms.push_back(t);
	t.carbon_name="NZ";
	t.proton_type=64;t.name="HZ1";t.cs_name="HZ";t.base_name="HZ";atoms.push_back(t);
	t.proton_type=64;t.name="HZ2";t.cs_name="HZ";t.base_name="HZ";atoms.push_back(t);
	t.proton_type=64;t.name="HZ3";t.cs_name="HZ";t.base_name="HZ";atoms.push_back(t);
		
	get_address("H")->cs_coil=8.24;
	get_address("N")->cs_coil=121.8;
	get_address("CA")->cs_coil=56.2;
	get_address("HA")->cs_coil=4.36;
	get_address("CB")->cs_coil=32.8;
	get_address("C")->cs_coil=176.0;
		
	get_address("H")->cs_mean=8.18;
	get_address("N")->cs_mean=121.00;
	get_address("CA")->cs_mean=57.00;
	get_address("HA")->cs_mean=4.26;
	get_address("CB")->cs_mean=32.78;
	get_address("C")->cs_mean=176.71;

	get_address("H")->cs_wang=8.17;
	get_address("N")->cs_wang=121.1;
	get_address("CA")->cs_wang=56.29;
	get_address("HA")->cs_wang=4.28;
	get_address("CB")->cs_wang=32.53;
	get_address("C")->cs_wang=176.15;

	coil_pre[0]=-0.26;
	coil_pre[1]=-0.09;
	coil_pre[2]=-0.07;
	coil_pre[3]= 0.08;
	coil_pre[4]=-0.02;
	coil_pre[5]=-0.00;

	coil_fol[0]=-0.13;
	coil_fol[1]=-0.13;
	coil_fol[2]= 0.08;
	coil_fol[3]= 0.01;
	coil_fol[4]=-0.04;
	coil_fol[5]=-0.03;
}
CLys::~CLys() {}



//MET

struct noeatoms CMet::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QG")==0)
	{
		tt.clear();
		tt.push_back(get("HG2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HG3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("ME")==0 || name.compare("QE")==0)
	{
		tt.clear();
		tt.push_back(get("HE1").index);
		tt.push_back(get("HE2").index);
		tt.push_back(get("HE3").index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}


void CMet::proton2(vector<struct proton> *sel)
{
	struct proton t;


	t.type=3;
	t.name="HE";
	get_proton(&t);
	sel->push_back(t);

	return;
}


void CMet::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("SD").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=5;
	t.x1=get("CB").index;
	t.x2=get("CG").index;
	t.x3=get("SD").index;
	t.x4=get("CE").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}

CMet::CMet() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='M';strcpy(ThreeLetterName,"MET");

	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="SD";t.cs_name="SD";t.base_name="CB";atoms.push_back(t);
	t.name="CE";t.cs_name="CE";t.base_name="CB";atoms.push_back(t);

	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=65;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=66;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="CG";
	t.proton_type=67;t.name="HG2";t.cs_name="HG2";t.base_name="HG";atoms.push_back(t);
	t.proton_type=68;t.name="HG3";t.cs_name="HG3";t.base_name="HG";atoms.push_back(t);
	
	t.carbon_name="CE";
	t.proton_type=3;t.name="HE1";t.cs_name="HE";t.base_name="HE";atoms.push_back(t);
	t.proton_type=3;t.name="HE2";t.cs_name="HE";t.base_name="HE";atoms.push_back(t);
	t.proton_type=3;t.name="HE3";t.cs_name="HE";t.base_name="HE";atoms.push_back(t);
	
	get_address("H")->cs_coil=8.35;
	get_address("N")->cs_coil=121.2;
	get_address("CA")->cs_coil=55.4;
	get_address("HA")->cs_coil=4.41;
	get_address("CB")->cs_coil=33.7;
	get_address("C")->cs_coil=174.6;

		
	get_address("H")->cs_mean=8.25;
	get_address("N")->cs_mean=120.07;
	get_address("CA")->cs_mean=56.15;
	get_address("HA")->cs_mean=4.40;
	get_address("CB")->cs_mean=32.95;
	get_address("C")->cs_mean=176.24;

	get_address("H")->cs_wang=8.22;
	get_address("N")->cs_wang=120.14;
	get_address("CA")->cs_wang=55.43;
	get_address("HA")->cs_wang=4.47;
	get_address("CB")->cs_wang=32.92;
	get_address("C")->cs_wang=175.94;

	coil_pre[0]= 0.69;
	coil_pre[1]= 0.10;
	coil_pre[2]= 0.09;
	coil_pre[3]= 0.10;
	coil_pre[4]= 0.05;
	coil_pre[5]= 0.05;

	coil_fol[0]=-0.02;
	coil_fol[1]= 0.19;
	coil_fol[2]= 0.22;
	coil_fol[3]= 0.06;
	coil_fol[4]= 0.01;
	coil_fol[5]=-0.01;
}
CMet::~CMet() {}


//PHE 
//special handling for aromatic ring's HD and HE atoms.
/*void CPhe::proton(vector<struct proton> *sel)
{
	struct proton t;
	int i;
	string oldname;

	oldname="XXX";
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).bb==0 && atoms.at(i).proton==1 && atoms.at(i).cs_name!=oldname)
		{
			oldname=atoms.at(i).cs_name;
			t.name=atoms.at(i).cs_name;
			if(get_proton(&t))
				sel->push_back(t);
		}
	}
	return;
}

void CPhe::proton3(vector<struct proton> *sel)
{
	proton(sel);
	return;
}*/


struct noeatoms CPhe::query(string name)
{
struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QD")==0)
	{
		tt.clear();
		tt.push_back(get("HD1").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HD2").index);
		t.atoms.push_back(tt);
		t.length=5.0;
	}
	else if(name.compare("QE")==0)
	{
		tt.clear();
		tt.push_back(get("HE1").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HE2").index);
		t.atoms.push_back(tt);
		t.length=5.0;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}

void CPhe::ring(vector<ring_group> *ring)
{
	ring_group t;

	t.x1=1; //type
	t.x2=get("CG").index;
	t.x3=get("CD1").index;
	t.x4=get("CE1").index;
	t.x5=get("CZ").index;
	t.x6=get("CE2").index;
	t.x7=get("CD2").index;

	ring->push_back(t);
	return ;
}

void CPhe::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("CD1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}


CPhe::CPhe() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='F';strcpy(ThreeLetterName,"PHE");
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD1";t.cs_name="CD";t.base_name="CB";atoms.push_back(t);
	t.name="CD2";t.cs_name="CD";t.base_name="CB";atoms.push_back(t);
	t.name="CE1";t.cs_name="CE";t.base_name="CB";atoms.push_back(t);
	t.name="CE2";t.cs_name="CE";t.base_name="CB";atoms.push_back(t);
	t.name="CZ";t.cs_name="CZ";t.base_name="CB";atoms.push_back(t);
	
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=69;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=70;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);

	t.carbon_name="CD1";
	t.proton_type=71;t.name="HD1";t.cs_name="HD";t.base_name="HD1";atoms.push_back(t);
	t.carbon_name="CD2";
	t.proton_type=71;t.name="HD2";t.cs_name="HD";t.base_name="HD2";atoms.push_back(t);
	t.carbon_name="CE1";
	t.proton_type=72;t.name="HE1";t.cs_name="HE";t.base_name="HE1";atoms.push_back(t);
	t.carbon_name="CE2";
	t.proton_type=72;t.name="HE2";t.cs_name="HE";t.base_name="HE2";atoms.push_back(t);
	t.carbon_name="CZ";
	t.proton_type=73;t.name="HZ";t.cs_name="HZ";t.base_name="HZ";atoms.push_back(t);

	get_address("H")->cs_coil=8.27;
	get_address("N")->cs_coil=120.1;
	get_address("CA")->cs_coil=57.2;
	get_address("HA")->cs_coil=4.65;
	get_address("CB")->cs_coil=40.2;
	get_address("C")->cs_coil=175.0;
		
	get_address("H")->cs_mean=8.34;
	get_address("N")->cs_mean=120.39;
	get_address("CA")->cs_mean=58.16;
	get_address("HA")->cs_mean=4.62;
	get_address("CB")->cs_mean=39.95;
	get_address("C")->cs_mean=175.48;

	get_address("H")->cs_wang=8.09;
	get_address("N")->cs_wang=119.72;
	get_address("CA")->cs_wang=57.46;
	get_address("HA")->cs_wang=4.59;
	get_address("CB")->cs_wang=39.41;
	get_address("C")->cs_wang=175.46;

	coil_pre[0]= 0.20;
	coil_pre[1]=-0.14;
	coil_pre[2]= 0.02;
	coil_pre[3]=-0.06;
	coil_pre[4]= 0.01;
	coil_pre[5]=-0.01;

	coil_fol[0]=-0.35;
	coil_fol[1]=-0.22;
	coil_fol[2]=-0.04;
	coil_fol[3]=-0.09;
	coil_fol[4]=-0.04;
	coil_fol[5]=-0.02;
}
CPhe::~CPhe() {}


//PRO

struct noeatoms CPro::query(string name)
{
struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QG")==0)
	{
		tt.clear();
		tt.push_back(get("HG2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HG3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QD")==0)
	{
		tt.clear();
		tt.push_back(get("HD2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HD3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}



void CPro::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	//dihe_index->push_back(t);

	t.type=3;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("CD").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=5;
	t.x1=get("CB").index;
	t.x2=get("CG").index;
	t.x3=get("CD").index;
	t.x4=get("n").index;
	//dihe_index->push_back(t);

	return;
}

CPro::CPro() 
{
    //CAminoacid::CAminoacid();
	bexploaded=0;
	atoms.clear();order_parameters.clear();
	OneLetterName='P';strcpy(ThreeLetterName,"PRO");

	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=1;
	t.proton=0;t.ambig=0;
	t.proton=0;t.name="N";t.cs_name="N";t.base_name="N";atoms.push_back(t);
	t.proton=0;t.name="CA";t.cs_name="CA";t.base_name="CA";atoms.push_back(t);
	t.proton=1;t.name="HA";t.cs_name="HA";t.base_name="HA";atoms.push_back(t);
	t.proton=0;t.name="C";t.cs_name="C";t.base_name="C";atoms.push_back(t);
	t.proton=0;t.name="O";t.cs_name="O";t.base_name="O";atoms.push_back(t);

	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;
	t.bb=0;
	t.proton=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD";t.cs_name="CD";t.base_name="CB";atoms.push_back(t);
	
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=74;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=75;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="CG";
	t.proton_type=76;t.name="HG2";t.cs_name="HG2";t.base_name="HG";atoms.push_back(t);
	t.proton_type=77;t.name="HG3";t.cs_name="HG3";t.base_name="HG";atoms.push_back(t);
	t.carbon_name="CD";
	t.proton_type=78;t.name="HD2";t.cs_name="HD2";t.base_name="HD";atoms.push_back(t);
	t.proton_type=79;t.name="HD3";t.cs_name="HD3";t.base_name="HD";atoms.push_back(t);


	get_address("CA")->cs_coil=62.6;
	get_address("HA")->cs_coil=4.44;
	get_address("CB")->cs_coil=31.9;
	get_address("C")->cs_coil=176.1;

		
	get_address("CA")->cs_mean=63.36;
	get_address("HA")->cs_mean=4.39;
	get_address("CB")->cs_mean=31.85;
	get_address("C")->cs_mean=176.76;

	get_address("CA")->cs_wang=63.24;
	get_address("HA")->cs_wang=4.41;
	get_address("CB")->cs_wang=31.81;
	get_address("C")->cs_wang=176.62;

	coil_pre[0]=-0.94;
	coil_pre[1]= 0.21;
	coil_pre[2]= 0.07;
	coil_pre[3]=-0.16;
	coil_pre[4]= 0.14;
	coil_pre[5]=-0.04;

	coil_fol[0]= 0.92;
	coil_fol[1]=-1.19;
	coil_fol[2]=-2.04;
	coil_fol[3]=-0.20;
	coil_fol[4]=-0.17;
	coil_fol[5]= 0.21;

}
CPro::~CPro() {}


//SER
void CSer::schbond(vector<bbhbond_group> *grp)
{
	bbhbond_group t;
	t.id=residue;
	t.code=OneLetterName;
	t.cpos=-1;
	t.opos=-1;
	t.npos=get("OG").index;
	t.hpos=get("HG").index;
	t.type=12; 

	grp->push_back(t);
	return;
}

struct noeatoms CSer::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}


void CSer::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("OG").index;	
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}



CSer::CSer() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='S';strcpy(ThreeLetterName,"SER");
	
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="OG";t.cs_name="OG";t.base_name="CB";atoms.push_back(t);
	
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=80;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=81;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="OG";
	t.proton_type=82;t.name="HG";t.cs_name="HG";t.base_name="HG";atoms.push_back(t);

	get_address("H")->cs_coil=8.36;
	get_address("N")->cs_coil=116.8;
	get_address("CA")->cs_coil=58.1;
	get_address("HA")->cs_coil=4.55;
	get_address("CB")->cs_coil=64.1;
	get_address("C")->cs_coil=174.2;

		
	get_address("H")->cs_mean=8.28;
	get_address("N")->cs_mean=116.27;
	get_address("CA")->cs_mean=58.76;
	get_address("HA")->cs_mean=4.47;
	get_address("CB")->cs_mean=63.79;
	get_address("C")->cs_mean=174.65;

	get_address("H")->cs_wang=8.22;
	get_address("N")->cs_wang=116.00;
	get_address("CA")->cs_wang=58.20;
	get_address("HA")->cs_wang=4.45;
	get_address("CB")->cs_wang=63.75;
	get_address("C")->cs_wang=174.41;

	coil_pre[0]= 1.16;
	coil_pre[1]=-0.10;
	coil_pre[2]= 0.11;
	coil_pre[3]=-0.06;
	coil_pre[4]= 0.02;
	coil_pre[5]= 0.01;

	coil_fol[0]= 0.30;
	coil_fol[1]= 0.10;
	coil_fol[2]= 0.10;
	coil_fol[3]= 0.14;
	coil_fol[4]= 0.07;
	coil_fol[5]= 0.02;
}
CSer::~CSer() {}


//THR
void CThr::schbond(vector<bbhbond_group> *grp)
{
	bbhbond_group t;
	t.id=residue;
	t.code=OneLetterName;
	t.cpos=-1;
	t.opos=-1;
	t.npos=get("OG1").index;
	t.hpos=get("HG1").index;
	t.type=12; 

	grp->push_back(t);
	return;
}



struct noeatoms CThr::query(string name)
{
	struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("MG")==0|| name.compare("QG2")==0)
	{
		tt.push_back(get("HG21").index);
		tt.push_back(get("HG22").index);
		tt.push_back(get("HG23").index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;;
}

void CThr::proton2(vector<struct proton> *sel)
{
	struct proton t;

	t.type=2;
	t.name="HG2";
	get_proton(&t);
	sel->push_back(t);


	return;
}


void CThr::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG2").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}

CThr::CThr() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='T';strcpy(ThreeLetterName,"THR");
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="OG1";t.cs_name="OG1";t.base_name="CB";atoms.push_back(t);
	t.name="CG2";t.cs_name="CG2";t.base_name="CB";atoms.push_back(t);
	
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=83;t.name="HB";t.cs_name="HB";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="OG1";
	t.proton_type=84;t.name="HG1";t.cs_name="HG1";t.base_name="HG1";atoms.push_back(t);

	t.carbon_name="CG2";
	t.proton_type=2;t.name="HG21";t.cs_name="HG2";t.base_name="HG2";atoms.push_back(t);
	t.proton_type=2;t.name="HG22";t.cs_name="HG2";t.base_name="HG2";atoms.push_back(t);
	t.proton_type=2;t.name="HG23";t.cs_name="HG2";t.base_name="HG2";atoms.push_back(t);
		
	get_address("H")->cs_coil=8.27;
	get_address("N")->cs_coil=114.6;
	get_address("CA")->cs_coil=60.9;
	get_address("HA")->cs_coil=4.55;
	get_address("CB")->cs_coil=69.7;
	get_address("C")->cs_coil=174.5;
		
	get_address("H")->cs_mean=8.24;
	get_address("N")->cs_mean=115.36;
	get_address("CA")->cs_mean=62.28;
	get_address("HA")->cs_mean=4.45;
	get_address("CB")->cs_mean=69.72;
	get_address("C")->cs_mean=174.57;


	get_address("H")->cs_wang=8.16;
	get_address("N")->cs_wang=113.88;
	get_address("CA")->cs_wang=61.30;
	get_address("HA")->cs_wang=4.44;
	get_address("CB")->cs_wang=68.92;
	get_address("C")->cs_wang=174.78;

	coil_pre[0]= 1.23;
	coil_pre[1]=-0.07;
	coil_pre[2]= 0.05;
	coil_pre[3]=-0.11;
	coil_pre[4]= 0.02;
	coil_pre[5]=-0.00;


	coil_fol[0]= 0.22;
	coil_fol[1]= 0.10;
	coil_fol[2]= 0.02;
	coil_fol[3]= 0.14;
	coil_fol[4]= 0.04;
	coil_fol[5]= 0.08;
	
}
CThr::~CThr() {}


//TRP
void CTrp::schbond(vector<bbhbond_group> *grp)
{
	bbhbond_group t;
	t.id=residue;
	t.code=OneLetterName;
	t.cpos=-1;
	t.opos=-1;
	t.npos=get("NE1").index;
	t.hpos=get("HE1").index;
	t.type=13; 

	grp->push_back(t);
	return;
}

struct noeatoms CTrp::query(string name)
{
struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}


void CTrp::ring(vector<ring_group> *ring)
{
	ring_group t;

	t.x1=4; //type
	t.x2=get("CG").index;
	t.x3=get("CD1").index;
	t.x4=get("NE1").index;
	t.x5=get("CE2").index;
	t.x6=get("CD2").index;
	ring->push_back(t);

	t.x1=5;
	t.x2=get("CD2").index;
	t.x3=get("CE2").index;
	t.x4=get("CZ2").index;
	t.x5=get("CH2").index;
	t.x6=get("CZ3").index;
	t.x7=get("CE3").index;

	ring->push_back(t);
	return ;
}


void CTrp::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("CD1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);
	
	return;
}



CTrp::CTrp() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='W';strcpy(ThreeLetterName,"TRP");
	
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD1";t.cs_name="CD1";t.base_name="CB";atoms.push_back(t);
	t.name="CD2";t.cs_name="CD2";t.base_name="CB";atoms.push_back(t);
	t.name="NE1";t.cs_name="NE1";t.base_name="CB";atoms.push_back(t);
	t.name="CE2";t.cs_name="CE2";t.base_name="CB";atoms.push_back(t);
	t.name="CE3";t.cs_name="CE3";t.base_name="CB";atoms.push_back(t);	
	t.name="CZ2";t.cs_name="CZ2";t.base_name="CB";atoms.push_back(t);
	t.name="CZ3";t.cs_name="CZ3";t.base_name="CB";atoms.push_back(t);
	t.name="CH2";t.cs_name="CH2";t.base_name="CB";atoms.push_back(t);
	
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=85;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=86;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);

	t.carbon_name="CD1";
	t.proton_type=87;t.name="HD1";t.cs_name="HD1";t.base_name="HD1";atoms.push_back(t);
	t.carbon_name="NE1";
	t.proton_type=88;t.name="HE1";t.cs_name="HE1";t.base_name="HE1";atoms.push_back(t);
	t.carbon_name="CE3";
	t.proton_type=89;t.name="HE3";t.cs_name="HE3";t.base_name="HE3";atoms.push_back(t);
	t.carbon_name="CZ2";
	t.proton_type=90;t.name="HZ2";t.cs_name="HZ2";t.base_name="HZ2";atoms.push_back(t);
	t.carbon_name="CZ3";
	t.proton_type=91;t.name="HZ3";t.cs_name="HZ3";t.base_name="HZ3";atoms.push_back(t);
	t.carbon_name="CH2";
	t.proton_type=92;t.name="HH2";t.cs_name="HH2";t.base_name="HH2";atoms.push_back(t);

	get_address("H")->cs_coil=8.19;
	get_address("N")->cs_coil=121.7;
	get_address("CA")->cs_coil=57.3;
	get_address("HA")->cs_coil=4.80;
	get_address("CB")->cs_coil=30.4;
	get_address("C")->cs_coil=175.5;

		
	get_address("H")->cs_mean=8.28;
	get_address("N")->cs_mean=121.58;
	get_address("CA")->cs_mean=57.71;
	get_address("HA")->cs_mean=4.67;
	get_address("CB")->cs_mean=29.98;
	get_address("C")->cs_mean=176.19;

	get_address("H")->cs_wang=7.97;
	get_address("N")->cs_wang=120.99;
	get_address("CA")->cs_wang=57.54;
	get_address("HA")->cs_wang=4.60;
	get_address("CB")->cs_wang=29.60;
	get_address("C")->cs_wang=175.87;

	coil_pre[0]= 0.97;
	coil_pre[1]=-0.46;
	coil_pre[2]= 0.00;
	coil_pre[3]= 0.28;
	coil_pre[4]=-0.08;
	coil_pre[5]=-0.05;


	coil_fol[0]=-0.59;
	coil_fol[1]=-0.33;
	coil_fol[2]=-0.06;
	coil_fol[3]=-0.03;
	coil_fol[4]=-0.10;
	coil_fol[5]=-0.03;

}
CTrp::~CTrp() {}



//TYR
void CTyr::schbond(vector<bbhbond_group> *grp)
{
	bbhbond_group t;
	t.id=residue;
	t.code=OneLetterName;
	t.cpos=-1;
	t.opos=-1;
	t.npos=get("OH").index;
	t.hpos=get("HH").index;
	t.type=12; 

	grp->push_back(t);
	return;
}


struct noeatoms CTyr::query(string name)
{

struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("QB")==0)
	{
		tt.clear();
		tt.push_back(get("HB2").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HB3").index);
		t.atoms.push_back(tt);
		t.length=1.8;
	}
	else if(name.compare("QD")==0)
	{
		tt.clear();
		tt.push_back(get("HD1").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HD2").index);
		t.atoms.push_back(tt);
		t.length=5.0;
	}
	else if(name.compare("QE")==0)
	{
		tt.clear();
		tt.push_back(get("HE1").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HE2").index);
		t.atoms.push_back(tt);
		t.length=5.0;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}


void CTyr::ring(vector<ring_group> *ring)
{
	ring_group t;

	t.x1=2; //type
	t.x2=get("CG").index;
	t.x3=get("CD1").index;
	t.x4=get("CE1").index;
	t.x5=get("CZ").index;
	t.x6=get("CE2").index;
	t.x7=get("CD2").index;
	ring->push_back(t);
	return ;
}

void CTyr::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);
	
	t.type=4;
	t.x1=get("CA").index;
	t.x2=get("CB").index;
	t.x3=get("CG").index;
	t.x4=get("CD1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);
	
	return;
}


CTyr::CTyr() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='Y';strcpy(ThreeLetterName,"TYR");

	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG";t.cs_name="CG";t.base_name="CB";atoms.push_back(t);
	t.name="CD1";t.cs_name="CD";t.base_name="CB";atoms.push_back(t);
	t.name="CD2";t.cs_name="CD";t.base_name="CB";atoms.push_back(t);
	t.name="CE1";t.cs_name="CE";t.base_name="CB";atoms.push_back(t);
	t.name="CE2";t.cs_name="CE";t.base_name="CB";atoms.push_back(t);
	t.name="CZ";t.cs_name="CZ";t.base_name="CB";atoms.push_back(t);
	t.name="OH";t.cs_name="OH";t.base_name="CB";atoms.push_back(t);
	
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=93;t.name="HB2";t.cs_name="HB2";t.base_name="HB";atoms.push_back(t);
	t.proton_type=94;t.name="HB3";t.cs_name="HB3";t.base_name="HB";atoms.push_back(t);

	t.carbon_name="CD";
	t.proton_type=95;t.name="HD1";t.cs_name="HD";t.base_name="HD1";atoms.push_back(t);
	t.proton_type=95;t.name="HD2";t.cs_name="HD";t.base_name="HD2";atoms.push_back(t);
	t.carbon_name="CE";
	t.proton_type=96;t.name="HE1";t.cs_name="HE";t.base_name="HE1";atoms.push_back(t);
	t.proton_type=96;t.name="HE2";t.cs_name="HE";t.base_name="HE2";atoms.push_back(t);
	
	t.carbon_name="OH";
	t.proton_type=97;t.name="HH";t.cs_name="HH";t.base_name="HH";atoms.push_back(t);

	get_address("H")->cs_coil=8.24;
	get_address("N")->cs_coil=120.0;
	get_address("CA")->cs_coil=57.6;
	get_address("HA")->cs_coil=4.73;
	get_address("CB")->cs_coil=39.4;
	get_address("C")->cs_coil=174.8;

		
	get_address("H")->cs_mean=8.30;
	get_address("N")->cs_mean=120.45;
	get_address("CA")->cs_mean=58.20;
	get_address("HA")->cs_mean=4.61;
	get_address("CB")->cs_mean=39.26;
	get_address("C")->cs_mean=175.45;

	get_address("H")->cs_wang=7.99;
	get_address("N")->cs_wang=119.37;
	get_address("CA")->cs_wang=57.64;
	get_address("HA")->cs_wang=4.56;
	get_address("CB")->cs_wang=38.78;
	get_address("C")->cs_wang=175.29;

	coil_pre[0]= 0.46;
	coil_pre[1]=-0.03;
	coil_pre[2]=-0.19;
	coil_pre[3]=-0.10;
	coil_pre[4]= 0.02;
	coil_pre[5]= 0.01;

	coil_fol[0]=-0.48;
	coil_fol[1]=-0.51;
	coil_fol[2]= 0.10;
	coil_fol[3]= 0.00;
	coil_fol[4]=-0.07;
	coil_fol[5]=-0.06;
}
CTyr::~CTyr() {}



//Val

void CVal::methyl_ambig(int flag)
{
	int i;
	for(i=0;i<(int)atoms.size();i++)
	{
		if(atoms.at(i).base_name=="HG1" && atoms.at(i).ambig==2)
		{
			if(atoms.at(i).cs_exp>0.0)
				atoms.at(i).cs_exp=999.0;
		}
		if(atoms.at(i).base_name=="HG2" && atoms.at(i).ambig==2)
		{
			if(atoms.at(i).cs_exp>0.0)
				atoms.at(i).cs_exp=999.0;
		}
	}
	return;
}

struct noeatoms CVal::query(string name)
{
struct noeatoms t;	
	vector<int> tt;

	t=CAminoacid::query(name);

	if(t.atoms.size()>0)
		;
	else if(name.compare("MG1")==0 || name.compare("QG1")==0)
	{
		tt.push_back(get("HG11").index);
		tt.push_back(get("HG12").index);
		tt.push_back(get("HG13").index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	else if(name.compare("MG2")==0 || name.compare("QG2")==0)
	{
		tt.push_back(get("HG21").index);
		tt.push_back(get("HG22").index);
		tt.push_back(get("HG23").index);
		t.atoms.push_back(tt);
		t.length=0.0;
	}
	else if(name.compare("QG")==0 || name.compare("QQG")==0)
	{
		tt.clear();
		tt.push_back(get("HG21").index);
		tt.push_back(get("HG22").index);
		tt.push_back(get("HG23").index);
		t.atoms.push_back(tt);
		tt.clear();
		tt.push_back(get("HG11").index);
		tt.push_back(get("HG12").index);
		tt.push_back(get("HG13").index);
		t.atoms.push_back(tt);
		t.length=5.0;
	}
	else
		cout<<"In "<<residue<<" "<<ThreeLetterName<<" no NOE match for "<<name.c_str()<<endl;
	return t;
}

void CVal::proton2(vector<struct proton> *sel)
{
	struct proton t;


	t.type=5;
	t.name="HG1";
	get_proton(&t);
	sel->push_back(t);

	t.type=6;
	t.name="HG2";
	get_proton(&t);
	sel->push_back(t);

	return;
}





void CVal::dihe(vector<dihe_group> * dihe_index)
{
	dihe_group t;
	t.code=OneLetterName;
	t.id=residue;

	t.type=3;
	t.x1=get("N").index;
	t.x2=get("CA").index;
	t.x3=get("CB").index;
	t.x4=get("CG1").index;
	t.bgood=0;if(t.x1>=0 && t.x2>=0 && t.x3>=0 && t.x4>=0) t.bgood=1; dihe_index->push_back(t);

	return;
}

CVal::CVal() 
{
    //CAminoacid::CAminoacid();
	OneLetterName='V';
	strcpy(ThreeLetterName,"VAL");
	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.cs_coil=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;
	t.bb=0;
	t.proton=0;t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.name="CG1";t.cs_name="CG1";t.base_name="CB";atoms.push_back(t);
	t.name="CG2";t.cs_name="CG2";t.base_name="CB";atoms.push_back(t);
	
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=98;t.name="HB";t.cs_name="HB";t.base_name="HB";atoms.push_back(t);
	t.carbon_name="CG1";
	t.proton_type=5;t.name="HG11";t.cs_name="HG1";t.base_name="HG1";atoms.push_back(t);
	t.proton_type=5;t.name="HG12";t.cs_name="HG1";t.base_name="HG1";atoms.push_back(t);
	t.proton_type=5;t.name="HG13";t.cs_name="HG1";t.base_name="HG1";atoms.push_back(t);
	t.carbon_name="CG2";
	t.proton_type=6;t.name="HG21";t.cs_name="HG2";t.base_name="HG2";atoms.push_back(t);
	t.proton_type=6;t.name="HG22";t.cs_name="HG2";t.base_name="HG2";atoms.push_back(t);
	t.proton_type=6;t.name="HG23";t.cs_name="HG2";t.base_name="HG2";atoms.push_back(t);

	get_address("H")->cs_coil=8.32;
	get_address("N")->cs_coil=121.8;
	get_address("CA")->cs_coil=61.4;
	get_address("HA")->cs_coil=4.30;
	get_address("CB")->cs_coil=32.8;
	get_address("C")->cs_coil=175.1;

		
	get_address("H")->cs_mean=8.28;
	get_address("N")->cs_mean=121.08;
	get_address("CA")->cs_mean=62.57;
	get_address("HA")->cs_mean=4.17;
	get_address("CB")->cs_mean=32.71;
	get_address("C")->cs_mean=175.68;

	get_address("H")->cs_wang=7.98;
	get_address("N")->cs_wang=119.91;
	get_address("CA")->cs_wang=62.00;
	get_address("HA")->cs_wang=4.13;
	get_address("CB")->cs_wang=32.35;
	get_address("C")->cs_wang=175.66;

	coil_pre[0]=2.77;
	coil_pre[1]=-0.00;
	coil_pre[2]=-0.16;
	coil_pre[3]=0.03;
	coil_pre[4]=0.17;
	coil_pre[5]=0.04;

	coil_fol[0]=-0.05;
	coil_fol[1]=-0.08;
	coil_fol[2]= 0.05;
	coil_fol[3]= 0.18;
	coil_fol[4]= 0.01;
	coil_fol[5]= 0.05;

}
CVal::~CVal() {}

//Missing residue
struct noeatoms CMiss::query(string name)
{
	struct noeatoms t;t.length=0.0;vector<int> tt;
	return t;
}

void CMiss::dihe(vector<dihe_group> * dihe_index)
{
	return;
}

void CMiss::process(vector<string> block)
{
	original_residue = 0;
	return;
}


CMiss::CMiss() 
{
    //CAminoacid::CAminoacid();
	resname="Missing";OneLetterName='X';strcpy(ThreeLetterName,"MIS");
	atoms.clear();
	order_parameters.clear();
	original_residue = 0;

	for(int i=0;i<6;i++)
	{
		coil_pre[i]=coil_fol[i]=0.0;
	}
}

CMiss::~CMiss() {}

// unknown residue

struct noeatoms CUnk::query(string name)
{
	struct noeatoms t;t.length=0.0;vector<int> tt;
	return t;
}


void CUnk::dihe(vector<dihe_group> * dihe_index)
{
	return;
}



void CUnk::process(vector<string> block)
{
	int i,j;
	int index;
	char c;
	string t;
	string atomname;
	string part;
	vector<string> subblock;
	struct Atom atom;

	original_residue = 0;

	for(i=0;i<(int)block.size();i++)
	{
		t=block[i];
		if(i==0)
		{
			head=t.substr(0,6);
			resname=t.substr(17,3);
		}
		atomname=t.substr(11,5);  //atom name here
		while((j=atomname.find(" ")) != string::npos)
		{
			atomname.replace(j, 1, "");
		}	
		part=t.substr(6,5); //atom index
		index=atoi(part.c_str());

		c=atomname.at(0);
		if(c>'0' && c<'9')  //if name started with number, move that number to the end
		{
			atomname.erase(0,1);
			atomname.append(1,c);
		}

		c=atomname.at(0);
		if(c=='H' && atomname.size()>1)  //started with H, but it is not HN
		{
			//do nothing for protons!
		}
		else  //heavy atoms or HN atom, just store them
		{
			atom.name=atomname;
			atom.cs_name=atomname;
			atom.base_name=atomname;
			atom.index=index;
			atom.bb=1;
			atom.proton=0;
			atom.cs_exp=atom.cs_pre=999.0;
		}
	}
	return;
}


CUnk::CUnk() 
{
    //CAminoacid::CAminoacid();
	resname="None";OneLetterName='X';strcpy(ThreeLetterName,"UNK");
	order_parameters.clear();

	struct Atom t;
	t.index=-1;t.cs_exp=999.0;t.cs_pre=999.0;t.x=-999.0;t.y=-999.0;t.z=-999.0;t.cs_coil=999.0;
	t.bb=0;
	t.proton=0;t.ambig=0;
	t.name="CB";t.cs_name="CB";t.base_name="CB";atoms.push_back(t);
	t.proton=1;
	t.carbon_name="CB";
	t.proton_type=1;t.name="HB1";t.cs_name="HB";t.base_name="HB";atoms.push_back(t);
	t.proton_type=1;t.name="HB2";t.cs_name="HB";t.base_name="HB";atoms.push_back(t);
	t.proton_type=1;t.name="HB3";t.cs_name="HB";t.base_name="HB";atoms.push_back(t);

	for(int i=0;i<6;i++)
	{
		coil_pre[i]=coil_fol[i]=0.0;
	}
}

CUnk::~CUnk() {}



//water or ligand as a fake residue


CLigand::CLigand() 
{
	atomindexs.clear();
	atomnames.clear();
	id=-1;
	resname="None";
}

CLigand::~CLigand() {}

void CLigand::process(vector<string> block)
{
	int i,j;
	int index;
	string t;
	string atomname;
	string part;
	vector<string> subblock;

	t=block.at(0);
	head=t.substr(0,6);
	resname=t.substr(17,3);

	original_residue=atoi(block.at(0).substr(22,4).c_str());
	
	for(i=0;i<(int)block.size();i++)
	{
		t=block[i];
		atomname=t.substr(11,5);  //atom name here
		while((j=atomname.find(" ")) != string::npos)
		{
			atomname.replace(j, 1, "");
		}
		part=t.substr(6,5); //atom index
		index=atoi(part.c_str());
		atomindexs.push_back(index);
		atomnames.push_back(atomname);
	}
	return;
}

void CLigand::heavycoor(vector<int> *t)
{
	t->insert(t->end(),atomindexs.begin(),atomindexs.end());
}


double CAminoacid::wishart[22][12]={
										{50.86,54.86,21.72,18.27,175.3,179.58,8.59,7.99,4.87,4.03,125.57,121.65},
										{54.63,59.05,32.36,30,175.04,178.11,8.57,8.03,4.85,4,122.6,118.99},
										{52.48,55.67,40.43,38.28,174.55,176.74,8.7,8.2,5.26,4.45,122.7,117.6},
										{53.4,57.04,42.78,40.5,175.15,178.07,8.56,8.05,5.01,4.44,123.82,119.9},
										{54.33,58.61,31.92,28.33,174.5,178.35,8.51,8.11,4.97,4.03,123.14,118.59},
										{55.55,59.3,32.45,29.2,175.01,178.46,8.66,8.32,4.76,3.99,123.52,119.89},
										{45.08,47.02,999,999,173.01,176.31,8.27,8.23,4.09,3.84,110.19,107.34},
										{54.8,59.62,32.2,29.91,173.8,176.83,8.76,8.03,5.07,4.06,121.65,118.09},
										{60,64.68,40.09,37.59,174.79,177.49,8.74,8.06,4.72,3.66,124.12,120.22},
										{53.94,57.54,44.02,41.4,175.16,178.42,8.63,8.02,4.85,4,125.69,120.18},
										{55.01,59.11,34.86,32.31,174.93,177.79,8.54,8.04,4.96,3.98,123.29,119.9},
										{54.1,58.45,34.34,31.7,174.64,177.76,8.43,8.05,4.94,4.03,121.67,118.69},
										{56.33,60.74,41.64,38.91,174.15,176.42,8.8,8.21,5.1,4.11,121.95,119.12},
										{62.79,65.52,32.45,31.08,176.41,178.34,999,999,4.72,4.13,999,999},
										{57.14,60.86,65.39,62.81,173.52,176.51,8.57,8.11,5.08,4.2,117.44,114.78},
										{61.1,65.89,70.82,68.64,173.47,176.62,8.5,8.1,4.81,4.02,118.09,115.3},
										{56.28,60.03,31.78,28.74,175.1,177.81,8.83,8.24,5.24,4.35,124.04,120.48},
										{56.56,61.07,40.79,38.38,174.65,177.05,8.69,8.1,5,4.14,122.55,119.67},
										{60.72,65.96,33.81,31.41,174.66,177.75,8.73,7.99,4.66,3.5,123.27,119.53},
										{57.64,62.86,29.48,26.99,173.86,177.42,9,8.22,5.18,4.16,123.27,117.4},
										{54.19,58.57,43.79,40.02,172.73,176.84,8.68,8.58,5.21,4.53,121.81,119.51},
										{999,999,999,999,999,999,999,999,999,999,999,999}
									};


double CAminoacid::wang_rc[21][6]={
{120.58, 175.52, 60.79, 38.43, 7.94, 4.18},
{119.91, 175.66, 62.00, 32.35, 7.98, 4.13},
{120.37, 176.00, 54.00, 40.78, 8.31, 4.62},
{118.50, 174.84, 53.00, 38.43, 8.35, 4.66},
{119.72, 175.46, 57.46, 39.41, 8.09, 4.59},
{118.92, 174.78, 55.74, 29.50, 8.18, 4.60},
{120.99, 175.87, 57.54, 29.60, 7.97, 4.60},
{119.37, 175.29, 57.64, 38.78, 7.99, 4.56},
{121.10, 176.15, 56.29, 32.53, 8.17, 4.28},
{121.57, 176.70, 54.77, 42.14, 8.06, 4.36},
{120.14, 175.94, 55.43, 32.92, 8.22, 4.47},
{119.82, 175.75, 55.89, 29.01, 8.20, 4.29},
{120.75, 176.01, 56.18, 30.36, 8.21, 4.26},
{120.62, 176.32, 56.66, 29.87, 8.36, 4.28},
{113.88, 174.78, 61.30, 68.92, 8.16, 4.44},
{118.10, 175.11, 58.24, 29.54, 8.10, 4.59},
{116.00, 174.41, 58.20, 63.75, 8.22, 4.45},
{123.82, 177.28, 52.46, 18.98, 8.09, 4.31},
{109.48, 174.01, 45.28, 0, 8.37, 3.97},
{0, 176.62, 63.24, 31.81, 0, 4.41},
{999,999,999,999,999,999}
};

double CAminoacid::wang_c1[21][6]={
{-2.21, 0.14, -0.01, -0.04, -0.07, -0.01},
{1.36, -0.00, 0.44, -0.18, 0.03, -0.00},
{-0.43, 0.07, 0.20, -0.07, -0.01, -0.03},
{-0.36, 0.04, 0.01, 0.04, -0.01, -0.01},
{0.20, -0.14, 0.02, -0.06, 0.01, -0.01},
{-0.43, -0.09, -0.26, 0.17, -0.10, 0.04},
{-0.05, -0.13, 0.16, 0.07, -0.01, -0.05},
{2.92, 0.11, -0.15, 0.26, 0.12, 0.04},
{-0.26, -0.09, -0.07, 0.08, -0.02, -0.00},
{-0.76, 0.13, -0.07, -0.03, -0.06, 0.00},
{0.69, 0.10, 0.09, 0.10, 0.05, 0.05},
{-0.76, -0.19, 0.18, -0.20, 0.01, -0.03},
{-0.94, 0.21, 0.07, -0.16, 0.14, -0.04},
{-0.09, 0.07, 0.13, 0.06, 0.01, -0.00},
{-0.45, -0.07, 0.00, 0.01, -0.03, 0.01},
{1.16, -0.10, 0.11, -0.06, 0.02, 0.01},
{1.23, -0.07, 0.05, -0.11, 0.02, -0.00},
{2.77, -0.00, -0.16, 0.03, 0.17, 0.04},
{0.97, -0.46, 0.00, 0.28, -0.08, -0.05},
{0.46, -0.03, -0.19, -0.10, 0.02, 0.01},
{0,0,0,0,0,0}
};

double CAminoacid::wang_c2[21][6]={
{-0.11, 0.05, 0.07, -0.09, -0.01, -0.03},
{-1.17, 0.10, 0.17, 0.21, 0.01, 0.03},
{0.23, -0.11, 0.28, 0.11, 0.04, -0.03},
{0.26, 0.14, 0.25, 0.06, 0.07, -0.04},
{-0.35, -0.22, -0.04, -0.09, -0.04, -0.02},
{0.13, 0.47, 0.12, -0.07, 0.06, -0.03},
{-0.09, -0.05, 0.22, -0.24, -0.01, -0.09},
{-0.20, -0.09, 0.03, 0.28, -0.01, 0.03},
{-0.13, -0.13, 0.08, 0.01, -0.04, -0.03},
{-0.49, 0.06, 0.10, -0.10, -0.07, -0.02},
{-0.02, 0.19, 0.22, 0.06, 0.01, -0.01},
{-0.03, -0.23, 0.24, -0.06, 0.04, -0.04},
{0.92, -1.19, -2.04, -0.20, -0.17, 0.21},
{-0.31, 0.10, 0.37, -0.13, 0.01, -0.05},
{-0.09, 0.16, 0.19, -0.12, -0.01, -0.02},
{0.30, 0.10, 0.10, 0.14, 0.07, 0.02},
{0.22, 0.10, 0.02, 0.14, 0.04, 0.08},
{-0.05, -0.08, 0.05, 0.18, 0.01, 0.05},
{-0.59, -0.33, -0.06, -0.03, -0.10, -0.03},
{-0.48, -0.51, 0.10, 0.00, -0.07, -0.06},
{0,0,0,0,0,0}
};
