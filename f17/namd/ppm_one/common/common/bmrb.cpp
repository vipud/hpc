#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <algorithm>

using namespace std;


#include "bmrb.h"

CBmrb::CBmrb()
{
	data.clear();
};

CBmrb::~CBmrb()
{};

void CBmrb::clear()
{
	data.clear();
}



vector<string> CBmrb::loadpdbpart(string bmrbname)
{
	int i;
	ifstream fin(bmrbname.c_str());
	string line;
	vector<string> part;
	string p;
	istringstream iss;
	int flag1,flag2,flag3,flag4;
	vector<string> codes;
	size_t pos1,pos2;
	string databasename,code,title;
	double percent,identity;
	int nres;

	flag1=flag2=flag3=flag4=0;
	physical_state="";
	oligomer_state="";
	component.clear();

	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		part.clear();
		while(iss>>p)
		{
			part.push_back(p);
		}

		if(part.size()<=0)  //blank line.
			continue; 

		if(part.at(0)=="_System_physical_state")
		{
			for(i=1;i<part.size();i++)
			{
				physical_state.append(part.at(i));
				physical_state.append(" ");
			}
		}
		if(part.at(0)=="_System_oligomer_state")
		{
			for(i=1;i<part.size();i++)
			{
				oligomer_state.append(part.at(i));
				oligomer_state.append(" ");
			}		
		}


		else if(part.size()==1)  //head informations
		{
			p=part.at(0);
			if(p.compare("stop_")==0)
			{
				flag1=0;
				flag2=0;
				flag3=0;
				flag4=0;
				continue;
			}
			else if(p.compare("loop_")==0)
			{
				flag1=1;
			}

			else if(flag1==1 && p.compare("_Sequence_homology_expectation_value")==0)
			{
				flag2=1;
			}

			else if(flag1==1 && p.compare("_Mol_system_component_name")==0)
				flag3=1;

			else if(flag1==1 && flag3==1 && p.compare("_Mol_label")==0)
				flag4=1;
		}

		else if(flag1==1 && flag3==1 && flag4==1)
		{
			component.push_back(line);
		}


		else if(flag1==1 && flag2==1)
		{
			if(part.at(0)=="PDB" && part.size()>=7 )
			{
				pos1=line.find("\"");
				pos2=line.find("\"",pos1+1);
				if(pos1!=string::npos && pos2!=string::npos)
				{
					//cout<<line<<endl;
					//cout<<pos1<<" "<<pos2<<endl;
					title=line.substr(pos1,pos2-pos1);
					line.erase(pos1,pos2-pos1+1);
					iss.clear();
					iss.str(line);
					part.clear();
					while(iss>>p)
					{
						part.push_back(p);
					}
				}
				
				databasename=part.at(0);
				code=part.at(1);
				percent=atof(part.at(2).c_str());
				nres=atoi(part.at(3).c_str());
				identity=atof(part.at(4).c_str());

				if(databasename.compare("PDB")==0 && percent>=90.00 && identity>=90.00)
				{
					codes.push_back(code);
				}
				
			}
		}

	}

	return codes;
}

void CBmrb::loaddetail(string name)
{
	int i,j;
	ifstream fin(name.c_str());
	string line,p;
	istringstream iss;
	vector< vector<double> > x;
	vector<double> y,t;

	if(getline(fin,line))
	{
		iss.clear();
		iss.str(line);

		while(iss>>p)
		{
			y.push_back(atof(p.c_str()));
		}
		x.push_back(y);
	
		while(getline(fin,line))
		{
			iss.clear();
			iss.str(line);

			y.clear();
			while(iss>>p)
			{
				y.push_back(atof(p.c_str()));
			}
			if(y.size()==x.at(0).size()) x.push_back(y);
			else break;
		}
	}


	for(i=0;i<x.at(0).size();i++)
	{
		t.clear();
		for(j=0;j<x.size();j++)
		{
			t.push_back(x.at(j).at(i));
		}
		detail.push_back(t);
	}

	return;
}
		

void CBmrb::attach_bmrb(CBmrb bmrb)
{
	string s1,s2;
	vector<int> out;
	double dtotal;
	int total;
	int i,j,n,m;
	vector<struct CBmrbline> empty;

	s1=getseq();
	s2=bmrb.getseq();

	out=Sequence::align(s1,s2);

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
	for(i=0;i<(int)data.size();i++)
	{
		for(j=0;j<bmrb.getsize();j++)
		{
			if(bmrb.getdata(j).res-total==data.at(i).res)
			{
				data.at(i).others.push_back(bmrb.getdata(j).block);
			}
		}

		if(data.at(i).others.size()==0)
			data.at(i).others.push_back(empty);
	}



	return;
}

void CBmrb::print(string name)
{
	int i,j,k,kk,index;
	string atom;
	FILE *fp=fopen(name.c_str(),"wt");

	index=0;
	for(i=0;i<data.size();i++)
	{	
		
		for(j=0;j<data.at(i).block.size();j++)
		{
			fprintf(fp,"%10d %10s ",data.at(i).res,data.at(i).name.c_str());
			atom=data.at(i).block.at(j).type;

			kk=-1;
			for(k=0;k<data.at(i).others.at(0).size();k++)
			{
				if((data.at(i).others.at(0).at(k).type).compare(atom)==0)
				{
					kk=k;
				}
			}
			if(kk==-1)
				fprintf(fp,"%s 999.0",atom.c_str()); 
			else
				fprintf(fp,"%s %f",atom.c_str(),data.at(i).others.at(0).at(kk).cs);

			for(k=0;k<detail.at(index).size();k++)
				fprintf(fp," %f",detail.at(index).at(k));
			fprintf(fp,"\n");
			index++;
		}
	}

	return;
}


void CBmrb::process(string bmrbname)
{
	ifstream fin(bmrbname.c_str());
	string line;
	vector<string> part;
	string p;
	int oldres;
	int n,index_atom,index_res,index_name,index_type,index_cs,index_ambig;
	int flag1,flag2,flag3;

	istringstream iss;
	int first=1;
	struct CBmrbline bmrbline;
	struct CBmrbdata bmrbdata;

	bmrbline.atom=-1;

	index_atom=index_res=index_name=index_type=index_cs=index_ambig=-1;
	flag1=flag2=flag3=0;
	first=1;
	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		part.clear();
		while(iss>>p)
		{
			part.push_back(p);
		}
		if(part.size()<=0)  //blank line.
			continue; 
		else if(part.size()==1)  //head informations
		{
			flag3=1;
			p=part.at(0);
			if(p.compare("stop_")==0)
			{
				flag1=0;
				n=0;
				continue;
			}
			else if(p.compare("loop_")==0)
			{
				flag1=1;
			}

			else if(flag1==1 && p.compare("_Atom_shift_assign_ID")==0)
			{
				index_atom=1;
				flag2=1;
				n=1;
			}
			else if(flag1==1 && p.compare("_Atom_chem_shift.ID")==0)
			{
				index_atom=1;
				flag2=1;
				n=1;
			}
			else if(flag1==1 && flag2==1)
			{
				n++;
				if(p.compare("_Residue_seq_code")==0)
					index_res=n;
				if(p.compare("_Residue_label")==0)
					index_name=n;
				if(p.compare("_Atom_name")==0)
					index_type=n;
				if(p.compare("_Chem_shift_value")==0)
					index_cs=n;
				if(p.compare("_Chem_shift_ambiguity_code")==0)
					index_ambig=n;


				if(p.compare("_Atom_chem_shift.Seq_ID")==0)
					index_res=n;
				if(p.compare("_Atom_chem_shift.Comp_ID")==0)
					index_name=n;
				if(p.compare("_Atom_chem_shift.Atom_ID")==0)
					index_type=n;
				if(p.compare("_Atom_chem_shift.Val")==0)
					index_cs=n;
				if(p.compare("_Atom_chem_shift.Ambiguity_code")==0)
					index_ambig=n;
			}

		}

		else if(flag1==1 && flag2==1)
		{
			if((int)part.size()>=index_ambig && index_atom!=-1 && index_res!=-1 && index_name!=-1 && index_type!=-1 && index_cs!=-1 && index_ambig!=-1 )
			{
				bmrbline.atom=atoi(part.at(index_atom-1).c_str());
				bmrbline.res=atoi(part.at(index_res-1).c_str());
				bmrbline.name=part.at(index_name-1);
				bmrbline.type=part.at(index_type-1);
				bmrbline.cs=atof(part.at(index_cs-1).c_str());
				bmrbline.ambig=atoi(part.at(index_ambig-1).c_str());
			}
		}


		else if(flag1+flag2+flag3==0) //read nothing, so this is a raw file
		{
			if(part.size()==8 ) // no head information v2.1
			{
				bmrbline.atom=atoi(part.at(0).c_str());
				bmrbline.res=atoi(part.at(1).c_str());
				bmrbline.name=part.at(2);
				bmrbline.type=part.at(3);
				bmrbline.cs=atof(part.at(5).c_str());
				bmrbline.ambig=atoi(part.at(7).c_str());
			}
			else if(part.size()==9 ) // no head information v2.1
			{
				bmrbline.atom=atoi(part.at(1).c_str());
				bmrbline.res=atoi( part.at(2).c_str() );
				bmrbline.name=part.at(3);
				bmrbline.type=part.at(4);
				bmrbline.cs=atof(part.at(6).c_str());
				bmrbline.ambig=atoi(part.at(8).c_str());
			}
		}
		//finish processing one line
		if(bmrbline.type.compare("HN")==0)
			bmrbline.type="H"; 

		if( (flag1==1 && flag2==1 && bmrbline.atom!=-1) || (flag1+flag2+flag3==0) )
		{
			if(first==1)
			{
				first=0;
				oldres=bmrbline.res;
			}
			if(bmrbline.res==oldres)
			{
				bmrbdata.block.push_back(bmrbline);
				bmrbdata.res=bmrbline.res;
				bmrbdata.name=bmrbline.name;
			}
			else if(bmrbline.res> oldres)
			{
				oldres=bmrbline.res;	
				data.push_back(bmrbdata);
				bmrbdata.block.clear();
				bmrbdata.block.push_back(bmrbline);
				bmrbdata.res=bmrbline.res;
				bmrbdata.name=bmrbline.name;
			}
		}
	}

	if(bmrbdata.block.size()>0)
	{
		data.push_back(bmrbdata);
	}

	adjust_first();
		
	return ;
}

string CBmrb::getseq()
{
	int i;
	char c;
	string s;
	int id;
	int base;


	if(data.size()>0)
		base=data.at(0).res-1;
	else
		base=0;

	for(i=0;i<(int)data.size();i++)
	{
		if(data.at(i).name.length()==1)
			c=data.at(i).name.at(0);
		else
		{
			c=Sequence::name2code(data.at(i).name);
		}
		id=data.at(i).res-base-1;
		if(id>=(int)s.size())
			s.resize(id+1,'X');
		s.at(id)=c;

	}
	return s;
}

string CBmrb::getseq(vector<int> &index)
{
	int i;
	char c;
	string s;
	

	for(i=0;i<(int)data.size();i++)
	{
		if(data.at(i).name.length()==1)
			c=data.at(i).name.at(0);
		else
		{
			c=Sequence::name2code(data.at(i).name);
		}
		s.push_back(c);
		index.push_back(data.at(i).res);

	}
	return s;
}


struct CBmrbdata CBmrb::getdata(int i)
{
	return data.at(i);
}

int CBmrb::run_shiftx(string name,int i)
{
	char buffer[100];
	string line;
	vector<string> part;
	string p;
	istringstream iss;
	int flag1,flag2,first,count;

	struct CBmrbline bmrbline;
	struct CBmrbdata bmrbdata;

	data.clear();

	//sprintf(buffer,"~/scratch/programs/shiftx/shiftx 1 %s shiftx%d.dat",name.c_str(),i);
	sprintf(buffer,"shiftx 1 %s shiftx%d.dat",name.c_str(),i);
	system(buffer);
	sprintf(buffer,"shiftx%d.dat",i);
	ifstream fin(buffer);

	flag1=flag2=0;
	first=1;
	count=0;
	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		part.clear();
		while(iss>>p)
		{
			part.push_back(p);
		}

		if(part.size()==0)
			continue;
		if(flag1==1 && part.at(0).compare("NUM")==0)
			break;


		if(part.at(0).compare("NUM")==0)
		{
			flag1=1;
			continue;
		}
		if(part.at(0).compare("---")==0 && flag1==1)
		{
			flag2=1;
			continue;
		}

		if(flag1==1 && flag2==1)
		{
			bmrbdata.block.clear();

			bmrbline.res=atoi(line.substr(0,5).c_str());
			bmrbline.name=Sequence::code2name(line.at(8));
			bmrbline.ambig=0;

			bmrbline.atom=++count;
			bmrbline.type="HA";
			if(bmrbline.name=="GLY") bmrbline.type="HA2";
			bmrbline.cs=atof(line.substr(10,7).c_str());
			if(bmrbline.cs>0.0) bmrbdata.block.push_back(bmrbline);

			bmrbline.atom=++count;
			bmrbline.type="H";
			bmrbline.cs=atof(line.substr(17,7).c_str());
			if(bmrbline.cs>0.0) bmrbdata.block.push_back(bmrbline);

			bmrbline.atom=++count;
			bmrbline.type="N";
			bmrbline.cs=atof(line.substr(24,9).c_str());
			if(bmrbline.cs>0.0) bmrbdata.block.push_back(bmrbline);

			bmrbline.atom=++count;
			bmrbline.type="CA";
			bmrbline.cs=atof(line.substr(33,8).c_str());
			if(bmrbline.cs>0.0) bmrbdata.block.push_back(bmrbline);

			bmrbline.atom=++count;
			bmrbline.type="CB";
			bmrbline.cs=atof(line.substr(41,8).c_str());
			if(bmrbline.cs>0.0) bmrbdata.block.push_back(bmrbline);

			bmrbline.atom=++count;
			bmrbline.type="C";
			bmrbline.cs=atof(line.substr(49,9).c_str());
			if(bmrbline.cs>0.0) bmrbdata.block.push_back(bmrbline);

			bmrbdata.res=bmrbline.res;
			bmrbdata.name=bmrbline.name;
			data.push_back(bmrbdata);
		}
	}
	adjust_first();
	return 0;
}


int CBmrb::run_shifts(string name,int i)
{
	char buffer[100];
	string line;
	vector<string> part;
	string p;
	istringstream iss;
	int flag1,flag2,first,count;
	int oldres;

	struct CBmrbline bmrbline;
	struct CBmrbdata bmrbdata;

	data.clear();

	//sprintf(buffer,"export SHIFTSHOME=~/scratch/programs/shifts-4.3");
	//system(buffer);
	sprintf(buffer,"cp %s t.pdb",name.c_str());
	system(buffer);
	//sprintf(buffer,"~/scratch/programs/shifts-4.3/bin/shifts  -qdb t");
	sprintf(buffer,"shifts  -qdb t");
	system(buffer);
	sprintf(buffer,"mv t.qdb shifts%d.qdb",i);
	system(buffer);
	sprintf(buffer,"shifts%d.qdb",i);
	ifstream fin(buffer);

	flag1=flag2=0;
	first=1;
	count=0;
	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		part.clear();
		while(iss>>p)
		{
			part.push_back(p);
		}

		if(part.size()==0)
			continue;

		if(part.at(0).compare("#PDB")==0)
		{
			flag1=1;
			continue;
		}
		if(part.at(0).find("#----")!=string::npos && flag1==1)
		{
			flag2=1;
			continue;
		}

		if(flag1==1 && flag2==1)
		{
			count++;
			bmrbline.atom=count;
			bmrbline.res=atoi(line.substr(12,4).c_str());
			bmrbline.name=line.substr(7,3);
			bmrbline.cs=atof(line.substr(78,7).c_str());
			bmrbline.ambig=0;
			bmrbline.type=line.substr(16,2);

			bmrbline.type.erase(remove(bmrbline.type.begin(), bmrbline.type.end(), ' '), bmrbline.type.end());
			transform(bmrbline.type.begin(), bmrbline.type.end(),bmrbline.type.begin(), ::toupper);


			if(bmrbline.type.compare("HN")==0)
				bmrbline.type="H";
			if(bmrbline.type.compare("CO")==0)
				bmrbline.type="C";
		
			if(first==1)
			{
				first=0;
				oldres=bmrbline.res;
			}
			if(bmrbline.res==oldres)
			{
				bmrbdata.block.push_back(bmrbline);
				bmrbdata.res=bmrbline.res;
				bmrbdata.name=bmrbline.name;
			}
			else if(bmrbline.res> oldres)
			{
				oldres=bmrbline.res;	
				data.push_back(bmrbdata);
				bmrbdata.block.clear();
				bmrbdata.block.push_back(bmrbline);
				bmrbdata.res=bmrbline.res;
				bmrbdata.name=bmrbline.name;
			}
		}
	}

	if(bmrbdata.block.size()>0)
	{
		data.push_back(bmrbdata);
	}

	adjust_first();
		
	return 0;
}

int CBmrb::run_sparta(string name,int i)
{
	char buffer[100];
	string line;
	vector<string> part;
	string p;
	istringstream iss;
	int flag1,flag2,first,count;
	int oldres;

	struct CBmrbline bmrbline;
	struct CBmrbdata bmrbdata;

	data.clear();
	
	sprintf(buffer,"sparta+ -in %s > /dev/null 2>/dev/null ",name.c_str());
	system(buffer);
	sprintf(buffer,"mv pred.tab pred%d.pdb",i);
	system(buffer);
	sprintf(buffer,"pred%d.pdb",i);
	ifstream fin(buffer);

	flag1=flag2=0;
	first=1;
	count=0;
	while(getline(fin,line))
	{
		iss.clear();
		iss.str(line);
		part.clear();
		while(iss>>p)
		{
			part.push_back(p);
		}

		if(part.size()==0)
			continue;

		if(part.at(0).compare("VARS")==0)
		{
			flag1=1;
			continue;
		}
		if(part.at(0).compare("FORMAT")==0)
		{
			flag2=1;
			continue;
		}

		if(flag1==1 && flag2==1)
		{
			count++;
			bmrbline.atom=count;
			bmrbline.res=atoi(line.substr(0,4).c_str());
			bmrbline.name=Sequence::code2name(line.at(8));
			bmrbline.cs=atof(line.substr(24,9).c_str());
			bmrbline.ambig=0;

			line=line.substr(10,4);
			iss.clear();
			iss.str(line);
			iss>>bmrbline.type;

			if(bmrbline.type.compare("HN")==0)
				bmrbline.type="H";
		
			if(first==1)
			{
				first=0;
				oldres=bmrbline.res;
			}
			if(bmrbline.res==oldres)
			{
				bmrbdata.block.push_back(bmrbline);
				bmrbdata.res=bmrbline.res;
				bmrbdata.name=bmrbline.name;
			}
			else if(bmrbline.res> oldres)
			{
				oldres=bmrbline.res;	
				data.push_back(bmrbdata);
				bmrbdata.block.clear();
				bmrbdata.block.push_back(bmrbline);
				bmrbdata.res=bmrbline.res;
				bmrbdata.name=bmrbline.name;
			}
		}
	}

	if(bmrbdata.block.size()>0)
	{
		data.push_back(bmrbdata);
	}

	adjust_first();
		
	return 0;
}

void CBmrb::adjust_first()
{
	int i;
	int base;

	if(data.size()>0)
	{
		base=data.at(0).res-1;
		for(i=0;i<(int)data.size();i++)
			data.at(i).res-=base;
	}
	
	return;
}