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
#include "lmmin.h"
#include "config.h"
#include "supply.h"
#include "bmrb.h"
#include  "aa.h"
#include "pdb.h"
#include "traj.h"
#include "mainbody.h"



// ///////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv)

{
	int begin,stop;
	bool bh,bdetail,bnew,bann,btest,bnew_pdb,blinear,bold;
	class CMainbody mainbody;
	string pdbname;
	string spartaname;
	string gmxname;
	string bmrbname;
	int para;
	int nconf;

	cout<<"PPM: an enemble-based chemical shifts predictor"<<endl;

	bh=bdetail=bnew=btest=bnew_pdb=bann=blinear=bold=0;

	CCommandline cmdline;
	vector<string> args,args2,args3;

	args.push_back("-h");		args2.push_back("no");						args3.push_back("print help informaiton then quit");
	args.push_back("-mode");	args2.push_back("ann");						args3.push_back("prediciton algorithm: ann(default) or linear");
	args.push_back("-pdb");		args2.push_back("pdb.pdb");					args3.push_back("input pdb file name");
	args.push_back("-bmrb");	args2.push_back("bmrb.dat");				args3.push_back("input experimental chemical shifts file in NMRSTAR format");
    args.push_back("-pre");		args2.push_back("bmrb_pre.dat");			args3.push_back("output filename for predicted chemical shifts in NMRSTAR format");
	args.push_back("-begin");	args2.push_back("0");						args3.push_back("Index of first snapshot to be used (start from 0)");
	args.push_back("-stop");	args2.push_back("0");						args3.push_back("Index of last snapshot to be used (0 means last snapshot)");
	args.push_back("-para");	args2.push_back("pdb");                     args3.push_back("Parameter set: pdb(ppm_one) or old(ppm)");
    args.push_back("-detail(s)");	args2.push_back("no");					args3.push_back("calculate chemical shifts for each snapshots, only applicable to ensemble predictions");

		
	cmdline.init(args,args2,args3);
	cmdline.pharse(argc,argv);


	spartaname="sparta.pdb";
	pdbname=cmdline.query("-pdb");
	bmrbname=cmdline.query("-bmrb");
	begin=atoi(cmdline.query("-begin").c_str());
	stop=atoi(cmdline.query("-stop").c_str());
    


	if(cmdline.query("-h").compare("yes")==0)		bh=1;
	if(cmdline.query("-mode").compare("linear")==0)		{blinear=1;bann=0;}
	if(cmdline.query("-mode").compare("ann")==0)		{blinear=0;bann=1;}


	if(cmdline.query("-para").compare("pdb")==0)		para=2;
	if(cmdline.query("-para").compare("old")==0)		bold=1;
    
    if(cmdline.query("-detail").compare("yes")==0)		bdetail=1;
    if(cmdline.query("-details").compare("yes")==0)		bdetail=1;



	cmdline.print();

	if(bh==1)
	{	
		exit(0);
	}


	nconf=mainbody.loadpdb(pdbname,gmxname);
	mainbody.set_range(begin,stop);
	mainbody.load(bmrbname);




	if(bold)  //previous generation PPM
	{
		cout<<"Prediction using the old PPM parameters\n";
		cout<<"Chemical shifts root-mean-square deviations (RMSDs) between predicted and experimental values:"<<endl;
        if(bdetail)
        {
            mainbody.predict_proton2();
            mainbody.predict_bb2();
        }
        else
        {
            mainbody.predict_proton();
            mainbody.predict_bb();
        }
		mainbody.print_prediction(cmdline.query("-pre"));
	}
	else if(blinear && para==2)
	{
		cout<<"Prediction using the linear model with static parameters set\n";
		cout<<"Chemical shifts root-mean-square deviations (RMSDs) between predicted and experimental values:"<<endl;
		mainbody.predict_bb_static_new();
		mainbody.predict_proton_static_new();
		mainbody.print_prediction(cmdline.query("-pre"));
	}


	else if(bann && para==2)
	{
		cout<<"Prediction using the ANN model with static parameters set\n";
		cout<<"Chemical shifts root-mean-square deviations (RMSDs) between predicted and experimental values:"<<endl;
		mainbody.predict_bb_static_ann();
		mainbody.predict_proton_static_new();
		mainbody.print_prediction(cmdline.query("-pre"));
	}

	else
	{
		cout<<"Unrecognized command line arguments!\n";
	}

	return 0;
}









