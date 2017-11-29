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

#include <unordered_map>



int main(int argc, char ** argv)

{

	cout << "PPM_ONE Chemical Shift Prediction" << endl;
	#ifdef _OPENACC
	cout << "With OpenACC" << endl;
	#endif

	class CMainbody mainbody;
	int begin,stop;
	string pdbname;
	string bmrbname;
	string pre;

	unordered_map<string, string> args;
	for(int i = 0; i < argc; i++)
	{
		if(argv[i][0] == '-')
		{
			string s1 = string(argv[i]);
			string s2;
			if(s1 == "-h"){
				s2 = " ";
			} else {
				s2 = string(argv[i+1]);
			}
			args.insert(pair<string, string>(s1, s2));
		}
	}
	
	pdbname = "pdb.pdb";
	bmrbname = "bmrb.dat";
	pre = "bmrb_pre.dat";
	begin = 0;
	stop = 0;

	if(args.find("-pdb") != args.end())
		pdbname = args["-pdb"];
	if(args.find("-bmrb") != args.end())
		bmrbname = args["-bmrb"];
	if(args.find("-begin") != args.end())
		begin = atoi(args["-begin"].c_str());
	if(args.find("-stop") != args.end())
		stop = atoi(args["-stop"].c_str());
	if(args.find("-h") != args.end())
		return 0;

	mainbody.loadpdb(pdbname);
	mainbody.set_range(begin,stop);
	mainbody.load(bmrbname);

	mainbody.predict_bb_static_ann();
	mainbody.predict_proton_static_new();
	mainbody.print_prediction(pre);

	return 0;

}









