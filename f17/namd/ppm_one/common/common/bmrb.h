#include "supply.h"


#ifndef BMRBFILE
#define BMRBFILE


struct CBmrbline
{
	int atom;
	int res;
	string name;
	string type;
	double cs;
	int ambig;
	bool used;
};

struct CBmrbdata
{
	vector<struct CBmrbline> block;
	int res;
	string name;
	vector< vector<struct CBmrbline> > others;
};

class CBmrb;

class CBmrb 
{
private:
protected:
		vector<CBmrbdata> data;
		string physical_state;
		string oligomer_state;
		vector<string> component;
		vector< vector<double> > detail;
		void adjust_first();
public:
		CBmrb(void);
		~CBmrb(void);
		void clear(void);
		void process(string);
		void loaddetail(string);
		string getseq();
		string getseq(vector<int> &index);
		struct CBmrbdata getdata(int);
		vector<string> loadpdbpart(string bmrbname);
		void attach_bmrb(class CBmrb);
		void print(string);

		//for CS prediction programs
		int run_sparta(string name,int i);
		int run_shifts(string name,int i);
		int run_shiftx(string name,int i);


		inline int getsize() {return data.size();};
		inline string get_physical() {return physical_state;};
		inline string get_oligomer() {return oligomer_state;};
		inline vector<string> get_component() {return component;};
		
};

#endif