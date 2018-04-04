// ///////////////////////////////////////////////////////////////////////////

#include "supply.h"
#include "bmrb.h"
#include  "aa.h"
#include <omp.h>


#ifndef PDB
#define PDB


struct dsspline
{
	int id;
	char code;
	char ss;
};

class CDssp
{
private:
		
protected:
public:
		vector<struct dsspline> data;
		CDssp(void);
		~CDssp(void);
		int loaddata(string);
		string getseq();
};

struct proteinblock
{
	vector<string> block;
	int iligand;
	string residue;
	int index;
};



class CPdb 
{
private:
		vector< struct proteinblock > blocks;
		vector<int> chains;
		vector<int> chain_block;
		vector<int> chain_ligand;
		vector<CAminoacid *> v;
		char *v_oneletternames;
		vector<CLigand *> ligand;
		string pdbfilename;
		string pdbseq;
		string headseq;
		int natom;
		int natom2;
		int nmiss;

		vector<double> x,y,z;  //coordinate
		vector<int> heavy,boundary,allcoor;
		vector<double> b,b2; //b-factor
		vector<string> atomname;

protected:
	void setup(CAminoacid *t, CLigand *tt, int iligand, string residue , int index_old,vector<string> block);


public:
		CAminoacid **v_arr;
		int v_size;
	//dssp stuff
		CDssp dssp;
	//nmr constrain stuff
		int loadnoedata(string,vector<struct noeline> *);
		int loadnoedata(string,string,vector <struct noeline> *);
		int actualload(vector<string> *, vector<string> *, vector <struct noeline> *);
		void outputnoe(string,vector <struct noeline> *);
		void inputnoe(string,vector <struct noeline> *);

		void print_gromacs(string,vector <struct diheline> *);
		void print_gromacs(string,vector <struct noeline> *);
		void print_print(FILE *,vector<int>,vector<int>,int,int,double,double,double,double);
		
		
		int loaddihecons(string ,vector<struct diheline> *);
		int dihecons_actualload(vector<string>*,vector<struct diheline> *);
		int load_exactnoe(string filename,vector <struct noeline> *nmrcons);

		//main stuff
		CPdb(void);
		~CPdb(void);
		void clear();
		int buildpdb(string);
		int loadpdb(string); // Updated for OpenACC
		int loadpdb_old(string);
		void getca(vector<int> *);
		vector<int> getselectca(vector<int>);
		void getdihe(vector<dihe_group> *, vector<int> *);
		void getdihe(vector<dihe_group> *);
		void getbbdihe(vector<dihe_group> *);
		void getbbdihe(vector<dihe_group> *, int *);
		void getbbdihe_nopro(vector<dihe_group> *, int *);
		void getring(vector<ring_group> *);
		void proton(vector<struct proton> *);
		void proton_acc(vector<struct proton> *); // Optimized Function
		void proton(vector<struct proton> *,int);
		void allproton(vector<struct proton> *);
		void allproton_acc(vector<struct proton> *); // Optimized Function
		void allproton3(vector<struct proton> *);
		void allproton3_acc(vector<struct proton> *); // Optimized Function
		void ani(vector<struct ani_group> *);
		void ani_acc(vector<struct ani_group> *); // Optimized Function
		void getbb(vector<struct bb_group> *);
        void getbb_assign(vector<struct bb_group> *);
		void bbhbond(vector<bbhbond_group> *);
		void schbond(vector<bbhbond_group> *);
		void bbnh(vector<nh_group> *);
		void bbco(vector<co_group> *);
		void ired(vector<struct ired> *);
		void loadred(vector<struct ired> *);
		void getred(int,vector<struct ired> *);
		void loadred(string);
		void clearred(void);
		void caha(vector<index_three> *);
		void heavycoor();
		char code(int);
		int  chain(int in);
		void name(int,char*);
		void print_prediction();
		void print_prediction(string);
////////////////////////////////////////////////////////////
		void print_debug(string);
////////////////////////////////////////////////////////////
		void attach_bbprediction(int,double*);
		void attach_bbprediction(int, double, double, double, double, double, double);
		void attach_protonprediction(int,string,double);
		int attach_bmrb(class CBmrb );
		double test_bmbr(class CBmrb bmrb);
		void attach_dssp();
		void attach_dssp(string);
		char getss(int);
		struct noeatoms query(int resid,string name);
		void printblock(vector<string>);
		void output(string filename);
		void caoutput(string filename);
		vector<int> getselect(string);
		vector<int> getselect(string, vector<int>);
        vector<int> getselect(vector<int> &,string);
		void clear_cs();
		void attach_coil(int);
		void attach_mean();
		void printpdb(char *filename, int flag=0);
		bool bisunknownmissing();
		void process_ambig(int flag);
		void attach_rmsf(vector<double>);
		void attach_rmsf(string);
		void print_rmsf(string);

		inline int getnatom(void) {return natom;};
		inline int getnatom2(void) {return natom2;};
		inline int getnres(void){return v.size();};
		inline vector<int> getheavy() {heavycoor();return heavy;};
		inline vector<int> getboundary() {heavycoor();return boundary;};
		inline void setcoor(vector<double> xx,vector<double> yy,vector<double> zz) {x=xx;y=yy;z=zz;};
		inline void setb(vector<double> bb){b=bb;};
		inline void setb2(vector<double> bb){b2=bb;};

		inline void setatomname(vector<string> name) {atomname=name;};
		inline string getseq(void) {return pdbseq;};
		inline void printpdb(string filename){printpdb(filename.c_str());};
		inline int getnmiss(void) {return nmiss;};

		inline vector<double> get_wishart(int id) {return v.at(id-1)->get_wishart();};
		inline int getvsize(void) {return v.size(); };
		inline char *getvoneletter(void) {return v_oneletternames;};

		int *code_pos;
};

struct modify
{
	string from;
	string to;
	char chain;
	int res;
};


class CPdb2 
{
private:
protected:
public:
	vector<class CPdb *> pdbs;
	CPdb2(void);
	~CPdb2(void);
		
	double loadpdb(string);

	inline vector<CPdb *> getpdbs(void) {return pdbs;};
	void clear(void);
};


#endif
