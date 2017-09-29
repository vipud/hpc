#ifndef RDC
#define RDC

#include "pdb.h"
#include "traj.h"

struct bbs2_pre
{
	int id;
	int h;
	int o;

	int self_begin;
	int self_stop;

	float s2;
	float experiment;
};


struct exp_rdc
{
	string name1,name2;
	string resname1,resname2;
	int atomid1,atomid2;
	int resid1,resid2;
	double vrdc;
	double pre_rdc;
};


class CNmr
{
private:
protected:
	class CPdb *pdb;
	class CTraj *traj;
	int nres,nconf,natom;
	bool bnew;

public:
	float pythag(float, float);
	float *myvector(int, int);
	void free_vector(float*, int, int);
	float** matrix(int, int, int, int);
	void free_matrix(float**, int, int, int, int);
	void svdcmp(float**, int, int, float*, float**);
	void svbksb(float**, int, int, float*, float**, float*, float*);
	void tred2(float **, int , float *, float *);
	void tqli(float *, float *, int , float **);
	int loadpdb(string);
	int loadpdb(CPdb *,CTraj *);
	CNmr();
	~CNmr();
};


class CRdc : public CNmr
{
private:
protected:

	int adj;
	ofstream fout;
	vector< vector<struct exp_rdc> > exp;
	
public:
	void fillin(void);
	float backcal(float *,float**, int, int, float*, float *,float*);
	float error();
    float error(int);
	void load(int adj);
	int load(string);
	void actualload(vector<string> *);
	
	CRdc();
	~CRdc();
};


class CS2 : public CNmr
{
private:
protected:
	
	int flag;
	int block_number,block_size;
	int ns2;
	vector<double> xx,yy,zz;
	vector<double> pre;
	vector<struct index_three> nh;
	vector<struct ired> red;

	
public:
	void init(int);
	void doit();
	void print(string);
	void contact_model(string,string);
	void methyl_contact(string outfile);
	CS2();
	~CS2();
};


#endif
