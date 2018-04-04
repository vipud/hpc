#ifndef ANNFILE
#define ANNFILE

#include "lmmin.h"


typedef struct
{
	int n_dim;
	int n_neuron;
	double *input;  //size is n_dim*m_dat
	double *y;      //size is m_dat
	double (*f)(int ndim, int nneuron, double *x, const double *p, int offset);
}data_struct_neuron;



typedef struct
{
	int n_dim;
	int n_neuron;
	double *input;  //size is n_dim*m_dat
	double *y;      //size is m_dat
	double (*f)(int n_neuron, double *x, const double *p);
}data_struct_mix;


typedef struct
{
	int n_dim;
	int n_neuron;
	int n_conf;

	double *input; //size is n_dim*n_conf*m_dat
	double *y;     //size of m_dat
	double (*f)(int n_dim, int n_conf, int n_neuron, double *x, const double *p);
}data_struct_md;



class CAnn
{
private:
	int n_neuron;
	int n_dat;
	int n_par;
	int n_conf;

	double *p;
	double *x;
	double *y;
	double *v_min;
	double *v_max;
	double y_min;
	double y_max;

	double *p_save_flat;
	int p_save_size;

	void mapminmax(void);
	void mapminmax_md(void);
	void xapplyminmax(void);
	void xapplyminmax(double *xx);
	#pragma acc routine vector
	void xapplyminmax_acc(double *xx);
	void xapplyminmax_md(void);
	bool loadx(string name);
	bool loadx_md(string name);
	bool loadx(vector<vector<double> > xx);
	bool loady(string name);
	void initp(int n);
	void train_it(double percent);
	void train_it_md(double percent);
	void push_p(void);
	void randperm(int n,int *perm);
	void savehead(string filename, int n);
	void savep(string filename,int n);

protected:

public:
	int n_dim;
	int train(int,string,string,int,double);
	int train_md(int,string,string,int,double);
	vector<double> predict(int,string,string,vector<vector< double> >);
	vector<double> predict_md(int,string,string,vector<vector< double> >);
	double predict_one(double *xx, int vec_size);
	#pragma acc routine vector
	double predict_one_acc(double *xx, int vec_size);
	double predict_one_md(int,vector<double>);
	double assess(string,string,string);
	double assess_md(string,string,string);
	void load(string filename);
	void loadp(double *); // Now with data directives
	void save(string filename);
	void close(void);
	
	CAnn();
	~CAnn();

	//callback functions.
	static double myfunc_neuron(int ndim, int nneuron, double *x, const double *p, int offset);
	#pragma acc routine seq
	static double myfunc_neuron_acc(int ndim, int nneuron, double *x, const double *p, int offset);
	static void   evaluation_neuron(const double *par, int n_dat, const void *pdata, double *fvect, int *user);

	static double myfunc_mix(int nneuron, double *x, const double *p);
	static void   evaluation_mix(const double *par, int n_dat, const void *pdata, double *fvect, int *user);

	static double myfunc_md(int ndim,int nconf,int nneuron, double *x,const double *p);
	static void evaluation_md(const double *par, int n_dat, const void *pdata, double *fvect, int *user);
};


#endif
