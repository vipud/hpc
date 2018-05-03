
#include "supply.h"
#include "bmrb.h"
#include  "aa.h"
#include "pdb.h"

#ifndef TRAJ
#define TRAJ

class CTraj 
{
private:

protected:
		double noedistance(vector<int> *, vector<int> *);
		double noedistance_frame(vector<int> *att1, vector<int> *att2,int n);

public:
		int nres;
		int natom,nframe;
		vector<double> x;
		vector<double> y;
		vector<double> z;
		double* x_arr;
		double* y_arr;
		double* z_arr;
		int x_size;
		int y_size;
		int z_size;

		vector<string> atomname;
		vector<double> rmsf;

		CTraj(void);
		~CTraj(void);

		inline void setnatom(int in){natom=in;};
		inline void setnres(int in){nres=in;};
		inline int getnframe(void) {return nframe;};
		virtual int loadcoor(string);
		int appendcoor(string);
		int set_range(int,int);
		virtual void clear(void);
		void dis_matrix(vector<int> *,vector<double> *);
		void getdihe(vector<dihe_group> *, vector<double> *);
		void getangle(vector<struct dihe_group> *, vector<double> *);

		void getring(vector<ring_group> *, vector<struct methyl_group> *, vector<double_five> *, enum methyl);
		void getring(vector<struct ring_group> *, vector<struct proton> *, vector<struct double_five> *);
		//void getring_acc(ring_group *, int, proton *, int, vector<struct double_five> *); // OpenACC version
		void getring_acc(ring_group *, int, proton *, int, double_five *, int); // OpenACC version
		void getring(vector<ring_group> *, vector<struct proton>* , vector< vector<struct double_five> > *);
		void getring(vector<ring_group> *, vector<struct nh_group>* , vector<struct double_five> *ring_effect);
		//void getring_acc(ring_group *, int, nh_group *, int, vector<struct double_five> *); // OpenACC version
		void getring_acc(ring_group *, int, nh_group *, int, double_five *, int); // OpenACC version
		void getring(vector<ring_group> *, vector<struct nh_group>* , vector< vector<struct double_five>  > *ring_effect);

		//void getani_acc(ani_group *, int, nh_group *, int, vector<struct double_four> *); // OpenACC version
		void getani_acc(ani_group *, int, nh_group *, int, double_four *, int); // OpenACC version
		void getani(vector<ani_group> *, vector<struct nh_group>* , vector< vector<double_four> > *);
		void getani(vector<ani_group> *, vector<struct methyl_group>* , vector<double_four> *, enum methyl);
		//void getani_acc(ani_group *, int, proton *, int, vector<double_four> *); // OpenACC version
		void getani_acc(ani_group *, int, proton *, int, double_four *, int); // OpenACC version
		void getani_acc(ani_group *, int, proton *, int, double_four *); // UNUSED
		void getani(vector<ani_group> *, vector<proton> *, vector<double_four> *);
		void getani(vector<ani_group> *, vector<struct proton>* , vector< vector<double_four> > *);
		void getani(vector<ani_group> *, vector<struct nh_group> *, vector<double_four> *);

		void getring_bb(vector<ring_group> *, vector<struct bb_group> *, vector<double_five> *,enum bb_carbon);

		void gethbond(vector<bbhbond_group> *bond,vector<ehbond> *effect);
		//void gethbond_acc(bbhbond_group *bond,int,vector<ehbond> *effect); // OpenACC version
		void gethbond_acc(bbhbond_group *bond,int,ehbond *, int); // OpenACC version
		void gethbond(vector<bbhbond_group> *hbond,vector<ehbond> *effect, double cutoff);
		void gethbond(vector<bbhbond_group> *bond,vector< vector<ehbond> > *effect);

		void gethbond2(vector<bbhbond_group> *hbond,vector< vector<ehbond> > *effect);
		void getschbond2(vector<struct proton> *protons, vector<struct bbhbond_group> *bb, vector< vector<ehbond> > *effect, vector< vector<eschbond> > *effect_sc);
		void evulatenmrcons(vector<struct noeline> *,double);
		void evaluatenmrcons_frame(vector<struct noeline> *nmrcons, double cutoff);
		void rmsd_matrix(vector< vector<double> > *,vector<int> *,int);
		void getvector(vector<struct index_three>,vector<double> *,vector<double> *,vector<double> *);
		void do_rmsf(void);

		void get_contact(vector<int> pos, int* used, int used_size, vector<float> * result);
		void get_contact(float rc,float shift, vector<int> pos, vector<int> used, vector<float> * result);
		//void get_all_contacts(vector<struct bb_group> *bb, vector<struct index_two> *index, int index_size, int *c2, int c2_size, float *results, int results_size);
		void get_all_contacts(bb_group *, int , index_two *, int, int *, int, float *, int);
		void get_all_contacts_double(vector<struct bb_group> *bb, vector<struct index_two> *index, int index_size, int *c2, int c2_size, double *results, int results_size);

		void getcoor(int,int,double *,double *,double *);
		void getcoor(vector<int>,int,vector<double> *,vector<double> *,vector<double> *);
		void getcoor(vector<int>,int,vector<float> *,vector<float> *,vector<float> *);
		void getcoor(vector<int>,vector<float> *xx,vector<float> *yy,vector<float> *zz);



};

//class to process crystal data
class smtry
{
	protected:
		double matrix[3][3];
		double motion[3];
	public:
		void getdata(double x[3][3],double y[3]);
		void trans(double x[3],double y[3]);
		void print();
};


class CTraj2: public CTraj
{
private:

protected:
	vector<smtry> crystal_sym_array;
	double translation[3][3];
	double scale[3][3];
	bool btrans;
	int ncry;
	bool bsmtry;
public:
	vector<double> b;  //b factors
	CTraj2(void);
	~CTraj2(void);
	void clear();
	void process_tran(void);
	int loadcoor(string);
	int select(vector<int>);
	int unitcell(void);
	void ninecells(void);
};

#endif
