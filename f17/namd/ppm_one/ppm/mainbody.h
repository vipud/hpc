#include "supply.h"
#include "bmrb.h"
#include  "aa.h"
#include "pdb.h"
#include "traj.h"

#ifndef MAINBODY
#define MAINBODY


class CDihe_process
{
	private:	
	protected:
		vector<dihe_group> *dihe_index;
		vector<double> *angle;
		vector<double> out;
		vector< vector<double> > out2;
		int *table;

		
		bool test_good(int id,int cut);

	public:
		int ndihe;
		int nframe;
		vector<int> num;
		vector<double> *dihe;
		static int hbs[18];
		bool test(int id, int t1, int t2);
		bool test_proton(int id,int type);
		void hb_expand(int);
        void hb_expand2(int);
		void ca(int);
		void ca2(int);

		bool ca_static_new(int id);
		void ca_ann(int id);
		void md_new(int id);
        void md_new_detail(int n,int id);
		void md_ann(int id,int n);
		void for_fit(int);
		void fit_angle(int id);
		void proton(int);
		void allproton(int);
        void allproton2(int);

		void process(int,int,int,int);
		void process2(int,int,int,int);
		void process_fit(int,int);
		void process_static_new(int id,int cut,int order,int order2);
		void process_md_sep(int n,int id,int cut, int order, int order2);

		void init(vector<int> innum,vector<double> *indihe);
		void init(vector<int> innum,vector<double> *indihe,vector<double> *inangle);
		void init(vector<int> innum,vector<double> *indihe,vector<dihe_group> *indihe_index);

		vector<int> pos(int in);
		vector<int> pos_angle(int in);

		vector<double> output(void);
		vector< vector<double> > output2(void);
		
		CDihe_process(void);
		~CDihe_process(void);
};

class CMainbody
{
private:
protected:
	int buffer[20];
	int ndihe;
	int nres;
	int natom;
	int nconf;
	int *sep_table;

	string pdb_name;
	
	bool bnew;
	class CBmrb bmrb;
	class CPdb *pdb;
	class CTraj *traj;
	class CDihe_process dihe_process;
	vector<struct methyl_group> select;
	vector<dihe_group> dihe_index;
	vector<int> dihe_num;
	vector<ring_group> ring_index;
	ring_group *ring_index_new;
	int ring_index_size;
	vector<ring_group> ring_index_internal;
	vector<ring_group> ring_index_external;
	vector<struct ani_group> anistropy;
	ani_group * anistropy_new;
	int anistropy_size;
	vector<double> dihe;
	vector<double> angle;
	vector<struct bb_group> bb;
	bb_group *bb_arr;
	int bb_size;
	vector<struct nh_group> bbnh;
	nh_group *bbnh_arr;
	int bbnh_size;
	vector<struct bbhbond_group> hbond;
	bbhbond_group *hbond_arr;
	int hbond_size;
	vector<struct proton> protons;
	vector<struct proton> allprotons;
	vector<struct proton> allprotons3;
	proton * allprotons3_new;
	int allprotons3_size;
	vector<int> heavy;

	double hill(double contact, double v1, double v2);

	void ha_protons_acc(struct proton *);
	void index_acc(struct index_two *, int);


public:

	//ppm old version
	static double c_c[5][366];
	static double c_h_add[5];

	//used for proton prediciton
	///static 
	static double c_all[251];
	static double c_sep[98][251];
	static int sep[19];


	//ppm_static for backbone
	static double static_c[6][1048];
	static double static_h[2][9];
	static double hill_para[6][2];


	//ppm_static_ann for backbone
	static double p_ann_ca[31120];
	static double p_ann_cb[31120];
	static double p_ann_co[31120];
	static double p_ann_n[31120];
	static double p_ann_h[33838];
	static double p_ann_ha[33838];







	//***************************************************//
	//functions
	int loadpdb(string);
	int loadpdb(string,string);
	int loadpdb(CPdb *,CTraj *);
	void load(string);
	int set_range(int,int);

	void clear(vector<struct bb_group> &bb);
	vector<bb_group> clear_acc(vector<bb_group> bb);
	void clear(vector<struct proton> &protons);
	vector<proton> clear_filter(vector<proton> protons);



	void predict_bb(void);
	void predict_bb2(void);

	void predict_bb_static_new(void);
	void predict_bb_static_ann(void);


	
	void predict_proton(void);
	void predict_proton2(void);

	void predict_proton_static_new(void);


	void cal_error();
	void compare(char *,vector< vector<double> > );
	void inline print_prediction(string name){pdb->print_prediction(name);};
	CMainbody();
	~CMainbody();
};


#endif
