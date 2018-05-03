//
#ifndef SUPPLY
#define SUPPLY

#define PI      3.1415926536

enum bb_carbon  {bb_ca,bb_cb,bb_co};
enum methyl {hydrogen,carbon};


struct dihe_group
{
	int id;
	char code;
	int type;
	bool bgood;


	int x1,x2,x3,x4;
};

struct index_three
{
	int x1,x2,x3;
};

struct index_two
{
	int x1,x2;
};

struct ring_group
{
	int x1;
	int x2,x3,x4,x5,x6,x7;

	int id;
	char code;
	bool bgood;
};

// Edited these data structures to start with initial values of 0.0
struct double_five
{
	double x[5] {0.0, 0.0, 0.0, 0.0, 0.0};
};

struct double_four
{
	double x[4] {0.0, 0.0, 0.0, 0.0};
};

struct double_three
{
	double x[3] {0.0, 0.0, 0.0};
};

struct methyl_group
{
	int id;
	char code;
	int type;
	int cpos;
	int hpos[6];
	double exp_h;
	double exp_c;
	double shifts_h;
	double shifts_c;
};


struct proton
{
	int id;
	char code;
	int type;
	string name;
	string cname;
	int hpos[6];
	int nh;
	int cpos;
	double exp;
	double exp_c;

	//used for HB2 HB3 ambiguity. 
	bool multy;
	int type2;
	string name2;
	string cname2;
	double exp1;
	double exp2;
};



struct bb_group
{
	int chain;  //chain index
	int id;     //residue index
	int id0;    //ID in PDB file
	char code;  //AA one letter name
	char ss;    //secondary structures

	int exploaded;           //experiment CS exist    
	bool previous_mut; 
	bool follow_mut;

	int hpos;
	int npos;
	int capos;
	int cbpos;
	int copos;
	int opos;
	int hapos;
	int hapos2;
 
	int follow_hpos;
	int follow_npos;
	double follow_exp_h;
	double follow_exp_n;
	double follow_exp_ca;

	double exp_ha;
	double exp_ha2;
	double exp_h;
	double exp_n;
	double exp_ca;
	double exp_cb;
	double exp_co;

	double pre_ca;
	double pre_cb;
	double pre_c;
	double pre_n;
	double pre_h;
	double pre_ha;

	vector<int> previous;
	vector<int> follows;
};

struct nh_group
{
	int id;
	char code;
	int hpos;
	int npos;

	double exp_h;
	double exp_n;

	double sparta_h;
	double sparta_n;
};

struct co_group
{
	int id;
	char code;

	int cpos;
	int opos;

	double exp_c;
	double exp_o;
};

struct ani_group
{
	int id;
	char code;
	int type;
	int pos[3];
};

struct bbhbond_group
{
	int id;
	char code;
	int npos;
	int hpos;
	int cpos;
	int opos;
	int type;
};

struct ehbond
{
	double n_length=0;
	double n_phi=0;
	double n_psi=0;

	double c_length=0;
	double c_phi=0;
	double c_psi=0;

};

struct eschbond
{
	double n_length;
	double n_phi;
	double n_psi;

	double c_length;
	double c_phi;
	double c_psi;

	int id;
	int type;
};

struct noeatoms
{
	vector< vector<int> > atoms;
	double length;
};

struct noeline
{
	int group;
	int id;
	int multi;
	int resid1,resid2;
	string atomname1,atomname2;
	string resname1,resname2; 
	int oldresid1,oldresid2;
	string oldatomname1,oldatomname2;
	double a,b,c;
	vector<double> obs;
	vector<int> pos1;
	vector<int> pos2;
	struct noeatoms index1,index2;
	double d;
	bool bvio;
};

struct diheline
{
	int id;
	vector<int>  resid;
	vector<string> resname;
	vector<string> atomname;
	vector< noeatoms > index;
	double upper;
	double lower;
	double middle;
	double delta;
};



class CCommandline
{
	private:
		int narg;
		vector<string> arguments;
		vector<string> parameters;
		vector<string> informations;
	protected:
	
	public:
		void pharse(int argc, char** argv);
		void init(vector<string>,vector<string>);
		void init(vector<string>,vector<string>,vector<string>);
		void print();
		string query(string);
		CCommandline();
		~CCommandline();
};




class CRmsd
	{
	private:
		float mov_com0;
		float mov_com1;
		float mov_com2;
		float ref_com0;
		float ref_com1;
		float ref_com2;
		float mov_to_ref0;
		float mov_to_ref1;
		float mov_to_ref2;
		float U00,U01,U02,U10,U11,U12,U20,U21,U22; 
		float R00;float R01;float R02;
		float R10;float R11;float R12;
		float R20;float R21;float R22;
		float vec00;float vec01;float vec02;
		float vec10;float vec11;float vec12;
		float vec20;float vec21;float vec22;
		float E0, residual;
		float RtR00,RtR01,RtR02,RtR10,RtR11,RtR12,RtR20,RtR21,RtR22;
		float left_eigenvec00,left_eigenvec01,left_eigenvec02,left_eigenvec10,left_eigenvec11,left_eigenvec12,left_eigenvec20,left_eigenvec21,left_eigenvec22;
		float right_eigenvec00,right_eigenvec01,right_eigenvec02,right_eigenvec10,right_eigenvec11,right_eigenvec12,right_eigenvec20,right_eigenvec21,right_eigenvec22;
		float eigenval0,eigenval1,eigenval2;


		void setup_rotation(float x[],float y[],float z[],
						float x0[],float y0[],float z0[], 
						int n_list);
		int jacobi3(int* n_rot);
		int diagonalize_symmetric();
		void get_rotation_matrix();
		int calculate_rotation_matrix();


	protected:
	public:
		float calculate_rotation_rmsd(float x[],float y[],float z[],float x0[],float y0[],float z0[],int n_list);
		float direct_rmsd(float x[],float y[],float z[],float x0[],float y0[],float z0[],int n_list);
		CRmsd();
		~CRmsd();
	};

namespace Sequence
{
	int code2pos(char code);
	string code2name(char c);
	void code2array(char code, int buffer[20]);
	void code2same(char code, int buffer[20]);
	char name2code(string in);
	vector<double> expand(char ,vector<double> *);
	vector<double> expand2(int ,vector<double> *);
	vector<int> align(string c1,string c2);
	vector<int> aligno(string c1,string c2, string &out1, string &out2,string &out3);
};

namespace ldw_math
{

	const double pi=3.14159265358979;
	double gaussrand(void);
	vector<int> cluster_pick2(int,vector<double>, vector<double>, double);
	double veclength(double x[3]);
	void cross(double z[3],double x[3],double y[3]);
	double dot(double x[3],double y[3]);

	double coor_to_angle(double x2,double y2,double z2,double x3,double y3,double z3,double x4,double y4,double z4);
	double coor_to_length(double x3,double y3,double z3,double x4,double y4,double z4);
	double coor_to_dihe(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double x4,double y4,double z4);

	void fit(vector<double> x,vector<double> y, vector<double> *z,double *a,double *b,double *rms,double *r);
	int dsvd(double a[6][3], int m, int n, double *w, double v[3][3]);
	int dsvd2(double *a, int m, int n, double *w, double v[3][3]);
	double area( double a, double b, double c );
	double PYTHAG(double a, double b);
	double mymax(double a, double b);
	double mysign(double a,double b);
	double effect(double x[6][3], int m, double ori[3], double p1[3]);
	void project(double ori[3], double p1[3], double p2[3]);
	void ring(double x[6][3], int m, double ori[3]);
	void regression_plane(double *x, int m, double ori[3]);


	void rotation_around_axis(double point[3], double ori[3], double theta);
}


#endif
