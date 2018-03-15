// ///////////////////////////////////////////////////////////////////////////

#include "bmrb.h"
#include "supply.h"


#ifndef AMINOACID
#define AMINOACID


struct S2
{
	string name1;
	string name2;
	double exp;
	double pre;
};

struct ired
{
	int pos;
	int id;
	int id0;
	int chain;
	int index1;
	int index2;

	char code;

	struct S2 s2;
};



struct Atom
{
	int index;			// index in PDB
	string name;		//name in PDB
	string cs_name;		// same cs_name if same chemical_shifts, used for methyl protons
	string base_name;  //carbon position name, like HB2 and HB3 have same base_name HB
	string carbon_name; //heavy atom name that this proton attached, might also be N or O, beside Carbon
	double cs_exp;
	double cs_pre;
	double cs_coil;
	double cs_mean;
	double cs_wang;
	double x,y,z;
	double rmsf;
	bool bb;
	bool proton;
	int proton_type;
	int ambig;         // ambigurious code when loading in BMRB file.
};



class CAminoacid 
{
private:	
protected:

	int  residue;   //residue index
	int  original_residue;
	int  previousc,followingn;
	char ss;  //secondary structure code
	//vector<struct Atom> atoms;
	vector<struct S2> order_parameters; 
	double pre_ca,pre_cb,pre_c,pre_h,pre_n,pre_ha;

	//random coil CS, corrections. etc
	static double wishart[22][12];
	static double wang_rc[21][6];
	static double wang_c1[21][6];
	static double wang_c2[21][6];

	int wishart_index();
	int wang_correct_index(char);
	int wang_rc_index(char);
	void set_wang(double *,int);
	

public:
	vector<struct Atom> atoms;
	Atom *atoms_arr;
	int atoms_size;
	struct Atom atom_nouse;
	bool bexploaded;
	char OneLetterName;
	char ThreeLetterName[4];
	int exploaded;
	int  chain;     //chain index
	double coil_pre[6];
	double coil_fol[6];

	CAminoacid(void);
	~CAminoacid(void);
	struct Atom get(const char* name);
	struct Atom* get_address(const char* name);
	void bbani(vector<struct ani_group> *);
	void bbhbond(vector<bbhbond_group> *);
	void bbnh(vector<struct nh_group> *);
	void bbco(vector<struct co_group> *);
	void caha(vector<struct index_three> *);
	void ired(vector<struct ired> *,int);
	void clearred(void);
	void loadred(struct ired);
	void bbdihe(vector<dihe_group> *);
	void sccoor(vector<int> *t);
	void bb(vector<bb_group> *);
	void follow_bb(vector<bb_group> *);
    void follow_bb_assign(vector<bb_group> *);
	void previous_bb(vector<bb_group> *);
	void loadexp(struct CBmrbdata data);
	void clearexp();
	void bbheavycoor(vector<int> *);
	void heavycoor(vector<int> *);
	void allcoor(vector<int> *t);
	int  get_proton(struct proton *t);
	int  get_proton3(struct proton *t);
	void setnterminal(void);
	void setcterminal(void);
	void proton3(vector<struct proton> *);
	void set_coil(int);
	void set_coil_wc(char pre, char fol);
	void set_mean(void);
	void set_mismatch(void);
	void remove_ambig(int);
	void combine_hsamec(int);

	void attach_rmsf(vector<double>);
	void print_rmsf(FILE *);
	void attach_bbprediction(double *);
	void attach_bbprediction(double pre_ca, double pre_cb, double pre_c, double pre_n, double pre_h, double pre_ha);
	void attach_protonprediction(string,double);
	void print_prediction(int *,FILE *);
	void print_bbprediction(FILE *);
	void print_protonprediction(FILE *);
	void printpdb(FILE *fp,vector<string> atomname,vector<double> x,vector<double> y,vector<double> z,vector<double> b,int &n);

	//CS
	vector<double> get_wishart();


	//virtual functions
	virtual void process(vector<string>); 
	virtual void ani(vector<struct ani_group> *);
	virtual void dihe(vector<dihe_group> * dihe_index);
	virtual void proton2(vector<struct proton> *);
	virtual void ring(vector<ring_group> *);
	virtual struct noeatoms query(string name);
	virtual void schbond(vector<bbhbond_group> *);
	virtual void methyl_ambig(int);
	virtual void proton(vector<struct proton> *);

	//inline functions
	inline int getca(void) {return get("CA").index;}
	inline void setfollowingn(int n) {followingn=n;}
	inline void setpreviousc(int n) {previousc=n;}
	inline void setresidue(int n) {residue=n;};
	inline char name(void) {return OneLetterName;}
	inline int base(void) {return get("N").index-1;}
	inline void setdssp(char c) {ss=c;};
	inline char getdssp(void) {return ss;};

	
	
};

// ///////////////////////////////////////////////////////////////////////////////////////
class CAla : public CAminoacid
{
private:
public:
	CAla(void);
	~CAla(void);
	void dihe(vector<dihe_group> * dihe_index);
	void proton2(vector<struct proton> *);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CArg : public CAminoacid
{
private:

public:
	CArg(void);
	~CArg(void);
	void ani(vector<struct ani_group> *anistropy);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);	
	void schbond(vector<bbhbond_group> *grp);
};

// ///////////////////////////////////////////////////////////////////////////////////////
class CAsn : public CAminoacid
{
private:

public:
	CAsn(void);
	~CAsn(void);
	void ani(vector<struct ani_group> *anistropy);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
	void schbond(vector<bbhbond_group> *grp);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CAsp : public CAminoacid
{
private:
public:
	CAsp(void);
	~CAsp(void);
	void ani(vector<struct ani_group> *anistropy);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
	void schbond(vector<bbhbond_group> *grp);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CCys : public CAminoacid
{
private:

public:
	CCys(void);
	~CCys(void);
	void dihe(vector<dihe_group> * dihe_index);
	void print_prediction(int *,FILE *);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CCyx : public CAminoacid
{
private:
public:
	CCyx(void);
	~CCyx(void);
	void dihe(vector<dihe_group> * dihe_index);
	void sccoor(vector<int> *t);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CGln : public CAminoacid
{
private:

public:
	CGln(void);
	~CGln(void);
	void ani(vector<struct ani_group> *anistropy);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
	void schbond(vector<bbhbond_group> *grp);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CGlu : public CAminoacid
{
private:
	
public:
	CGlu(void);
	~CGlu(void);
	void ani(vector<struct ani_group> *anistropy);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
	void schbond(vector<bbhbond_group> *grp);
};
// ///////////////////////////////////////////////////////////////////////////////////////



class CGly : public CAminoacid
{
private:

public:
	CGly(void);
	~CGly(void);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CHis : public CAminoacid
{
private:

public:
	CHis(void);
	~CHis(void);
	void dihe(vector<dihe_group> * dihe_index);
	void ring(vector<ring_group> *);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CIle : public CAminoacid
{
private:

public:
	CIle(void);
	~CIle(void);
	void dihe(vector<dihe_group> * dihe_index);
	void proton2(vector<struct proton> *sel);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CLeu : public CAminoacid
{
private:

public:
	CLeu(void);
	~CLeu(void);
	void dihe(vector<dihe_group> * dihe_index);
	void proton2(vector<struct proton> *sel);
	struct noeatoms  query(string name);
	void methyl_ambig(int);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CLys : public CAminoacid
{
private:

public:
	CLys(void);
	~CLys(void);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CMet : public CAminoacid
{
private:
	
public:
	CMet(void);
	~CMet(void);
	void dihe(vector<dihe_group> * dihe_index);
	void proton2(vector<struct proton> *sel);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CPhe : public CAminoacid
{
private:

public:
	CPhe(void);
	~CPhe(void);
	void dihe(vector<dihe_group> * dihe_index);
	void ring(vector<ring_group> *);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CPro : public CAminoacid
{
private:

public:
	CPro(void);
	~CPro(void);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CSer : public CAminoacid
{
private:
public:
	CSer(void);
	~CSer(void);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
	void schbond(vector<bbhbond_group> *grp);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CThr : public CAminoacid
{
private:

public:
	CThr(void);
	~CThr(void);
	void dihe(vector<dihe_group> * dihe_index);
	void proton2(vector<struct proton> *sel);
	struct noeatoms  query(string name);
	void schbond(vector<bbhbond_group> *grp);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CTrp : public CAminoacid
{
private:

public:
	CTrp(void);
	~CTrp(void);
	void dihe(vector<dihe_group> * dihe_index);
	void ring(vector<ring_group> *);
	struct noeatoms  query(string name);
	void schbond(vector<bbhbond_group> *grp);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CTyr : public CAminoacid
{
private:
	
public:
	CTyr(void);
	~CTyr(void);
	void dihe(vector<dihe_group> * dihe_index);
	void ring(vector<ring_group> *ring);
	struct noeatoms  query(string name);
	void schbond(vector<bbhbond_group> *grp);
};
// ///////////////////////////////////////////////////////////////////////////////////////
class CVal : public CAminoacid
{
private:

public:
	CVal(void);
	~CVal(void);
	void dihe(vector<dihe_group> * dihe_index);
	void proton2(vector<struct proton> *sel);
	struct noeatoms  query(string name);
	void methyl_ambig(int);
};

class CMiss : public CAminoacid
{
private:
	string head;
	string resname;
public:
	CMiss(void);
	~CMiss(void);
	void process(vector<string>);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
};

class CUnk : public CAminoacid
{
private:
	string head;
	string resname;
public:
	CUnk(void);
	~CUnk(void);
	void process(vector<string>);
	void dihe(vector<dihe_group> * dihe_index);
	struct noeatoms  query(string name);
};

class CLigand
{
private:
	vector<int>  atomindexs;
	vector<string>  atomnames;
	string resname;
	int id;
	int  residue;   //residue index
	int  original_residue;
	string head; //HETATM or ATOM
public:
	CLigand(void);
	~CLigand(void);
	void heavycoor(vector<int> *);
	void process(vector<string> block);
	inline void setresidue(int n) {residue=n;};
};

#endif
