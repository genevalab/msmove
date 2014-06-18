#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define STATE1 '1'
#define STATE2 '0'
#define SITESINC 10 
#define SEGINC 80 

/* Data structures */
struct devent {
	double time;
	int popi;
	int popj;
	double pmove;
	double paramv;
	double **mat;
	char detype;
	struct devent *nextde;
};

struct c_params {
	int npop;
	int nsam;
	int *config;
	double **mig_mat;
	double r;
	int nsites;
	double f;
	double track_len;
	double *size;
	double *alphag;
	struct devent *deventlist;
};

struct m_params {
	double theta;
	int segsitesin;
	int treeflag;
	int timeflag;
	int mfreq;
};

struct params { 
	struct c_params cp;
	struct m_params mp;
};

struct node {
	int abv;
	int ndes;
	double time;
};

struct segl {
	int beg;
	struct node *ptree;
	int next;
};

/* Function prototypes */
void getpars(int, char*[], int*);
int gensam(char**, double*, double*, double*, int*);
char **cmatrix(int, int);
struct segl *segtre_mig(struct c_params*, int*, int*);
double ttime(struct node*, int);
double ttimemf(struct node*, int, int);
void prtree(struct node*, int);
int make_gametes(int, int,  struct node*, double, int, int, char**);
void ndes_setup(struct node*, int);
int biggerlist(int, char**);
int locate(int, double, double, double*);
void addtoelist(struct devent*, struct devent*); 
void argcheck(int, int, char**);
void free_eventlist(struct devent*, int);
int usage(void);
int pickb(int, struct node*, double);
int pickbmf(int, int, struct node*, double);
void parens(struct node*, int*, int*, int);
int tdesn(struct node*, int, int);
int pick2(int, int*, int*);
int ordran(int, double[]);
int mnmial(int, int, double[], int[]);
int order(int, double[]);
int ranvec(int, double[]);
int poisso(double);
double gasdev(double, double);
int re(int);
int xover(int, int, int);
int cinr(int, int);
int cleftr(int);
int ca(int, int, int, int);
int isseg(int, int, int*);
int pick2_chrom(int, int[], int*, int*);
int links(int);
int dbinom(int, double);
void init_genrand(unsigned long);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void) ;