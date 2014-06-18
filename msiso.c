/* msiso.c
 * author: Daniel Garrigan
 * last updated: Jan. 30 2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BASE 10
#define BINOM(x) ((x) * ((x) - 1)) / 2
#define NANALYSES 6

/* Define data structures */
typedef struct {
	int segsites;     /* Number of segregating sites */
	double piTotal;   /* Average pairwise differences (pi) for entire sample */
	double *piW;	  /* Pi within each population */
	double *piVarW;   /* Variance of pairwise differences within populations */
	double *piWmax;   /* Maximum pairwise difference within populations */
	double *piB;	  /* Pi between each pair of populations */
	double *piVarB;   /* Variance in pairwise diffs between population pairs */
	double *piBmin;   /* Minimum pairwise differences in population pairs */
	double psi;       /* Wakeley's 1996 ad-hoc psi statistic */
} stats_t;

typedef int (*analysis_f)(stats_t*, const int*, const int, char**);

/* Function prototypes */
int print_header(const int, FILE*);
int print_result(stats_t, const int, FILE*, int);
int pi_total(stats_t*, const int *, const int, char**);
int pi_within(stats_t*, const int*, const int, char**);
int pi_varw(stats_t*, const int*, const int, char**);
int pi_between(stats_t*, const int*, const int, char**);
int pi_varb(stats_t*, const int*, const int, char**);
int calc_psi(stats_t*, const int*, const int, char **);
double frequency(const char, const int, const int, const int, char**);
char **cmatrix(const int);
int biggerlist(const unsigned int, char**);
int free_stats(stats_t);

/* Globally scoped variables */
int maxsites = 1000;
int nsam;


int main(int argc, char *argv[])
{
	int i, howmany, count, npops, *config;
	int migflag = 0;
	double *posit;
	char **list, line[1001], dum[20], *mscanline;
	stats_t stats;
	FILE *fin, *fout;
	analysis_f af[NANALYSES] = {&pi_total, &pi_within, &pi_varw, &pi_between,
						        &pi_varb, &calc_psi};

	/* Set input stream */
	fin = stdin;

	/* Read and parse the ms command line */
	fgets(line, sizeof(line), fin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);

	/* Scan for sample size configuration */
	mscanline = strstr(line, "-I");

	/* If present parse into config array, else config=nsam */
	if (mscanline)
	{
		mscanline = strchr(mscanline, ' ');
		npops = (int)strtol(mscanline, &mscanline, BASE);
		if ( !(config = (int *)malloc((size_t)npops * sizeof(int))) )
			perror("malloc error in main(). config");
		for (i = 0; i < npops; i++)
			config[i] = (int)strtol(mscanline, &mscanline, BASE);
	}
	else
	{
		npops = 1;
		if ( !(config = (int *)malloc((size_t)npops * sizeof(int))) )
			perror("alloc error in main(). config");
		config[0] = nsam;
	}

	/* Gobble line with rng seeds */
	fgets(line, sizeof(line), fin);

	/* Allocate memory for chromosomes and site positions */
	list = cmatrix(maxsites + 1);
	if ( !(posit = (double *)malloc((size_t)maxsites * sizeof(double))) )
		perror("malloc error in main(). posit");

	/* Initialize replicate counter */
	count = 0;

	/* Set output stream */
	fout = stdout;
	
	/* Print header for output */
	print_header(npops, fout);

	/* Iterate through repliates */
	while (howmany-count++)
	{

		/* Read lines until data */
		while (line[0] != '/')
			fgets(line, sizeof(line), fin);

		if (line[2] == '*')
			migflag = 1;

		/* Get segregating sites */
		fscanf(fin, "  segsites: %d", &stats.segsites);

		/* Reallocate if too many sites for current matrix */
		if (stats.segsites >= maxsites)
		{
			maxsites = stats.segsites + 10;
			posit = (double*)realloc(posit, (size_t)maxsites * sizeof(double));
			biggerlist(maxsites, list);
		}

		/* Gobble positions line */
	    fgets(line, sizeof(line), fin);

	    /* Read in positions and chromosome matrix */
	   	if (stats.segsites > 0)
		{
			fscanf(fin, " positions:");
			for (i = 0; i < stats.segsites; i++)
				fscanf(fin, " %lf", posit + i);
			for (i = 0; i < nsam; i++)
				fscanf(fin, " %s", list[i]);
		}

		/* Calculate statistics */
		for (i = 0; i < NANALYSES; i++)
			af[i](&stats, config, npops, list);

		/* Print results to stdout */
		print_result(stats, npops, fout, migflag);

		// reset migflag
		migflag = 0;

		/* Take out the within-loop garbage */
		free_stats(stats);
	}

	free(posit);
	free(config);

	return 0;
}


/* Printer header for output */
int print_header(const int npops, FILE *fout)
{
	int i;
	
	for (i = 0; i < npops; i++)
		fprintf(fout, "Pi[%d]    \tPiWsd[%d] \tPiWmax[%d]\t", i+1, i+1, i+1);
	for (i = 0; i < npops * (npops - 1) / 2; i++)
		fprintf(fout, "PiB      \tPiBmin   \tpiBsd    \t");
	fprintf(fout, "Psi\t    mig\n");
	
	return 0;
}

/* Print statistics to stream fout */
int print_result(stats_t stats, const int npops, FILE *fout, int migflag)
{
	int i;
	for (i = 0; i < npops; i++)
		fprintf(fout, "%-9.5lf\t%-9.5lf \t%-9.5lf  \t", stats.piW[i], sqrt(stats.piVarW[i]), stats.piWmax[i]);
	for (i = 0; i < npops * (npops - 1) / 2; i++)
		fprintf(fout, "%-9.5lf\t%-9.5lf\t%-9.5lf\t", stats.piB[i], stats.piBmin[i], sqrt(stats.piVarB[i]));
	fprintf(fout, "%-9.5lf\t", stats.psi);
	if (migflag)
		fprintf(fout, "1\n");
	else
		fprintf(fout, "0\n");
	return 0;
}


/* Calculate pi for the entire sample */
int pi_total(stats_t *stats, const int *config, const int npops, char **list)
{
	int i;
	double fq;

	stats->piTotal = 0;

	for (i = 0; i < stats->segsites; i++)
	{
		fq = frequency('1', i, 0, nsam, list);
		stats->piTotal += 2.0 * fq * (nsam - fq) / (nsam * (nsam - 1.0));
	}

	return 0;
}


/* Calculate average within population pi */
int pi_within(stats_t *stats, const int *config, const int npops, char **list)
{
	int i, j, first, last;
	double fq, n;

	/* Allocate memory for piW */
	if ( !(stats->piW = (double *)calloc((size_t)npops, sizeof(double))) )
		perror("calloc error in pi_within()");

	for (i = 0; i < stats->segsites; i++) {
		first = 0;
		last = 0;
		for (j = 0; j < npops; j++)
		{
			n = (double)config[j];
			last += config[j];
			fq = frequency('1', i, first, last, list);
			if (n > 1)
				stats->piW[j] += 2.0 * fq * (n - fq) / (n * (n - 1.0));
			first += config[j];
		}
	}

	return 0;
}


/* Calculate variance of pi within populations */
int pi_varw(stats_t *stats, const int *config, const int npops, char **list)
{
	int i, j, m, n;
	int first = 0, last = 0;
	int dx;

	/* Allocate memory for piVar */
	if ( !(stats->piVarW = (double *)calloc((size_t)npops, sizeof(double))) )
		perror("calloc error in pi_varw().piVarW");
	if ( !(stats->piWmax = (double *)calloc((size_t)npops, sizeof(double))) )
		perror("calloc error in pi_varw().piWmax");

	for (j = 0; j < npops; j++)
	{
		last += config[j];
		for (m = first; m < last - 1; m++)
		{
			for (n = m + 1; n < last; n++)
			{
				dx = 0;
				for (i = 0; i < stats->segsites; i++)
					dx += list[m][i] == list[n][i] ? 0 : 1;
				stats->piWmax[j] = dx > stats->piWmax[j] ? dx : stats->piWmax[j];
				stats->piVarW[j] += (stats->piW[j] - dx) * (stats->piW[j] - dx);
			}
		}
		stats->piVarW[j] *= 2.0 / (config[j] * config[j]);
		first += config[j];
	}

	return 0;
}


/* Calculate average between population pi */
int pi_between(stats_t *stats, const int *config, const int npops, char **list)
{
	int i, j, k, m;
	int first1, last1, first2, last2;
	int bc = BINOM(npops);
	double f1, f2, n1, n2;

	/* Allocate memory for piB */
	if ( !(stats->piB = (double *)calloc((size_t)bc, sizeof(double))) )
		perror("calloc error in pi_between()");

	for (i = 0; i < stats->segsites; i++)
	{
		m = 0;
		first1 = 0;
		last1 = 0;
		for (j = 0; j < npops - 1; j++)
		{
			n1 = config[j];
			last1 += config[j];
			first2 = last1;
			for (k = j + 1; k < npops; k++)
			{
				last2 = first2 + config[k];
				n2 = config[k];
				f1 = frequency('1', i, first1, last1, list);
				f2 = frequency('1', i, first2, last2, list);
				stats->piB[m] += (f1 * (n2 - f2) + f2 * (n1 - f1)) / (n1 * n2);
				m++;
				first2 = last2;
			}
			first1 += config[j];
		}
	}

	return 0;
}


/* Calculate variance of pi between populations */
int pi_varb(stats_t *stats, const int *config, const int npops, char **list)
{
	int i, j, k, m, n, p = 0;
	int first1 = 0, last1 = 0, first2, last2;
	int bc = BINOM(npops);
	int dxy;
	double n1, n2;

	/* Allocate memory for piVar and piBmin */
	if ( !(stats->piVarB = (double *)calloc((size_t)bc, sizeof(double))) )
		perror("calloc error in pi_varb().piVarB");
	if( !(stats->piBmin = (double*)malloc((size_t)bc * sizeof(double))) )
		perror("malloc error in pi_varb().piBmin");

	/* Initialize min_piB */
	for (i = 0; i < bc; i++)
		stats->piBmin[i] = 1e8;

	for (j = 0; j < npops - 1; j++)
	{
		n1 = config[j];
		last1 += config[j];
		first2 = last1;
		for (k = j + 1; k < npops; k++)
		{
			last2 = first2 + config[k];
			n2 = config[k];
			for (m = first1; m < last1; m++)
			{
				for (n = first2; n < last2; n++)
				{
					dxy = 0;
					for (i = 0; i < stats->segsites; i++)
						dxy += list[m][i] == list[n][i] ? 0 : 1;
					stats->piBmin[p] = dxy < stats->piBmin[p] ? dxy : stats->piBmin[p];
					stats->piVarB[p] += (dxy - stats->piB[p]) * (dxy - stats->piB[p]);
				}
			}
			stats->piVarB[p] /= n1 * n2;
			++p;
			first2 = last2;
		}
		first1 += config[j];
	}

	return 0;
}


/* Calculate the psi statistics of Wakeley (1996) */
/* NOTE: needs stats->piTotal */
int calc_psi(stats_t *stats, const int *config, int npops, char **list)
{
	stats->psi = 1.0 / (nsam * (nsam - 1));
	stats->psi *= config[0] * (config[0] - 1) * 
				  (sqrt(stats->piVarW[0]) / stats->piW[0]) +
	       		  config[1] * (config[1] - 1) * 
	       		  (sqrt(stats->piVarW[1]) / stats->piW[1]) +
	       		  2 * config[0] * config[1] * 
	       		  (sqrt(stats->piVarB[0]) / stats->piTotal);

	return 0;
}


/* Calculate allele frequency at a given site*/
double frequency(const char allele, const int site, const int first, 
			     const int last, char **list)
{
	int i, count = 0;

	for (i = first; i < last; i++)
		count += list[i][site] == allele ? 1 : 0;

	return (double)count;
}


/* Allocates space for chromosomes (character strings) */
char **cmatrix(const int nsites)
{
	int i;
	char **m;

	if ( !(m = (char**)malloc((size_t)nsam * sizeof(char*))) )
		perror("malloc error in cmatrix()");

	for (i = 0; i < nsam; i++)
	{
		if( !(m[i] = (char*)malloc((size_t)nsites * sizeof(char))) )
		 	perror("malloc error in cmatrix(). 2");
	}

	return m;
}


/* Resize the memory block for chromosome alignment */
int biggerlist(const unsigned int nmax, char **list)
{
	int i;

	maxsites = nmax;

	for (i = 0; i < nsam; i++)
	{
		list[i] = (char*)realloc(list[i], (size_t)maxsites * sizeof(char));
		if (list[i] == NULL)
			perror("realloc error in biggerlist()");
	}

	return 0;
}


/* Free the statistics vectors */
int free_stats(stats_t stats)
{
	free(stats.piW);
	free(stats.piB);
	free(stats.piVarB);
	free(stats.piVarW);
	free(stats.piBmin);
	free(stats.piWmax);

	return 0;
}

