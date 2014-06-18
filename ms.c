#include "ms.h"
#include <time.h>
#include <string.h>

unsigned int maxsites = SITESINC;
int count, ntbs;
double *posit;
double segfac;
struct params pars;

int main(int argc, char *argv[])
{
	int i, k, howmany, segsites;
	int migflag;
	unsigned long int seed;
	double probss, tmrca, ttot;
	char **list, **tbsparamstrs;
	FILE *pf;

	ntbs = 0;
	tbsparamstrs = (char**)malloc((size_t)argc*sizeof(char*));

	for (i=0; i < argc; i++)
		printf("%s ", argv[i]);
	printf("\n");
	for (i=0; i < argc; i++)
		tbsparamstrs[i] = (char*)malloc(30*sizeof(char));
	for (i=1; i < argc; i++)
		if (strcmp(argv[i], "tbs") == 0)
			argv[i] = tbsparamstrs[ntbs++];

	count = 0;

	if (ntbs > 0)
		for (k = 0; k < ntbs; k++)
			scanf(" %s", tbsparamstrs[k]);
	getpars(argc, argv, &howmany);

	/* Initialize MT19937 rng from system clock */
	seed = (unsigned long int)time(0);
	init_genrand(seed);

	pf = stdout;
	fprintf(pf, "%ld\n", seed);

	if (pars.mp.segsitesin ==  0) 
	{
	     list = cmatrix(pars.cp.nsam, maxsites+1);
	     posit = (double*)malloc((size_t)maxsites*sizeof(double));
		 if (posit == NULL)
			 perror("Cannot allocate memory for posit");
	}
	else 
	{
	     list = cmatrix(pars.cp.nsam, pars.mp.segsitesin+1);
	     posit = (double*)malloc((size_t)pars.mp.segsitesin*sizeof(double));
		 if (posit == NULL)
			 perror("Cannot allocate memory for posit");
	     if (pars.mp.theta > 0) 
		 {
		    segfac = 1.0;
		    for (i=pars.mp.segsitesin; i > 1; i--)
				segfac *= i;
		 }
	}

	while (howmany-count++)
	{
		if ((ntbs > 0) && (count > 1)) 
		{
			for (k=0; k < ntbs; k++)
				if (scanf(" %s", tbsparamstrs[k]) == EOF)
					exit(EXIT_FAILURE);
			getpars(argc, argv, &howmany);
		}

		// generate coalescent history
		segsites = gensam(list, &probss, &tmrca, &ttot, &migflag);
		fprintf(pf, "\n//");
		if (migflag == 1)
			fprintf(pf, "*");
		if (ntbs > 0)
		{
			for (k=0; k < ntbs; k++)
				printf("\t%s", tbsparamstrs[k]);
		}
		printf("\n");
  		if (pars.mp.timeflag)
			fprintf(pf, "time:\t%lf\t%lf\n", tmrca, ttot);
        if ((segsites > 0) || (pars.mp.theta > 0))
		{
			if ((pars.mp.segsitesin > 0) && (pars.mp.theta > 0)) 
				fprintf(pf, "prob: %g\n", probss);
			fprintf(pf, "segsites: %d\n", segsites);
			if (segsites > 0)
				fprintf(pf, "positions: ");
			for (i=0; i < segsites; i++)
				fprintf(pf, "%6.4lf ", posit[i]);
			fprintf(pf, "\n");
			if (segsites > 0) 
				for (i=0; i < pars.cp.nsam; i++)
					fprintf(pf, "%s\n", list[i]); 
		}
	}
	return 0;
}

int gensam(char **list, double *pprobss, double *ptmrca, double *pttot, int *migflag) 
{
	int nsegs, i, k, seg, ns, start, end, len, segsit;
	int *ss, segsitesin, nsites, nsam, mfreq;
	double nsinv, tseg, tt, theta, es, *pk;
	struct segl *seglst;

	nsites = pars.cp.nsites;
	nsinv = 1.0/nsites;
	*migflag = 0;
	seglst = segtre_mig(&(pars.cp), &nsegs, migflag);
	nsam = pars.cp.nsam;
	segsitesin = pars.mp.segsitesin;
	theta = pars.mp.theta;
	mfreq = pars.mp.mfreq;

	if (pars.mp.treeflag) 
	{
	  	ns = 0;
	    for (seg=0, k=0; k < nsegs; seg = seglst[seg].next, k++)
		{
			if ( (pars.cp.r > 0) || (pars.cp.f > 0) )
			{
				end = k < nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1;
				start = seglst[seg].beg;
				len = end-start+1;
				fprintf(stdout, "[%d]", len);
			}
			prtree(seglst[seg].ptree, nsam);
			if ((segsitesin == 0) && (theta == 0) && (pars.mp.timeflag == 0)) 
				free(seglst[seg].ptree);
		}
	}

	if (pars.mp.timeflag)
	{
		tt = 0;
		for (seg=0, k=0; k < nsegs; seg = seglst[seg].next, k++)
		{ 
			if (mfreq > 1) 
				ndes_setup(seglst[seg].ptree, nsam);
			end = k < nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1;
			start = seglst[seg].beg;
			if ((nsegs == 1) || ((start <= nsites/2) && (end >= nsites/2)))
				*ptmrca = (seglst[seg].ptree+2*nsam-2)->time;
			len = end-start+1;
			tseg = len / (double)nsites;
			if (mfreq == 1)
				tt += ttime(seglst[seg].ptree, nsam) * tseg;
			else 
				tt += ttimemf(seglst[seg].ptree, nsam, mfreq)*tseg;
			if ((segsitesin == 0) && (theta == 0)) 
				free(seglst[seg].ptree);
		}
		*pttot = tt;
	}	
	
    if ((segsitesin == 0) && (theta > 0))
	{
		ns = 0;
		for (seg=0, k=0; k < nsegs; seg = seglst[seg].next, k++)
		{ 
			if (mfreq > 1)
				ndes_setup(seglst[seg].ptree, nsam);
			end = k < nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1;
			start = seglst[seg].beg;
			len = end-start+1;
			tseg = len*theta/nsites;
			if (mfreq == 1)
				tt = ttime(seglst[seg].ptree, nsam);
			else
				tt = ttimemf(seglst[seg].ptree, nsam, mfreq);
			segsit = poisso(tseg*tt);
			if ((unsigned int)(segsit+ns) >= maxsites)
			{
				maxsites = segsit+ns+SITESINC;
				posit = (double*)realloc(posit, (size_t)maxsites*sizeof(double));
				if (posit == NULL)
					perror("Realloc failure for posit");
				biggerlist(nsam, list); 
			}
			make_gametes(nsam, mfreq, seglst[seg].ptree, tt, segsit, ns, list);
			free(seglst[seg].ptree);
			locate(segsit, start * nsinv, len * nsinv, posit + ns);   
			ns += segsit;
		}
	}
	else
		if (segsitesin > 0)
		{
			pk = (double*)malloc((size_t)nsegs*sizeof(double));
			ss = (int*)malloc((size_t)nsegs*sizeof(int));
			if ((pk == 0) || (ss == 0))
				perror("malloc error. gensam.2");
			tt = 0;
			for (seg=0, k=0; k < nsegs; seg = seglst[seg].next, k++)
			{ 
				if (mfreq > 1)
					ndes_setup(seglst[seg].ptree, nsam);
				end = k < nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1;
				start = seglst[seg].beg;
				len = end-start+1;
				tseg = len/(double)nsites;
				if (mfreq == 1)
					pk[k] = ttime(seglst[seg].ptree, nsam)*tseg;
				else 
					pk[k] = ttimemf(seglst[seg].ptree, nsam, mfreq)*tseg;
				tt += pk[k];
			}
			if (theta > 0)
			{ 
				es = theta*tt;
				*pprobss = exp(-1.0*es)*pow(es, (double)segsitesin)/segfac;
			}
			if (tt > 0)
			{
				for (k=0; k < nsegs; k++)
					pk[k] *= 1.0/tt;
				mnmial(segsitesin, nsegs, pk, ss);
			}
			else
				for (k=0; k < nsegs; k++)
					ss[k] = 0;
			ns = 0;
			for (seg=0, k=0; k < nsegs; seg = seglst[seg].next, k++)
			{ 
				end = k < nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1;
				start = seglst[seg].beg;
				len = end-start+1;
				tseg = len/(double)nsites;
				make_gametes(nsam, mfreq, seglst[seg].ptree, tt*pk[k]/tseg, ss[k], ns, list);
				free(seglst[seg].ptree);
				locate(ss[k], start*nsinv, len*nsinv, posit+ns);
				ns += ss[k];
			}
			free(pk);
			free(ss);
		}
		
		for (i=0; i < nsam; i++)
			list[i][ns] = '\0';
		return ns;
}

void ndes_setup(struct node *ptree, int nsam)
{
	int i;
	for (i=0; i < nsam; i++)
		(ptree + i)->ndes = 1;
	for (i=nsam; i < 2*nsam-1; i++)
		(ptree+i)->ndes = 0;
	for (i=0; i < 2*nsam-2; i++)
		(ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes;
}

int biggerlist(int nsam, char **list)
{
	int i;
	for (i=0; i < nsam; i++)
	{
		list[i] = (char*)realloc(list[i], (size_t)maxsites*sizeof(char));
		if (list[i] == NULL)
			perror("realloc error. bigger");
	}
	return 0;
}
	   
char **cmatrix(int nsam, int len)
{
	int i;
	char **m;
	if (!(m = (char**)malloc((size_t)nsam*sizeof(char*))))
		perror("alloc error in cmatrix");
	for (i=0; i < nsam; i++)
		if (!(m[i] = (char*)malloc((size_t)len*sizeof(char))))
			perror("alloc error in cmatric. 2");
	return m;
}

int locate(int n, double beg, double len, double *ptr)
{
	int i;
	ordran(n, ptr);
	for (i=0; i < n; i++)
		ptr[i] = beg+ptr[i]*len;
	return 0;
}

void getpars(int argc, char *argv[], int *phowmany)
{
	int arg, i, j, sum, pop, argstart, npop, npop2, pop2;
	double migr, mij, psize, palpha;
	char ch3;
	struct devent *ptemp, *pt;
	FILE *pf;
	
	if (count == 0)
	{
		if (argc < 4)
		{
			fprintf(stderr, "Too few command line arguments\n");
			usage();
		}
		pars.cp.nsam = atoi(argv[1]);
		if (pars.cp.nsam <= 0)
		{
			fprintf(stderr, "First argument error. nsam <= 0. \n");
			usage();
		}
		*phowmany = atoi(argv[2]);
		if (*phowmany <= 0)
		{
			fprintf(stderr, "Second argument error. howmany <= 0. \n");
			usage();
		}
		pars.cp.r = pars.mp.theta = pars.cp.f = 0;
		pars.cp.track_len = 0;
		pars.cp.npop = npop = 1;
		pars.cp.mig_mat = (double**)malloc(sizeof(double*));
		pars.cp.mig_mat[0] = (double*)malloc(sizeof(double));
		pars.cp.mig_mat[0][0] = 0;
		pars.mp.segsitesin = 0;
		pars.mp.treeflag = 0;
 		pars.mp.timeflag = 0;
		pars.mp.mfreq = 1;
		pars.cp.config = (int*)malloc((size_t)(pars.cp.npop+1)*sizeof(int));
		(pars.cp.config)[0] = pars.cp.nsam;
		pars.cp.size= (double*)malloc((size_t)pars.cp.npop*sizeof(double));
		(pars.cp.size)[0] = 1.0;
		pars.cp.alphag = (double*)malloc((size_t)pars.cp.npop*sizeof(double));
		(pars.cp.alphag)[0] = 0;
		pars.cp.nsites = 2;
	}
	else
	{
		npop = pars.cp.npop;
		free_eventlist(pars.cp.deventlist, npop);
	}
  	pars.cp.deventlist = 0;
	arg = 3;

	while (arg < argc)
	{
		if (argv[arg][0] != '-')
		{ 
			fprintf(stderr, " argument should be -%s ?\n", argv[arg]); 
			usage();
		}
		switch (argv[arg][1])
		{
			case 'f' :
				if (ntbs > 0)
				{
					fprintf(stderr, " can't use tbs args and -f option.\n");
					exit(EXIT_FAILURE);
				}
				arg++;
				argcheck(arg, argc, argv);
				pf = fopen(argv[arg], "r");
				if (pf == 0)
				{
					fprintf(stderr, " no parameter file %s\n", argv[arg]);
					exit(EXIT_FAILURE);
				}
				arg++;
				argc++;
				argv = (char**)malloc((size_t)(argc+1)*sizeof(char*));
				argv[arg] = (char*)malloc(20*sizeof(char));
				argstart = arg;
				while (fscanf(pf, " %s", argv[arg]) != EOF)
				{
					arg++;
					argc++;
					argv = (char**)realloc(argv, (size_t)argc*sizeof(char*));
					argv[arg] = (char*)malloc(20*sizeof(char));
				}
				fclose(pf);
				argc--;
				arg = argstart;
				break;
			case 'r' : 
				arg++;
				argcheck(arg, argc, argv);
				pars.cp.r = atof(argv[arg++]);
				argcheck(arg, argc, argv);
				pars.cp.nsites = atoi(argv[arg++]);
				if (pars.cp.nsites < 2)
				{
					fprintf(stderr, "with -r option must specify both rec_rate and nsites > 1\n");
					usage();
				}
				break;		
			case 'c' : 
				arg++;
				argcheck(arg, argc, argv);
				pars.cp.f = atof(argv[arg++]);
				argcheck(arg, argc, argv);
				pars.cp.track_len = atof(argv[arg++]);
				if (pars.cp.track_len < 1.0)
				{
					fprintf(stderr, "with -c option must specify both f and track_len>0\n");
					usage();
				}
				break;		
			case 't' : 
				arg++;
				argcheck(arg, argc, argv);
				pars.mp.theta = atof(argv[arg++]);
				break;
			case 's' : 
				arg++;
				argcheck(arg, argc, argv);
				pars.mp.segsitesin = atoi(argv[arg++]);
				break;
			case 'F' : 
				arg++;
				argcheck(arg, argc, argv);
				pars.mp.mfreq = atoi(argv[arg++]);
				if ((pars.mp.mfreq < 2) || (pars.mp.mfreq > pars.cp.nsam/2))
				{
					fprintf(stderr, " mfreq must be >= 2 and <= nsam/2.\n");
					usage();
				}
				break;
			case 'T' : 
				pars.mp.treeflag = 1;
				arg++;
				break;
			case 'L' : 
				pars.mp.timeflag = 1;
				arg++;
				break;
			case 'I' : 
			    arg++;
			    if (count == 0)
				{
					argcheck(arg, argc, argv);
					pars.cp.npop = atoi(argv[arg]);
					pars.cp.config = (int*)realloc(pars.cp.config, (size_t)pars.cp.npop*sizeof(int));
					npop = pars.cp.npop;
				}
			    arg++;
			    for (i=0; i < pars.cp.npop; i++)
				{
					argcheck(arg, argc, argv);
					pars.cp.config[i] = atoi(argv[arg++]);
				}
			    if (count == 0)
				{
					pars.cp.mig_mat = (double**)realloc(pars.cp.mig_mat, (size_t)pars.cp.npop*sizeof(double*));
					pars.cp.mig_mat[0] = (double*)realloc(pars.cp.mig_mat[0], (size_t)pars.cp.npop*sizeof(double));
					for(i=1; i < pars.cp.npop; i++)
						pars.cp.mig_mat[i] = (double*)malloc((size_t)pars.cp.npop*sizeof(double));
					pars.cp.size = (double*)realloc(pars.cp.size, (size_t)pars.cp.npop*sizeof(double));
					pars.cp.alphag = (double*)realloc(pars.cp.alphag, (size_t)pars.cp.npop*sizeof(double));
			        for (i=1; i < pars.cp.npop; i++) {
						(pars.cp.size)[i] = (pars.cp.size)[0];
						(pars.cp.alphag)[i] = (pars.cp.alphag)[0];
					}
				}
				if ((arg < argc) && (argv[arg][0] != '-'))
				{
					argcheck(arg, argc, argv);
					migr = atof(argv[arg++]);
				}
			    else
					migr = 0;
				for (i=0; i < pars.cp.npop; i++) 
					for (j=0; j < pars.cp.npop; j++)
						pars.cp.mig_mat[i][j] = migr/(pars.cp.npop-1);
				for (i=0; i < pars.cp.npop; i++)
					pars.cp.mig_mat[i][i] = migr;
				break;
			case 'm' :
				if (npop < 2)
				{
					fprintf(stderr, "Must use -I option first.\n");
					usage();
				}
				if (argv[arg][2] == 'a')
				{
					arg++;
				    for (pop=0; pop < npop; pop++)
						for (pop2=0; pop2 < npop; pop2++)
						{
							argcheck(arg, argc, argv);
							pars.cp.mig_mat[pop][pop2] = atof(argv[arg++]);
						}
					for (pop=0; pop < npop; pop++)
					{
						pars.cp.mig_mat[pop][pop] = 0;
						for (pop2=0; pop2 < npop; pop2++)
							if (pop2 != pop)
								pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2];
					}	
				}
			    else {
					arg++;
					argcheck(arg, argc, argv);
					i = atoi(argv[arg++])-1;
					argcheck(arg, argc, argv);
					j = atoi(argv[arg++])-1;
					argcheck(arg, argc, argv);
					mij = atof(argv[arg++]);
					pars.cp.mig_mat[i][i] += mij -  pars.cp.mig_mat[i][j];
					pars.cp.mig_mat[i][j] = mij;
				}
				break;
			case 'n' :
				if (npop < 2)
				{
					fprintf(stderr, "Must use -I option first.\n");
					usage();
				}
			    arg++;
			    argcheck(arg, argc, argv);
			    pop = atoi(argv[arg++])-1;
			    argcheck(arg, argc, argv);
			    psize = atof(argv[arg++]);
			    pars.cp.size[pop] = psize;
			   break;
			case 'g' :
				if (npop < 2)
				{
					fprintf(stderr, "Must use -I option first.\n");
					usage();
				}
			    arg++;
			    argcheck(arg, argc, argv);
			    pop = atoi(argv[arg++])-1;
			    if (arg >= argc)
				{
					fprintf(stderr, "Not enough arg's after -G.\n");
					usage();
				}
			    palpha = atof(argv[arg++]);
			    pars.cp.alphag[pop] = palpha;
			   break;
			case 'G' :
			    arg++;
			    if (arg >= argc)
				{
					fprintf(stderr, "Not enough arg's after -G.\n");
					usage();
				}
			    palpha = atof(argv[arg++]);
			    for (i=0; i < pars.cp.npop; i++) 
					pars.cp.alphag[i] = palpha;
			   break;
			case 'e' :
			    pt = (struct devent*)malloc(sizeof(struct devent));
			    pt->detype = argv[arg][2];
			    ch3 = argv[arg][3];
			    arg++;
			    argcheck(arg, argc, argv);
			    pt->time = atof(argv[arg++]);
			    pt->nextde = 0;
			    if (pars.cp.deventlist == 0) 
					pars.cp.deventlist = pt;
			    else if (pt->time < pars.cp.deventlist->time)
				{ 
					ptemp = pars.cp.deventlist;
					pars.cp.deventlist = pt;
					pt->nextde = ptemp;
				}	
			    else
					addtoelist(pt, pars.cp.deventlist);
			    switch (pt->detype)
				{
					case 'N' :
						argcheck(arg, argc, argv);
						pt->paramv = atof(argv[arg++]);
						break;
					case 'G' :
						if (arg >= argc)
						{
							fprintf(stderr, "Not enough arg's after -eG.\n");
							usage();
						}
						pt->paramv = atof(argv[arg++]);
						break;
					case 'M' :
						argcheck(arg, argc, argv);
						pt->paramv = atof(argv[arg++]);
						break;
					case 'n' :
						argcheck(arg, argc, argv);
						pt->popi = atoi(argv[arg++])-1;
						argcheck(arg, argc, argv);
						pt->paramv = atof(argv[arg++]);
						break;
					case 'g' :
						argcheck(arg, argc, argv);
						pt->popi = atoi(argv[arg++])-1;
						if (arg >= argc)
						{
							fprintf(stderr, "Not enough arg's after -eg.\n");
							usage();
						}
						pt->paramv = atof(argv[arg++]);
						break;
					case 's' :
						argcheck(arg, argc, argv);
						pt->popi = atoi(argv[arg++])-1;
						argcheck(arg, argc, argv);
						pt->paramv = atof(argv[arg++]);
						break;
					case 'm' :
						if (ch3 == 'a')
						{
							pt->detype = 'a';
							argcheck(arg, argc, argv);
							npop2 = atoi(argv[arg++]);
							pt->mat = (double**)malloc((size_t)npop2*sizeof(double*));
							for (pop=0; pop < npop2; pop++)
							{
								(pt->mat)[pop] = (double*)malloc((size_t)npop2*sizeof(double));
								for (i=0; i < npop2; i++)
								{
									if (i == pop)
										arg++;
									else
									{
										argcheck(arg, argc, argv);
										(pt->mat)[pop][i] = atof(argv[arg++]);
									}
								}
							}
							for (pop=0; pop < npop2; pop++)
							{
								(pt->mat)[pop][pop] = 0;
								for (pop2=0; pop2 < npop2; pop2++)
									if (pop2 != pop)
										(pt->mat)[pop][pop] += (pt->mat)[pop][pop2];
							}
						}
						else
						{
							argcheck(arg, argc, argv);
							pt->popi = atoi(argv[arg++])-1;
							argcheck(arg, argc, argv);
							pt->popj = atoi(argv[arg++])-1;
							argcheck(arg, argc, argv);
							pt->paramv = atof(argv[arg++]);
						}
						break;
					case 'j' :
						argcheck(arg, argc, argv);
						pt->popi = atoi(argv[arg++])-1;
						argcheck(arg, argc, argv);
						pt->popj = atoi(argv[arg++])-1;
						break;
					case 'v' :
						argcheck(arg, argc, argv);
						pt->popi = atoi(argv[arg++])-1;
						argcheck(arg, argc, argv);
						pt->popj = atoi(argv[arg++])-1;
						argcheck(arg, argc, argv);
						pt->pmove = atof(argv[arg++]);
						break;
					default :
						fprintf(stderr, "e event\n");
						usage();
				}
				break;
			default :
				fprintf(stderr, " option default\n");
				usage();
		}
	}
	if ((pars.mp.theta == 0) && (pars.mp.segsitesin == 0) && (pars.mp.treeflag == 0) && (pars.mp.timeflag == 0))
	{
		fprintf(stderr, " either -s or -t or -T option must be used. \n");
		usage();
		exit(EXIT_FAILURE);
	}

	sum = 0;
	for (i=0; i < pars.cp.npop; i++)
		sum += (pars.cp.config)[i];
	if (sum != pars.cp.nsam)
	{
		fprintf(stderr, " sum sample sizes != nsam\n");
		usage();
		exit(EXIT_FAILURE);
	}
}

void argcheck(int arg, int argc, char *argv[])
{
	if ((arg >= argc ) || (argv[arg][0] == '-'))
	{
		fprintf(stderr, "not enough arguments after %s\n", argv[arg-1]);
		fprintf(stderr, "For usage type: msmove<return>\n");
		exit(EXIT_FAILURE);
	}
}
	
int usage(void)
{
	fprintf(stderr, "usage: msmove nsam howmany \n");
	fprintf(stderr, "  Options: \n"); 
	fprintf(stderr, "\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
	fprintf(stderr, "\t -s segsites   ( fixed number of segregating sites)\n");
	fprintf(stderr, "\t -T          (Output gene tree.)\n");
	fprintf(stderr, "\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
	fprintf(stderr, "\t -r rho nsites     (rho here is 4Nc)\n");
	fprintf(stderr, "\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
	fprintf(stderr, "\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n"); 
	fprintf(stderr, "\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");      
	fprintf(stderr, "\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");    
	fprintf(stderr, "\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n"); 
	fprintf(stderr, "\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n"); 
	fprintf(stderr, "\t\t -n i size_i   (popi has size set to size_i*N0 \n");
	fprintf(stderr, "\t\t -g i alpha_i  (If used must appear after -M option.)\n"); 
	fprintf(stderr, "\t   The following options modify parameters at the time 't' specified as the first argument:\n");
	fprintf(stderr, "\t -eG t alpha  (Modify growth rate of all pop's.)\n");     
	fprintf(stderr, "\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");    
	fprintf(stderr, "\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n"); 
	fprintf(stderr, "\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n"); 
	fprintf(stderr, "\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");  
	fprintf(stderr, "\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");    
	fprintf(stderr, "\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
	fprintf(stderr, "\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");    
	fprintf(stderr, "\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
	fprintf(stderr, "\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
	fprintf(stderr, "\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
	fprintf(stderr, "\t\t  size, alpha and M are unchanged.)\n");
	fprintf(stderr, "\t -ev t i j x  ( Move lineages from pop i into pop j with probability x\n");
	fprintf(stderr, "\t\t  size, alpha and M are unchanged.\n");
	fprintf(stderr, "\t\t  If less than x lineages are present in pop i, then all lineages are moved.)\n"); 
	fprintf(stderr, "\t  -f filename     ( Read command line arguments from file filename.)\n");  
	fprintf(stderr, " See msdoc.pdf for explanation of these parameters.\n");
	exit(EXIT_FAILURE);
}

void addtoelist(struct devent *pt, struct devent *elist) 
{
	struct devent *plast, *pevent, *ptemp;

	pevent = elist;
	while ((pevent != 0) && (pevent->time <= pt->time))
	{
		plast = pevent;
		pevent = pevent->nextde;
	}
	ptemp = plast->nextde;
	plast->nextde = pt;
	pt->nextde = ptemp;
}

void free_eventlist(struct devent *pt, int npop)
{
	int pop;
	struct devent *next;
   
	while (pt != 0)
	{
		next = pt->nextde;
		if (pt->detype == 'a')
		{
			for (pop=0; pop < npop; pop++)
				free((pt->mat)[pop]);
			free(pt->mat);
		}
		free(pt);
		pt = next;
	}
}
	
int make_gametes(int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list)
{
	int tip, j, node;

	for (j=ns; j < ns+newsites; j++)
	{
		if (mfreq == 1)
			node = pickb(nsam, ptree, tt);
		else
			node = pickbmf(nsam, mfreq, ptree, tt);
		for (tip=0; tip < nsam; tip++)
		{
			if (tdesn(ptree, tip, node))
				list[tip][j] = STATE1;
			else
				list[tip][j] = STATE2;
		}
	}
	return 0;
}

double ttime(struct node *ptree, int nsam)
{
	int i;
	double t = (ptree+2*nsam-2)->time;

	for (i=nsam; i < 2*nsam-1; i++)
		t += (ptree+i)->time;
	
	return t;
}

double ttimemf(struct node *ptree, int nsam, int mfreq)
{
	int i;
	double t = 0;

	for (i=0; i < 2*nsam-2; i++)
		if (((ptree+i)->ndes >= mfreq) && ((ptree+i)->ndes <= nsam - mfreq))
			t += (ptree+(ptree+i)->abv)->time - (ptree+i)->time;
	
	return t;
}

void prtree(struct node *ptree, int nsam)
{
	int i;
	int *descl = (int*)malloc((size_t)(2*nsam-1)*sizeof(int));
	int *descr = (int*)malloc((size_t)(2*nsam-1)*sizeof(int));

	for (i=0; i < 2*nsam-1; i++)
		descl[i] = descr[i] = -1;
	for (i=0; i <  2*nsam-2; i++)
	{
		if (descl[(ptree+i)->abv] == -1)
			descl[(ptree+i)->abv] = i;
		else
			descr[(ptree+i)->abv] = i;
	}
	parens(ptree, descl, descr, 2*nsam-2);
	free(descl);
	free(descr);
}

void parens(struct node *ptree, int *descl, int *descr, int noden)
{
	double time;

	if (descl[noden] == -1)
		fprintf(stdout, "%d:%5.3lf", noden+1, (ptree+((ptree+noden)->abv))->time);
	else
	{
		fprintf(stdout, "(");
		parens(ptree, descl, descr, descl[noden]);
		fprintf(stdout, ",");
		parens(ptree, descl, descr, descr[noden]);
		if ((ptree+noden)->abv == 0)
			fprintf(stdout, ");\n"); 
		else
		{
			time = (ptree+(ptree+noden)->abv)->time - (ptree+noden)->time;
			fprintf(stdout, "):%5.3lf", time);
		}
	}
}

int pickb(int nsam, struct node *ptree, double tt)
{
	int i;
	double y = 0;
	double x = genrand_real3()*tt;

	for (i=0; i < 2*nsam-2; i++)
	{
		y += (ptree+(ptree+i)->abv)->time - (ptree+i)->time;
		if (y >= x)
			return i;
	}

	return 2*nsam-3;
}

int pickbmf(int nsam, int mfreq, struct node *ptree, double tt)
{
	int i;
	int lastbranch = 0;
	double x = genrand_real3()*tt;
	double y = 0;

	for (i=0; i < 2*nsam-2; i++)
	{
		if (((ptree+i)->ndes >= mfreq) && ((ptree+i)->ndes <= nsam - mfreq))
		{
			y += (ptree+(ptree+i)->abv)->time - (ptree+i)->time;
			lastbranch = i;
		}
		if (y >= x)
			return i;
	}

	return lastbranch;
}

int tdesn(struct node *ptree, int tip, int node)
{
	int k;

	for (k=tip; k < node; k = (ptree+k)->abv);

	if (k == node)
		return 1;
	else
		return 0;
}

int pick2(int n, int *i, int *j)
{
	*i = (int)(n*genrand_real3());
	while ((*j = (int)(n*genrand_real3())) == *i);
	return 0;
}

int ordran(int n, double pbuf[])
{
	ranvec(n, pbuf);
	order(n, pbuf);
	return 0;
}

int dbinom(int n, double p)
{
	int i, s;
	for (i=0, s=0; i < n; i++) 
		if (genrand_real3() <= p)
			++s;
	return s;
}

int mnmial(int n, int nclass, double p[], int rv[])
{
	int i, j;
	double x, s;

	for (i=0; i < nclass; i++)
		rv[i] = 0;
	for (i=0; i < n; i++)
	{
		x = genrand_real3();
		j = 0;
		s = p[0];
		while ((x > s) && (j < (nclass-1)))
			s += p[++j];
		rv[j]++;
	}
	return j;
}

int order(int n, double pbuf[])
{
	int gap, i, j;
	double temp;

	for (gap=n/2; gap > 0; gap /= 2)
		for (i=gap; i < n; i++)
			for (j=i-gap; j >= 0 && pbuf[j] > pbuf[j+gap]; j -= gap)
			{
				temp = pbuf[j];
				pbuf[j] = pbuf[j+gap];
				pbuf[j+gap] = temp;
			}
	return 0;
}

int ranvec(int n, double pbuf[])
{
	int i;
	for (i=0; i < n; i++)
		pbuf[i] = genrand_real3();
	return 0;
}

int poisso(double u)
{
	int i = 1;
	double cump, ru, p;

	if (u > 30.0)
	{
		i = (int)(0.5+gasdev(u,u));
		if (i < 0)
			return 0;
		else
			return i;
	} 
	ru = genrand_real3();
	p = exp(-u);
	if (ru < p)
		return 0;
	cump = p;
	while (ru > (cump += (p *= u/i)) )
		i++;
	return i;
}

double gasdev(double m, double v)
{
	static int iset = 0;
	static double gset;
	double fac, r, v1, v2;

	if (iset == 0)
	{
		do
		{
			v1 = 2.0*genrand_real3()-1.0;
			v2 = 2.0*genrand_real3()-1.0;
			r = v1*v1+v2*v2;
		} while (r >= 1.0);

		fac = sqrt(-2.0*log(r)/r);
		gset = v1*fac;
		iset = 1;

		return m+sqrt(v)*v2*fac;
	} 
	else
	{
		iset = 0;
		return m+sqrt(v)*gset;
	}
}
