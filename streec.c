#include "ms.h"

struct seg {
	int beg;
	int end;
	int desc;
};

struct chromo {
	int nseg;
	int pop;
	struct seg *pseg;
};

extern int flag;
int nchrom, begs, nsegs;
long int nlinks;
static int *nnodes = 0;  
double t, cleft, pc, lnpc;
static unsigned int seglimit = SEGINC;
static unsigned int maxchr;
struct node *ptree1, *ptree2;
static struct chromo *chrom = 0;
static struct segl *seglst = 0;

struct segl *segtre_mig(struct c_params *cp, int *pnsegs, int *migflag) 
{
	int i, j, k, dec, pop, pop2, c1, c2, ind, rchrom;
	int migrant, source_pop, *config, flagint;
	int eflag, cpop, ic, nsam, npop, nsites, *inconfig;
	double sum, x, ttemp, rft, clefta, tmin, p;
	double prec, cin, prect, mig, ran, coal_prob, rdum, arg;
	double r, f, rf, track_len, **migm;
	double *size, *alphag, *tlast;
	char ev;
	struct devent *nextevent;
	
	nsam = cp->nsam;
	npop = cp->npop;
	nsites = cp->nsites;
	inconfig = cp->config;
	r = cp->r;
	f = cp->f;
	track_len = cp->track_len;
	migm = (double**)malloc((size_t)npop*sizeof(double*));
	for (i=0; i < npop; i++)
	{
		migm[i] = (double*)malloc((size_t)npop*sizeof(double));
		for (j=0; j < npop; j++)
			migm[i][j] = (cp->mig_mat)[i][j];
	}
	nextevent = cp->deventlist;
	
	if (chrom == NULL)
	{
		maxchr = nsam+20;
		chrom = (struct chromo*)malloc((size_t)maxchr*sizeof(struct chromo));
		if (chrom == NULL)
			perror("malloc error. segtre");
	}
	if (nnodes == NULL)
	{
		nnodes = (int*)malloc((size_t)seglimit*sizeof(int));
		if (nnodes == NULL)
			perror("malloc error. segtre_mig");
	}
	if (seglst == NULL)
	{
		seglst = (struct segl*)malloc((size_t)seglimit*sizeof(struct segl));
		if (seglst == NULL)
			perror("malloc error. segtre_mig.c 2");
	}
	config = (int*)malloc((size_t)(npop+1)*sizeof(int));
	if (config == NULL)
		perror("malloc error. segtre.");
	size = (double*)malloc((size_t)npop*sizeof(double));
	if (size == NULL)
		perror("malloc error. segtre.");
	alphag = (double*)malloc((size_t)npop*sizeof(double));
	if (alphag == NULL)
		perror("malloc error. segtre.");
	tlast = (double*)malloc((size_t)npop*sizeof(double));
	if (tlast == NULL)
		perror("malloc error. segtre.");
	for (pop=0; pop < npop; pop++)
	{
		config[pop] = inconfig[pop];
		size[pop] = (cp->size)[pop];
		alphag[pop] = (cp->alphag)[pop];
		tlast[pop] = 0;
	}
	for (pop=ind=0; pop < npop; pop++)
		for (j=0; j < inconfig[pop]; j++, ind++)
		{
			chrom[ind].nseg = 1;
			if (!(chrom[ind].pseg = (struct seg*)malloc(sizeof(struct seg))))
				perror("calloc error. se1");
			(chrom[ind].pseg)->beg = 0;
			(chrom[ind].pseg)->end = nsites-1;
			(chrom[ind].pseg)->desc = ind;
			chrom[ind].pop = pop;
		}
	seglst[0].beg = 0;
	if (!(seglst[0].ptree = (struct node*)calloc((size_t)(2*nsam), sizeof(struct node))))
		perror("calloc error. se2");
	nnodes[0] = nsam-1;
	nchrom = nsam;
	nlinks = (long int)(nsam*(nsites-1));
	nsegs = 1;
	t = 0;
	r /= nsites-1;
	if (f > 0)
		pc = (track_len-1.0)/track_len;
	else
		pc = 1.0;
	lnpc = log(pc);
	cleft = nsam * (1.0-pow(pc, (double)(nsites-1)));
	if (r > 0)
		rf = r*f;
	else
		rf = f/(nsites-1);
	rft = rf*track_len;
	flagint = 0;

	while (nchrom > 1)
	{
		prec = nlinks*r;
		cin = nlinks*rf;
		clefta = cleft*rft;
		prect = prec+cin+clefta;
		mig = 0;
		for (i=0; i < npop; i++)
			mig += config[i]*migm[i][i];
		if ((npop > 1) && (mig == 0) && (nextevent == 0))
		{
			i = 0;
			for (j=0; j < npop; j++) 
				if (config[j] > 0)
					i++;
			if (i > 1)
			{
				fprintf(stderr, " Infinite coalescent time. No migration.\n");
				exit(EXIT_FAILURE);
		   }
		}
		eflag = 0;
		if (prect > 0)
		{
			while ((rdum = genrand_real3()) == 0);
			ttemp = -1.0*log(rdum)/prect;
			if ((eflag == 0) || (ttemp < tmin))
			{
				tmin = ttemp;
				ev = 'r';
				eflag = 1;
			}
		}
		if (mig > 0)
		{
			while ((rdum = genrand_real3()) == 0);
			ttemp = -1.0*log(rdum)/mig;
			if ((eflag == 0) || (ttemp < tmin))
			{
				tmin = ttemp;
				ev = 'm';
				eflag = 1;
			}
		}
		for (pop=0; pop < npop; pop++)
		{
			coal_prob = (double)config[pop]*(config[pop]-1);
			if (coal_prob > 0)
			{
				while ((rdum = genrand_real3()) == 0);
				if (alphag[pop] == 0)
				{
					ttemp = -1.0*log(rdum)*size[pop]/coal_prob;
					if ((eflag == 0) || (ttemp < tmin))
					{
						tmin = ttemp;
						ev = 'c';
						eflag = 1;
						cpop = pop;
					}
				}
				else
				{
					arg = 1.0-alphag[pop]*size[pop]*exp(-1.0*alphag[pop]*(t-tlast[pop]))*log(rdum)/coal_prob;
					if (arg > 0)
					{
						ttemp = log(arg)/alphag[pop];
						if ((eflag == 0) || (ttemp < tmin))
						{
							tmin = ttemp;
							ev = 'c';
							eflag = 1;
							cpop = pop;
						}
					}
				}
			}
		}
		if ((eflag == 0) && (nextevent == 0))
		{
			fprintf(stderr, " Infinite time to next event. Negative growth rate in last time interval or non-communicating subpops.\n");
			exit(EXIT_FAILURE);
		}
		if (((eflag == 0) && (nextevent != 0)) || ((nextevent != 0) && ((t+tmin) >= nextevent->time)))
		{
			t = nextevent->time;
			switch (nextevent->detype)
			{
				case 'N' :
					for (pop=0; pop < npop; pop++)
					{
						size[pop]= nextevent->paramv;
						alphag[pop] = 0;
					}
					nextevent = nextevent->nextde;
					break;
				case 'n' :
					size[nextevent->popi] = nextevent->paramv;
					alphag[nextevent->popi] = 0;
					nextevent = nextevent->nextde;
					break;
				case 'G' :
					for (pop=0; pop < npop; pop++)
					{
						size[pop] = size[pop]*exp(-1.0*alphag[pop]*(t-tlast[pop]));
						alphag[pop]= nextevent->paramv;
						tlast[pop] = t;
					}
					nextevent = nextevent->nextde;
					break;
				case 'g' :
					pop = nextevent->popi;
					size[pop] = size[pop]*exp(-1.0*alphag[pop]*(t-tlast[pop]));
					alphag[pop]= nextevent->paramv;
					tlast[pop] = t;
					nextevent = nextevent->nextde;
					break;
				case 'M' :
					for (pop=0; pop < npop; pop++)
						for (pop2=0; pop2 < npop; pop2++)
							migm[pop][pop2] = (nextevent->paramv)/(npop-1.0);
					for (pop = 0; pop < npop; pop++)
						migm[pop][pop] = nextevent->paramv;
					nextevent = nextevent->nextde;
					break;
				case 'a' :
					for (pop=0; pop < npop; pop++)
						for (pop2=0; pop2 < npop; pop2++)
							migm[pop][pop2] = (nextevent->mat)[pop][pop2];
					nextevent = nextevent->nextde;
					break;
				case 'm' :
					i = nextevent->popi;
					j = nextevent->popj;
					migm[i][i] += nextevent->paramv - migm[i][j];
					migm[i][j]= nextevent->paramv;
					nextevent = nextevent->nextde;
					break;
				case 'j' :
					i = nextevent->popi;
					j = nextevent->popj;
					config[j] += config[i];
					config[i] = 0;
					for (ic=0; ic < nchrom; ic++)
						if (chrom[ic].pop == i)
							chrom[ic].pop = j;
					for (k=0; k < npop; k++)
					{
						if (k != i)
						{
							migm[k][k] -= migm[k][i];
							migm[k][i] = 0;
						}
					}
					nextevent = nextevent->nextde;
					break;
				case 'v' :
					i = nextevent->popi;
					j = nextevent->popj;
					/* start new code here-- old code below */
					config[i] = 0;
					for (ic=0; ic < nchrom; ic++)
					{
						if (chrom[ic].pop == i)
						{
							if (genrand_real3() > nextevent->pmove)
								config[i]++;
							else
							{
								chrom[ic].pop = j;
								config[j]++;
								*migflag = 1;
							}
						}
					}
					nextevent = nextevent->nextde;
					break;
				case 's' :
					i = nextevent->popi;
					p = nextevent->paramv;
					npop++;
					config = (int*)realloc(config, (size_t)npop*sizeof(int));
					size = (double*)realloc(size, (size_t)npop*sizeof(double));
					alphag = (double*)realloc(alphag, (size_t)npop*sizeof(double));
					tlast = (double*)realloc(tlast, (size_t)npop*sizeof(double));
					tlast[npop-1] = t;
					size[npop-1] = 1.0;
					alphag[npop-1] = 0;
					migm = (double**)realloc(migm, (size_t)npop*sizeof(double*));
					for (j=0; j < npop-1; j++)
						migm[j] = (double*)realloc(migm[j], (size_t)npop*sizeof(double));
					migm[npop-1] = (double*)malloc((size_t)npop * sizeof(double));
					for (j=0; j < npop; j++)
						migm[npop-1][j] = migm[j][npop-1] = 0;
					config[npop-1] = 0;
					config[i] = 0;
					for (ic=0; ic < nchrom; ic++)
					{
						if (chrom[ic].pop == i)
						{
							if (genrand_real3() < p)
								config[i]++;
							else
							{
								chrom[ic].pop = npop-1;
								config[npop-1]++;
							}
						}
					}
					nextevent = nextevent->nextde;
					break;
			}
		}
		else
		{
			t += tmin;
			if (ev == 'r')
			{
				if ((ran = genrand_real3()) < prec/prect)
				{
					rchrom = re(nsam);
					config[chrom[rchrom].pop] += 1;
				}
				else if (ran < (prec+clefta)/prect)
				{
					rchrom = cleftr(nsam);
					config[chrom[rchrom].pop] += 1;
				}
				else
				{
					rchrom = cinr(nsam, nsites);
					if (rchrom >= 0)
						config[chrom[rchrom].pop] += 1;
				}
			}
			else if (ev == 'm')
			{
				x = mig*genrand_real3();
				sum = 0;
				for (i=0; i < nchrom; i++)
				{
					sum += migm[chrom[i].pop][chrom[i].pop];
					if (x < sum)
						break;
				}
				migrant = i;
				x = genrand_real3()*migm[chrom[i].pop][chrom[i].pop];
				sum = 0;
				for (i=0; i < npop; i++)
				{
					if (i != chrom[migrant].pop)
					{
						sum += migm[chrom[migrant].pop][i];
						if (x < sum)
							break;
					}
				}
				source_pop = i;
				config[chrom[migrant].pop] -= 1;
				config[source_pop] += 1;
				chrom[migrant].pop = source_pop;
			}
			else
			{
				pick2_chrom(cpop, config, &c1, &c2);
				dec = ca(nsam, nsites, c1, c2);
				config[cpop] -= dec;
			}
		}
	}  
	*pnsegs = nsegs;
	free(config); 
	free(size);
	free(alphag);
	free(tlast);
	for (i=0; i < npop; i++)
		free(migm[i]);
	free(migm);
	return seglst;
}

int re(int nsam)
{
	int el, lsg, lsgm1, ic, is;
	int spot = (int)(nlinks * genrand_real3()) + 1;
	struct seg *pseg;

	for (ic=0; ic < nchrom; ic++)
	{
		lsg = chrom[ic].nseg;
		lsgm1 = lsg-1;
		pseg = chrom[ic].pseg;
		el = ((pseg+lsgm1)->end) - (pseg->beg);
		if (spot <= el)
			break;
		spot -= el;
	}
	is = pseg->beg+spot-1;
	xover(nsam, ic, is);
	return ic;	
}

int cleftr(int nsam)
{
	int is;
	int ic = -1;
	double x, len;
	double sum = 0;
	struct seg *pseg;

	while ((x = cleft * genrand_real3()) == 0);
	while (sum < x)
		sum +=  1.0-pow(pc, links(++ic));
	pseg = chrom[ic].pseg;
	len = links(ic);
	is = (int)(pseg->beg+floor(1.0+log(1.0-(1.0-pow(pc, len))*genrand_real3())/lnpc)-1);
	xover(nsam, ic, is);
	return ic;
}

int cinr(int nsam, int nsites)
{
	int len, el, lsg, lsgm1, ic, is, endic;
	int spot = (int)(nlinks * genrand_real3()) + 1;
	struct seg *pseg;
	
	for (ic=0; ic < nchrom; ic++)
	{
		lsg = chrom[ic].nseg;
		lsgm1 = lsg-1;
		pseg = chrom[ic].pseg;
		el = ((pseg+lsgm1)->end)-(pseg->beg);
		if (spot <= el)
			break;
		spot -= el;
	}
	is = pseg->beg+spot-1;
	endic = (pseg+lsgm1)->end;
	xover(nsam, ic, is);
	len = (int)floor(1.0+log(genrand_real3())/lnpc);
	if (is+len >= endic)
		return ic;  
	if (is+len < (chrom[nchrom-1].pseg)->beg)
	{
		ca(nsam, nsites, ic, nchrom-1);
		return -1;
	}
	xover(nsam, nchrom-1, is+len);
	ca(nsam, nsites, ic, nchrom-1);
	return ic;	

}

int xover(int nsam, int ic, int is)
{
	int i, k, lsgm1, newsg, jseg, in;
	int lsg = chrom[ic].nseg;	
	struct seg *pseg = chrom[ic].pseg;
	struct seg *pseg2;
	double len = (pseg+lsg-1)->end - pseg->beg;
	
	cleft -= 1.0-pow(pc, len);

	for (jseg=0; is >= (pseg + jseg)->end; jseg++);
	if (is >= (pseg+jseg)->beg)
		in = 1;
	else
		in = 0;
	newsg = lsg-jseg;
	nchrom++;
	if ((unsigned int)nchrom >= maxchr)
	{
		maxchr += 20;
		chrom = (struct chromo*)realloc(chrom, (size_t)maxchr*sizeof(struct chromo));
		if (chrom == NULL)
			perror("malloc error. segtre2");
	}
	if (!(pseg2 = chrom[nchrom-1].pseg = (struct seg*)calloc((size_t)newsg, sizeof(struct seg))))
		perror(" alloc error. re1");
	chrom[nchrom-1].nseg = newsg;
	chrom[nchrom-1].pop = chrom[ic].pop;
	pseg2->end = (pseg+jseg)->end;
	if (in)
	{
		pseg2->beg = is+1;
		(pseg+jseg)->end = is;
	}
	else
		pseg2->beg = (pseg+jseg)->beg;
	pseg2->desc = (pseg+jseg)->desc;
	for (k=1; k < newsg; k++)
	{
		(pseg2+k)->beg = (pseg+jseg+k)->beg;
		(pseg2+k)->end = (pseg+jseg+k)->end;
		(pseg2+k)->desc = (pseg+jseg+k)->desc;
	}
	lsg = chrom[ic].nseg = lsg-newsg+in;
	lsgm1 = lsg-1;
	nlinks -= pseg2->beg - (pseg+lsgm1)->end;
	len = (pseg+lsgm1)->end - (pseg->beg);
	cleft += 1.0-pow(pc, len);
	len = (pseg2+newsg-1)->end - pseg2->beg;
	cleft += 1.0-pow(pc, len);
	if (!(chrom[ic].pseg = (struct seg*)realloc(chrom[ic].pseg, (size_t)lsg*sizeof(struct seg))))
		perror("realloc error. re2");
	if (in)
	{
		begs = pseg2->beg;
		for (i=0, k=0; (k < nsegs-1) && (begs > seglst[seglst[i].next].beg-1); i=seglst[i].next, k++);
		if (begs != seglst[i].beg)
		{
			if ((unsigned int)nsegs >= seglimit)
			{
				seglimit += SEGINC;
				nnodes = (int*)realloc(nnodes, (size_t)seglimit*sizeof(int));
				if (nnodes == NULL)
					perror("realloc error. 1. segtre_mig.c");
				seglst = (struct segl*)realloc(seglst,(size_t)seglimit*sizeof(struct segl));
				if (seglst == NULL)
					perror("realloc error. 2. segtre_mig.c");
			}
			seglst[nsegs].next = seglst[i].next;
			seglst[i].next = nsegs;
			seglst[nsegs].beg = begs;
			if (!(seglst[nsegs].ptree = (struct node*)calloc((size_t)(2*nsam), sizeof(struct node))))
				perror("calloc error. re3.");
			nnodes[nsegs] = nnodes[i];
			ptree1 = seglst[i].ptree;
			ptree2 = seglst[nsegs].ptree;
			nsegs++;
			for (k=0; k <= nnodes[i]; k++)
			{
				(ptree2+k)->abv = (ptree1+k)->abv;
				(ptree2+k)->time = (ptree1+k)->time;
			}
		}
	}
	return ic;
}

int ca(int nsam, int nsites, int c1, int c2)
{
	int yes1, yes2, seg, start, end, desc, k;
	int seg1 = 0;
	int seg2 = 0;
	int tseg = -1;
	struct seg *pseg;
	struct node *ptree;

	if (!(pseg = (struct seg*)calloc((size_t)nsegs, sizeof(struct seg)))) 
		perror("alloc error.ca1");

	for (seg=0, k=0; k < nsegs; seg = seglst[seg].next, k++)
	{
		start = seglst[seg].beg;
		yes1 = isseg(start, c1, &seg1);
		yes2 = isseg(start, c2, &seg2);
		if (yes1 || yes2)
		{
			tseg++;
			(pseg+tseg)->beg = seglst[seg].beg;
			end = k < nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1;
			(pseg+tseg)->end = end;
			if (yes1 && yes2)
			{
				nnodes[seg]++;
				if (nnodes[seg] >= 2*nsam-2)
					tseg--;
				else
					(pseg+tseg)->desc = nnodes[seg];
				ptree = seglst[seg].ptree;
				desc = (chrom[c1].pseg+seg1)->desc;
				(ptree+desc)->abv = nnodes[seg];
				desc = (chrom[c2].pseg+seg2)-> desc;
				(ptree+desc)->abv = nnodes[seg];
				(ptree+nnodes[seg])->time = t;
			}
			else
				(pseg+tseg)->desc = yes1 ? (chrom[c1].pseg+seg1)->desc : (chrom[c2].pseg+seg2)->desc;
		}
	}
	nlinks -= links(c1);
	cleft -= 1.0-pow(pc, (double)links(c1));
	free(chrom[c1].pseg);
	if (tseg < 0)
	{
		free(pseg);
		chrom[c1].pseg = chrom[nchrom-1].pseg;
		chrom[c1].nseg = chrom[nchrom-1].nseg;
		chrom[c1].pop = chrom[nchrom-1].pop;
		if (c2 == nchrom-1)
			c2 = c1;
		nchrom--;
	}
	else
	{
		if(!(pseg = (struct seg*)realloc(pseg, (size_t)(tseg+1)*sizeof(struct seg))))
			perror(" realloc error. ca1");
		chrom[c1].pseg = pseg;
		chrom[c1].nseg = tseg+1;
		nlinks += links(c1);
	   	cleft += 1.0-pow(pc, (double)links(c1));
	}
	nlinks -= links(c2);
	cleft -= 1.0-pow(pc, (double)links(c2));
	free(chrom[c2].pseg);
	chrom[c2].pseg = chrom[nchrom-1].pseg;
	chrom[c2].nseg = chrom[nchrom-1].nseg;
	chrom[c2].pop = chrom[nchrom-1].pop;
	nchrom--;
	if (tseg < 0)
		return 2;
	else 
		return 1;
}

int isseg(int start, int c, int *psg)
{
	int ns = chrom[c].nseg;
	struct seg *pseg = chrom[c].pseg;

	for(; ((*psg) < ns) && ((pseg+(*psg))->beg <= start); ++(*psg))
		if ((pseg+(*psg))->end >= start)
			return 1;
	return 0;
}
	
int pick2_chrom(int pop, int config[], int *pc1, int *pc2)
{
	int c1, c2, cs, cb, i, count;
	
	pick2(config[pop], &c1, &c2);
	cs = c1 > c2 ? c2 : c1;
	cb = c1 > c2 ? c1 : c2;
	i = count = 0;
	for(;;)
	{
		while (chrom[i].pop != pop) 
			i++;
		if (count == cs)
			break;
		count++;
		i++;
	}
	*pc1 = i;
	i++;
	count++;
	for(;;)
	{
		while (chrom[i].pop != pop)
			i++;
		if (count == cb)
			break;
		count++;
		i++;
	}
	*pc2 = i;
	return 0;
}
	
int links(int c)
{
	int ns = chrom[c].nseg-1;
	return (chrom[c].pseg+ns)->end - (chrom[c].pseg)->beg;
}

