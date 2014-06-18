#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#define BASE 10
#define LINE_SIZE 1000
#define BINOM(x) ( (x) * ((x) - 1) / 2 )

// Function prototypes
int pairDiffs(int, int, int*, char**, int*);
int addToBin(int, int*, int*);
int growMatrix(int, unsigned int, char**);
char **allocMatrix(int, int);

// Globally scoped variables
unsigned int maxSites = 1000;
unsigned int maxBins = 1000;
unsigned int maxDiffs = 0;
int between_flag = 0;

int main(int argc, char *argv[])
{
	int i;                          // Generic iterator
	int nSamples;                   // Total number of samples
	int nReplicates;                // Number of simulated replicates
	int nPops;                      // Total number of populations sampled
	int *config;                    // Array of population sample sizes
	int segSites;                   // Total number of segregating sites
	int count = 0;                  // Count of current replicate data set
	int nPairs;                     // Number of pairwise comparisons
	int *diffs;                     // Array of pairwise differences
	int *bin;                       // Bins for storing counts
	double *posit;                  // Relative position of each site
	char **dna;                     // Matrix of simulated DNA sequences
	char line[LINE_SIZE];           // Input line
	char slashline[LINE_SIZE];      // Contents of the slash line
	char cmdline[LINE_SIZE];        // The ms command line
	char dum[20];                   // A dummy string
	char astr[100];                 // Another dummy string
	char *mscanline;
	int opt = 0;
	int long_index = 0;
	static struct option long_options[] = {
	  {"between", no_argument,  0,  'b' },
	  {0,         0,            0,   0  }
	};

	// Process command line arguments
	while ( (opt = getopt_long(argc, argv, "b", long_options, &long_index)) != -1 ) 
	{
		switch(opt)
		{
			case 'b' :
				between_flag = 1;
				break;
			default :
				break;
		}
	}
	
	
	// Get ms command line and pull down parameters
	fgets(cmdline, LINE_SIZE, stdin);
	fprintf(stdout, "# %s", cmdline);
	sscanf(cmdline, " %s  %d %d", dum, &nSamples, &nReplicates);
	mscanline = strstr(cmdline, "-I");

	// Get population sample configuration
	if ( mscanline && between_flag )
	{
		mscanline = strchr(mscanline, ' ');
		nPops = (int)strtol(mscanline, &mscanline, BASE);
		if ( nPops > 2 )
		{
			fprintf(stderr, "mismatch currently only works with two pops\n");
			exit(1);
		}
		if ( !(config = (int *)malloc((size_t)nPops * sizeof(int))) )
			perror("malloc error for config");
		for (i = 0; i < nPops; ++i)
			config[i] = (int)strtol(mscanline, &mscanline, BASE);
		nPairs = config[0] * config[1];
	} else if ( !mscanline && between_flag )
	{
		fprintf(stderr, "--between switch used, but only one pop found\n");
		exit(1);
	} else if ( mscanline && !between_flag )
	{
		fprintf(stderr, "two pops detected, but no --between switch used.\n");
		exit(1);
	} else
	{
		nPops = 1;
		config = NULL;
		nPairs = BINOM(nSamples);
	}
	
	// Gobble line with rng seed
	fgets(line, LINE_SIZE, stdin);
	
	// Allocate memory for variables
	dna = allocMatrix(nSamples, maxSites + 1);
	posit = (double *)malloc((size_t)maxSites * sizeof(double));
	if ( posit == NULL )
		perror("Allocation failure for posit");
	diffs = (int *)malloc((size_t)nPairs * sizeof(int));
	if ( diffs == NULL )
		perror("Allocation failure for diffs");
	bin = (int *)malloc((size_t)1000 * sizeof(int));
	if ( bin == NULL )
		perror("Allocation failure for bin");
	memset(bin, 0, (size_t)maxBins * sizeof(int));
	
	// Iterate through replicates
	while ( nReplicates - count++ )
	{
	
		// Move to the next replicate data set
		do 
		{
			if ( fgets(line, LINE_SIZE, stdin) == NULL ) 
			{
				if ( feof(stdin) )
					return 0;
				else if ( ferror(stdin) )
					perror("Error reading file");
			}
			if ( line[0] == '/' )
				strcpy(slashline, line + 2);
		} while ( (line[0] != 's') && (line[0] != 'p') );
		
		// Get the number of segregating sites for the replicate
		sscanf(line, "  segsites: %d", &segSites);
		
		// Reallocate memory for the sequence matrix, if needed
		if ( segSites >= maxSites )
		{
			maxSites = segSites + 10;
			posit = (double *)realloc(posit, (size_t)maxSites * sizeof(double));
			if (posit == NULL)
				perror("realloc error for posit");
			growMatrix(nSamples, maxSites, dna);
		}
		
		// If variation is present-- run analyses
		if ( segSites > 0 )
		{
			fscanf(stdin, " %s", astr);
			
			// Read in relative site positions
			for (i = 0; i < segSites; ++i)
				fscanf(stdin, " %lf", posit + i);
				
			// Read in sequence matrix
			for (i = 0; i < nSamples; ++i)
				fscanf(stdin, " %s", dna[i]);
				
			// Calculate pairwise differences
			pairDiffs(nSamples, segSites, config, dna, diffs);

			// Add to histogram
			addToBin(nPairs, diffs, bin);
		}
		// End of replicate
	}		

	// Generate program output
	for (i = 0; i <= maxDiffs; ++i)
		fprintf(stdout, "%d\t%1.8E\n", i, (double)bin[i] / (nPairs * nReplicates));

	free(diffs);
	free(bin);
	free(posit);
	for (i = 0; i < nSamples; ++i)
		free(dna[i]);
	free(dna);

	return 0;
}


// Function to calculate pairwise differences between chromosomes

int pairDiffs(int n, int s, int *config, char **dna, int *d)
{
	int i, j, k;
	int idx;
	int e1, e2;

	// Initialize differences to zero
	if (config) {
		memset(d, 0, (size_t)(config[0] * config[1]) * sizeof(int));
		e1 = config[0];
	} else
	{
		memset(d, 0, (size_t)BINOM(n) * sizeof(int));
		e1 = n - 1;
	}
		
	// Count the number of pairwise differences
	for (i = 0, idx = 0; i < e1; ++i)
	{
		if (config)
			e2 = config[0];
		else
			e2 = i + 1;
		for (j = e2; j < n; ++j)
		{
			for (k = 0; k < s; ++k)
				if ( dna[i][k] != dna[j][k] )
					++d[idx];
			++idx;
		}
	}

	return 0;
}


// Add pairwise differences per replicate to histogram bins

int addToBin(int np, int *diffs, int *bin)
{
	int i;

	for (i = 0; i < np; ++i) 
	{
		if ( diffs[i] > maxDiffs )
			maxDiffs = (unsigned int)diffs[i];
		if ( diffs[i] >= maxBins ) 
		{
			unsigned int oldBins = maxBins;
			maxBins = diffs[i] + 10;
			bin = (int *)realloc(bin, (size_t)maxBins * sizeof(int));
			if ( bin == NULL )
				perror("Realloc failure for bin");
			memset(bin + oldBins, 0, (size_t)(maxBins - oldBins) * sizeof(int));
		} 
		++bin[diffs[i]];
	}
	
	return 0;
}


// Function to allocate memory for a two dimensional matrix

char **allocMatrix(int n, int l)
{
	int i;
	char **dna;
	
	if ( !(dna = (char **)malloc((size_t)n * sizeof(char *))) )
		perror("alloc error in allocMatrix()");
	for (i = 0; i < n; ++i) 
	{
		if ( !(dna[i] = (char *)malloc((size_t)l * sizeof(char))) )
			perror("alloc error in allocMatrix().2");
	}
	
	return dna;
}


// Function to reallocate the size of a two dimensional matrix

int growMatrix(int n, unsigned int nMax, char **dna)
{
	int i;
	
	maxSites = nMax;
	for (i = 0; i < n; ++i)
	{
		dna[i] = (char *)realloc(dna[i], (size_t)maxSites * sizeof(char));
		if ( dna[i] == NULL )
			perror("realloc error in growMatrix()");
	}
	
	return 0;
}                        
	
