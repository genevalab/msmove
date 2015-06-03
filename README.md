msmove
======

A modified version of Hudson's coalescent simulator ms allowing for finer control and tracking of migrant genealogies.

####To build msmove
```
make
```
The binary will be located in a directory named gccRelease.

####To build msiso
```
gcc -O2 -Wall -o msiso msiso.c -lm
```
======

#### Usage

All commands and options associated with Hudson's ms are available in msmove, see the ms documentation for more details on these commands. msmove includes an additional command which allows for control of the timing and rate of migration. Simulations in which migration has occurred are reported along with the standard ms output.

```
msmove nsam howmany [other ms options options] -ev t i j x
```
This command moves lineages at time t from pop i into pop j with probability x. Population size, alpha, and M are unchanged. If less than i*x lineages are present in pop i, then all lineages are moved.

The output is identical to ms output with one addition. If an introgression event occurred in a particular replicate a '*' is added after the '//' that indicates the start of a set of results. 

For example:
```
msmove 20 100 -t 10 -I 2 10 10 -r 100 10000 -ej 1 2 1 -ev  0.01 2 1 0.01
```

should produce around 10-20 replicates for which migration has occurred

======
####Citation

Users should cite the original description of ms. The manual and source for ms are available [here](https://webshare.uchicago.edu/xythoswfs/webview/fileManager?stk=7EB767DDE194BCF&entryName=%2Fusers%2Frhudson1%2FPublic%2Fms.folder&msgStatus=).   
   * Hudson, R. R. (2002) Generating samples under a Wright-Fisher neutral model. Bioinformatics. 18:337-8.

The citation and  DOI number for msmove is:   
  * Garrigan, Daniel; Geneva, Anthony J (2014): msmove. figshare. http://dx.doi.org/10.6084/m9.figshare.1060474
