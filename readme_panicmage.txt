MANUAL FOR panicmage

INSTALLATION  (LINUX)

extract all files to any directory

For panicmage the GNU Scientific Library (GSL) -- development package  (libgsl-dev) is needed, please install it.
In addition the following packages are needed:
ginac-tools
libginac-dev


compile with 

"g++ panicmage.c -lm -lgsl -lgslcblas -lcln -lginac -o panicmage"


you may now run panicmage from commandline with

./image [TREEFILE] [GFS_FILE] [INT] ... [OPTIONS]


[TREEFILE] 
  Has to be the tree your analysis is based on.
  Panicmage accepts only files in Newick format.
  The tree must contain the distances between the nodes and should be rooted.
  Only ultrametric trees will give meaningful results!
  I.e. the distances between each pair of individuals have to equal.
  The names/ID's for the individuals can be any integer, but 1,...,No. of individuals.
  It is recommend to use ID's between 1000 and 9999.
  A simple tree for 3 individuals might look like this:
  (1001:0.6,(1003:0.22,1002:0.22):0.38);
  For more infos have a look at http://en.wikipedia.org/wiki/Newick_format.

[GFS_FILE] 
  This file should contain the gene frequency data.
  Format should be simple text file and contain the number of genes per frequency class seperated by whitespaces.
  First number has to be the total number of unique genes in the population.
  Second number has to be the number of genes which are present in exactly two individuals.
  .....The last number has to be the number of genes present in each of the individuals (core genes).

[INT]
  Enter the number of individuals.

[OPTIONS]

  -g [FLOAT] -> set the number of generations up to the most recent common ancestor (MRCA) for your sample. 
		This optional information is necessary to compute the pansize and per generation rates.
		Number is treated in millions, so "-g 1.234" equals 1234000 generations up to the MRCA.

  -d  -> less details are printed

  -p  -> the tree (scaled to the estimated height) is printed

  -n  -> no neutraliy test is done. Normally panicmage will test whether neutrality holds or not. If this option is set the neutrality test is skipped.

  -s  -> test for sampling bias is done. Normally panicmage will not test whether the gene frequency spectrum is typical for a neutral evolving population,
	 where samples have not been drawn randomly. If this option is set sampling bias is tested.
  
  -q [INT] -> set the quantity of runs for the tests. Standard value is 10.000 runs.

  -a  -> all tests and estimations are skipped, only the theoretical results for a population with given parameters are shown. You have to set the parameters theta and rho
  
  -b  -> no rescaling of the supported tree is done


Setting the parameters to custom values:
  -t [FLOAT]  -> manually set the value for the parameter theta. If theta is set no estimation is done for any parameter.

  -r [FLOAT]  -> manually set the value for the parameter rho. If rho is set no estimation is done. The parameter rho has to be larger than zero.

  -c [FLOAT]  -> manually set the value for the number of core genes. If core is set no estimation is done for any parameter. 
Note:
if no further options than -t  [FLOAT] -r [FLOAT] -c [FLOAT] are set the neutrality test will be done with the given parameters.
Otherwise it depends on the additionally set options:
-n : the given parameters are only used for the theoretical results; same for -a
-s : the given parameters are used for the neutrality test. For the sampling bias test the given parameters are rescaled.
-n -s : the given paramters are used for the sampling bias test.



  

  
