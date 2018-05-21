<img src="source/panicmage.png" width="350">
  
# panicmage


Short summary of panicmage features
===================================

pan-genomes
-----------

The software panicmage simulates distributed genomes according to the
*Infinitely Many Genes Model* and can
estimate model parameters from observed data. Moreover panicmage
computes a p-value for the observed gene frequency spectrum and a given
genealogy in the *Infinitely Many Genes Model*. So far p-values for
neutral evolution can be computed. In addition panicmage can account for
sampling bias under neutral evolution. Given

-   a proxy for the true clonal genealogy of the sample

-   the gene frequency spectrum of the sample

-   the sample size n

-   the number of generations to the most recent common ancestor (MRCA)
    of the sample in millions (optional)

panicmage estimates gene gain and gene loss rates and the number of core
genes according to the *Infinitely Many Genes Model*. In addition
panicmage computes the p-value of the observed gene frequency spectrum
under neutral evolution. Small p-values hint at global selective
forces or population growth while very high p-values may appear if
gene transfer and/or recombination occur frequently. Finally panicmage
provides estimates for the following statistics in a neutral evolving
population:

-   the average number of genes per individual

-   the average total number of genes present in 2, n, 1000 and 10000
    individuals for a random clonal genealogy, resp.

-   the average total number of genes present in 2 or n individuals
    for the given clonal genealogy, resp.

-   the average number of new genes to be found in the next (n+1-th)
    sequenced individual

If in addition the number of generations to the MRCA is given panicmage
can also estimate:

-   The size of the pangenome, that is the total number of different
    genes present in the whole population.

-   The size of the persistant pangenome, that is the total number of
    different genes present in at least 1% of the whole population.

-   The per individual per generation gene gain rate

-   The per individual per gene per generation loss rate

CRISPR spacer arrays
--------------------

panicmage can be used to estimate the insertion and deletion rates of
CRISPR spacers, as outlined in @Baumdicker2018. If you intend to do so
use the spacer frequency spectrum instead of the gene frequency spectrum
and use the `-z` option.

Installing panicmage
====================

panicmage is written in C and C++ and designed for Linux. Currently
there is no manual for Windows or Mac.

dependecies
-----------

For panicmage the GNU Scientific Library (GSL) – development package
(libgsl0-dev) is needed, please install it, e.g. by typing

     sudo apt-get install libgsl0-dev

compiling
---------

In order to install panicmage extract panicmage.zip to any directory,
e.g. by typing

     unzip panicmage.zip

there is a binary file *panicmage* which might already work for your
system. Otherwise compile it again with:

    g++ panicmage.c -lm -lgsl -lgslcblas -o panicmage

In addition the flag -std=c++11 is required for gcc < 6.0. You may now
run panicmage from commandline. Type

    ./panicmage

to see basic usage.

panicmage requires at least 3 input parameters which have to be prompted
in the correct order.

    ./panicmage [TREEFILE] [GFS\_FILE] [INT] ... [OPTIONS]

Running panicmage
=================

[TREEFILE] Has to be the tree your analysis is based on
and represents the clonal genealogy of your sample. panicmage accepts
only files in Newick format. For more infos on the Newick format have a
look at <http://en.wikipedia.org/wiki/Newick_format>. The tree must
contain the distances between the nodes and should be rooted. Note that
only ultrametric bifurcating trees will give meaningful results! I.e.
the distances between each pair of individuals have to equal. However
panicmage does not check whether your tree is ultrametric, so be
careful. The names/ID’s for the individuals can be any string or
integer, but 1,...,no. of individuals. It is recommend to use ID’s
between 1000 and 9999. A simple tree for 3 individuals might look like
this: (101:0.6,(103:0.22,102:0.22):0.38);

[GFS\_FILE] This file should contain the gene frequency
data. The correct format is a simple text file which contains the number
of genes per frequency class seperated by whitespaces. To be precise,
the first number has to be the total number of unique genes in the
population. The second number has to be the number of genes which are
present in exactly two individuals. And so on such that the last number
is the number of genes present in all individuals (i.e. the number of
core genes).

[INT] Enter the number of individuals in your sample.
This number should equal the number of leafs in your treefile!

optional option

-g [FLOAT] Set the number of generations up to the most
recent common ancestor (MRCA) for your sample. Optional parameter but,
enables panicmage to estimate several statistics of interest (e.g. total
pangenome size, per generation gene gain and loss rates). Number is
treated in millions, so “-g 1.234” equals 1234000 generations up to the
MRCA.

For further options and instructions have a look at the ReadMe_Panicmage.pdf and the ReadMe_Panicmage.txt files.