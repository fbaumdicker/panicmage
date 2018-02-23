/***** panicmage.c ***********************************************
 * Description:  panicmage estimates the parameters of the infinitely many genes model for given datasets
 * and performs statistical tests on neutral genome evolution and sampling bias.
 * Expected values for various values of interest may be computed (Pangenome size, expected number of new genes in the next individual etc.)
 * How to cite: if you use panicmage please cite
 * Baumdicker, F., Hess, W. R., & Pfaffelhuber, P. (2010). The diversity of a distributed genome in bacterial populations. The Annals of Applied Probability, 20(5), 1567–1606.
 * or
 * Baumdicker, F., Hess, W. R., & Pfaffelhuber, P. (2012). The infinitely many genes model for the distributed genome of bacteria. Genome Biology and Evolution, 4(4), 443–456. doi:10.1093/gbe/evs016
 * this is Version 1.1
 * Author: Franz Baumdicker, baumdicker@stochastik.uni-freiburg.de
 * Author: Peter Pfaffelhuber, pp@stochastik.uni-freiburg.de
 *****************************************************************/
// compile with g++ panicmage.c -lm -lgsl -lgslcblas -lcln -lginac -o panicmage


#include <gsl/gsl_multimin.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <ginac/ginac.h>
using namespace GiNaC;
#include "source/treestructure.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "source/treefunctions.h"
#include "source/treeprints.h"
#include "source/readingtree.h"
#include "source/treestatistics.h"
#include "source/theoreticalvalues.h"
#include "source/treesamplingbias.h"
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <iostream>
#include "source/treesymbolic.h"


int g_includecoreflag = 0;









//parse the commandline for options


  int get_args(int argc, char** argv, int* notest_flag, int* samplingbias_flag,  float* theta_value, float* rho_value, float* core_value, int* estimate_flag, int* details_flag, int* skipall_flag, int* runs_input, int* printtree_flag, int* pansize_flag, float* millgenerationstoMRCA, int* scale_flag)
{
    int i;

    /* Start at i = 4 to skip the command name and the required inputs. */

    for (i = 4; i < argc; i++) {

	/* Check for a switch (leading "-"). */

	if (argv[i][0] == '-') {

	    /* Use the next character to decide what to do. */

	    switch (argv[i][1]) {

		case 'n':	*notest_flag = 1;
				break;
				
		case 's':	*samplingbias_flag = 1;
				break;
				
		case 't':	if (i+1 >= argc) return -1;
				*theta_value = atof(argv[++i]);
				*estimate_flag = 0;
				break;
				
		case 'r':	if (i+1 >= argc) return -1;
				*rho_value = atof(argv[++i]);
				*estimate_flag = 0;
				break;
				
		case 'c':	if (i+1 >= argc) return -1;
				*core_value = atof(argv[++i]);
				*estimate_flag = 0;
				break;
				
		case 'q':	if (i+1 >= argc) return -1;
				*runs_input = atoi(argv[++i]);
				break;
				
		case 'd':	*details_flag = 0;
				break;
				
		case 'a':	*skipall_flag = 1;
				break;
								
		case 'b':	*scale_flag = 0;
				break;
				
		case 'z':	g_includecoreflag = 1;
				break;
				
		case 'p':	*printtree_flag = 1;
				break;
				
		case 'g':	if(i+1 >= argc) return -1;
				*millgenerationstoMRCA = atof(argv[++i]);
				*pansize_flag = 1;
				break;

		default:	fprintf(stderr,
				"Unknown option %s\n", argv[i]);
	    }
	}
	else fprintf(stderr,"Options should start with \"-\". You wrote: %s\n", argv[i]);
    }
    return 0;
}




// define global random number generator
const gsl_rng_type * T;
gsl_rng * r;





int main( int argc, char* argv[]){
    
 
    // this defines the global variable for the paremeter rho (the gene loss rate)
    symbol x("rhoS");
   
  /* create a generator chosen by the 
          environment variable GSL_RNG_TYPE */
   gsl_rng_env_setup();
     
   T = gsl_rng_default;
   r = gsl_rng_alloc (T);
  
  
  
// reading the parameters from commandline
  
//   printf("Progname is:\t %s\n", argv[0]);
//   printf("Number of arguments is:\t %d\n", argc);
  if (argc < 4){
    printf("missing parameters, expected input is:\n./panicmage [TREEFILE] [GFS_FILE] [INT] ... [OPTIONS]\nfor more info have a look at the readme files\n");
    return -1;
  }

printf("Filename Newick-Tree:\t\t %s\n", argv[1]);
printf("Filename Gene Frequency Data:\t %s\n", argv[2]);
  
//////////////////////*reading the data *////////////////////////////////////////
    char dateiname[1000];
    char treeline[100000];
    
    strcpy(dateiname,argv[1]);
    
    FILE *inputdatei;
    inputdatei = fopen(dateiname,"r");
    
    if (inputdatei != NULL){
      fscanf (inputdatei, "%s", treeline);
//    fclose (inputdatei);   // if this is included there are strange errors in the parsenewick routine!!
    }
    else{ printf("ERROR: no such file found\n"); return -1;}


    int leaves;    
    leaves = atoi(argv[3]);

    strcpy(dateiname,argv[2]);
    
    FILE *inputdatei_gfs;
    inputdatei_gfs = fopen(dateiname,"r");
   
    float input_gfs[leaves];
    char gfsline[5000]; 


    
    if (inputdatei_gfs != NULL){
      int idex = 0;
//       printf("Gene frequency line is: ");
      while (fscanf (inputdatei_gfs, "%s", gfsline) != EOF){
//       printf ("%s\t", gfsline);
      input_gfs[idex] = atof(gfsline);
      if (idex >= leaves) {printf("\nWARNING: number of gene frequencies exceeds number of leaves - additional info is ignored\n");}
      idex++;
      }
//    fclose (inputdatei_gfs);   // if this is included there are strange errors in the parsenewick routine!!
    }
    else{
      printf("ERROR: no such file found\n");
      return -1;
    }
////////////////////*END reading the data */////////////////////////////////
  





///////////////////* several options are checked *//////////////////////////
  
  
//   -n  -> no neutraliy test is done
//   -n -s  -> samplingbias is tested -> only neutrality is not tested
//   -s  -> test for neutrality and for neutrality + sampling bias is done
//   -t theta -> no estimation is done
//   -r rho  -> no estimation is done
//   -c core -> no estimation is done
//   -d less details are printed 
//   -q quantity of runs for the tests is set to custom number
//   -a skip all but the results for a typical population
//   -p print tree
//   -b scaling of the tree is skipped
  

    // Set defaults for all parameters
    int notest_flag = 0, samplingbias_flag = 0, estimate_flag = 1, skipall_flag = 0, details_flag = 1, printtree_flag = 0, pansize_flag = 0, oldalgo_flag = 0, scale_flag = 1;
    float theta_input = 0.0, rho_input = 0.0, core_input = 0.0, millgenstoMRCA_input = 0.0;
    int runs_input = 0;
    
    
    if ( get_args(argc, argv, &notest_flag, &samplingbias_flag, &theta_input, &rho_input, &core_input, &estimate_flag, &details_flag, &skipall_flag, &runs_input, &printtree_flag, &pansize_flag, &millgenstoMRCA_input, &scale_flag ) == -1 ){
      printf("something was wrong with your options, maybe you wrote '-t200' instead of '-t 200' ?\n");
      return -1;
    }


  int runs = 10000;
  if (runs_input != 0) runs = runs_input;


if (details_flag == 1){
    if (argc == 4) printf("no options specified.\n");
    if (argc > 4){
      printf("options set:\n");
      if (skipall_flag == 0){
	if (notest_flag == 1) printf("Neutrality test will be skipped\n");
	if (notest_flag == 0) printf("Neutrality test will be done\n");
	if (samplingbias_flag == 1) printf("Test for sampling bias will be done\n");
      }else{
	printf("skipping all tests and estimations\n");
      }
      if (estimate_flag == 0)  printf("Parameters are not estimated. User sets parameters to:\n theta: %f\t rho: %f\t core: %f \n", theta_input, rho_input, core_input );  
      if(pansize_flag == 1) printf("The size of the pangenome will be computed\nUser sets Number of generations up to the MRCA to %f million generations.\nEffective population size will be estimated\n", millgenstoMRCA_input);
    }
    // print inputs:
    printf("\nNumber of leaves is %d\n", leaves);

    printf("\nGene frequency spectrum is given by:\n");
    printgfs(input_gfs,leaves,1,1); 

    printf("\nTree is:\n");
    printf ("%s\n", treeline);
    // end printing the inputs
}
//////////////* END reading parameters  *//////////////////////////////
  

  
  
  
  
  
/////////////////*parsing the newick tree*//////////////////////////////// 
Node * intree;

int n = 2*leaves-1; 
// intree = (Node *)malloc(n*sizeof(Node));
intree = new Node[n];
Node *listofleav[leaves];

int leaveornot = 0, place;
static int leaveid = 0;
static int innerid = 0;

char s1[100000];

strcpy(s1,treeline);
int last = strlen(s1)-1;

// start check root
if ( s1[last] == ';' ){
  s1[last] = '\0';
}
else{printf("ERROR: ';' in newick tree is missing\n"); return -1;}
// end check root


// this function parses the newick tree recursive (calls parsenewick(child1) and parsenewick(child2) )
parsenewick(s1,-1,n,intree,-1);
// if (printtree_flag == 1){
//  printf("\n--------tree-data---rooted----------------\n"); 
//  printrootedTree(intree+n-1);
//  printf("-------------------------------------------\n");
// }

findmaxdepth(intree,leaves);

// remove the root
if (details_flag == 1) printf("removing the root node:");
removeroot(intree+n-1);
if (details_flag == 1) printf("\tdone\n");

// root tree at first leave node 
if (details_flag == 1) printf("rooting tree at first leave node:");
rootTree(intree,NULL);
if (details_flag == 1) printf("\tdone\n");

// if (printtree_flag == 1){
//  printf("\n-----given--tree---------------------------\n");
//  printrootedTree(intree);
//  printf("-------------------------------------------\n");
// }

float depth;
depth = findmaxdepth(intree,leaves);

if (details_flag == 1) printf("Tree height is: %.2f\t", treeheight(intree,leaves));
  ////////////////////*END parsing the newick tree*//////////////////////////




  /////////* estimate the real tree height and scale the tree to that size*/////////////

float wantedheight;
wantedheight = expectedtreeheight(leaves);

float treeparts[leaves-1], estimatedheight, scalingfactor, neutraltreeheight;

if (scale_flag == 1){
  
  if (details_flag == 1){
    printf("Scale factor is thus: %.4f\n", (wantedheight*2.)/depth);
    printf("Scaling the tree to it's expected length...\n");
  }

  scalerootedtree(intree,(wantedheight*2.)/depth);

  if (details_flag == 1){
    printf("tree height is for the moment: %.2f\t which equals the expected height of a tree with %d leaves\n",treeheight(intree,leaves), leaves);
    printf("\n\n");
  }





  
  rootTree(intree,NULL);
  if (details_flag == 1) printf("Estimating the real tree height based on the newick file:\n");
  scalingfactor = estimate_treeheight_short(intree,leaves,treeparts,0);
  estimatedheight = scalingfactor*treeheight(intree,leaves);
  neutraltreeheight = estimatedheight;
  if (details_flag == 1){
    printf("Faktor is: %.2f\n", scalingfactor);
    printf("Estimated Height is: %.2f\n", estimatedheight);
  }
  rootTree(intree,NULL);
  scalerootedtree(intree,scalingfactor);
}
rootTree(intree,NULL);
if (details_flag == 1){
  printf("Tree height is now: %.2f\n",treeheight(intree,leaves));
  if (scale_flag == 0){
    printf("Tree has not been scaled\n");
  }
}

//////////////////////* END estimate the real tree height and scale the tree to that size////////////////////////////////

if (printtree_flag == 1){
 printf("\n-----scaled--tree---------------------------\n");
 printrootedTree(intree);
 printf("-------------------------------------------\n");
}




// // starting to test the symbolic computation
// 
// rootTree(intree,NULL);
// printf("rooted the tree...\n");
// unprob_symb(intree);
// printf("unprobed the tree symbolic\n");
// computeprobsall_symbolic(intree, leaves, x);
// 
// rootTree(intree,NULL);
// printf("rooted the tree...\n");
// unprob(intree);
// printf("unprobed the tree nonsymbolic\n");
// computeprobsall(intree, leaves, 0.5);
// 
// if (printtree_flag == 1){
//  printf("\n-----scaled--tree---------------------------\n");
//  rootTree(intree,NULL);
//  printrootedTree(intree);
//  printf("-------------------------------------------\n");
// }
// 
// 
// //   rootTree(intree,NULL);
// //   unprob_symb(intree);
// //   computeprobsall_symbolic(intree,leaves,x);
// //   rootTree(intree,NULL);
//    ex *theogfs_symb;
// // //   theogfs = (float *)malloc(anzahl*sizeof(float)); // this is old c code and does no longer work with c++ and ginac
//    theogfs_symb = new ex[leaves];
// //   treegfs_symbolic_slow(intree,leaves,theogfs_symb,x);
// // 
// // //printgfs_thin_symb(theogfs_symb,leaves);
// // 
//  int inc = 0;
// // for (inc = 0; inc < leaves; inc++){
// //  theogfs_symb[inc] =  theogfs_symb[inc].subs(x == 0.5).evalf();
// // }
// // printgfs_thin_symb(theogfs_symb,leaves);
// 
// 
// 
// //   rootTree(intree,NULL);
// //   unprob(intree);
// //   computeprobsall(intree,leaves,0.5);
// //   rootTree(intree,NULL);
//    float *theogfs_numeric;
// // //   theogfs = (float *)malloc(anzahl*sizeof(float));
//    theogfs_numeric = new float[leaves];
// //   treegfs(intree,leaves,theogfs_numeric,0.5);
// //   printgfs(theogfs_numeric,leaves,1.0,0.5);
// // 
// // 
// // printf("Done\n");
// 
// 
// 
// // try to compute the gfs in a clever way:
// rootTree(intree,NULL);
// initialize_tree(intree,leaves);
// check_probs(intree,leaves,x);
// 
// printf("_____________________________\n");
// comp_pkfhs(intree,leaves,x);
// check_probs(intree,leaves,x);
// 
//    ex *theogfs_symb_fast;
//    theogfs_symb_fast = new ex[leaves];
//   treegfs_symbolic_fast(intree,leaves,theogfs_symb_fast,x);
//   for (inc = 0; inc < leaves; inc++){
//     theogfs_symb_fast[inc] =  theogfs_symb_fast[inc].subs(x == 0.5).evalf();
//   }
//   printgfs_thin_symb(theogfs_symb_fast,leaves);
//   printf("_--------------------------------------\n");
//   printgfs_thin_symb(theogfs_symb,leaves);
//   printf("_--------------------------------------\n");
//    printgfs(theogfs_numeric,leaves,1.0,0.5);
//   
//   
//   for (inc = 0; inc < leaves; inc++){
//    theogfs_numeric[inc] = to_double(ex_to<numeric>(theogfs_symb_fast[inc].subs(x == 0.5).evalf()));
//   }
//   
//   printf("_--------------------------------------\n");
//   printgfs(theogfs_numeric,leaves,1.0,1.0);
// 
// 
// 
// //return 0;
// 
// // end to test the symbolic computation



  // initialize parameters
  float theta_hat , rho_hat, core_hat;
  float pvalue,pvalue_samplingbias;

  theta_hat = 0.;
  rho_hat = 0.;
  core_hat = 0.;
  
  
  //////////////////////* estimate theta and rho *///////////////////////////////////////////////
  Params *para;
  para = (struct Params *) malloc(sizeof(Params) + leaves * sizeof(double));
  para->anzahl = leaves;
  para->tree = intree;

  
  Params_symbolic *paraS;
  paraS = new Params_symbolic;
  paraS->anzahl = leaves;
  paraS->rhoS = x;
  

  // the numbers for the gfs from input
  int i;
  for(i=0; i<leaves; i++){
    para->datagfs[i] = input_gfs[i];
    paraS->datagfs[i] = input_gfs[i];
  }
  
    // -t or -r or -c is set to custom value
  if (estimate_flag == 0){
    printf("Setting theta ,rho and core:\n");
    theta_hat = theta_input;
    rho_hat = rho_input;
    core_hat = core_input; 
  }
  
  if (skipall_flag == 1) printf("theta = %f \t rho = %f \t core = %f\n" , theta_hat , rho_hat, core_hat); 
  
  
  float *estimated_gfs_theo, *est_gfs_giventree_theo;
  estimated_gfs_theo     = (float *)malloc(leaves*sizeof(float));
  est_gfs_giventree_theo = (float *)malloc(leaves*sizeof(float));
  

if (skipall_flag == 0) {


//   ex *theogfs_symb_fast;
//   theogfs_symb_fast = new ex[leaves];
    
    
  // only estimate if -t or -r or -c is not set to custom value
  if (estimate_flag == 1){
    printf("Estimating theta and rho...this may take some time\n");
    if (oldalgo_flag == 1){
      estimate(&theta_hat,&rho_hat,para);
    }
    else{
// // // // //         // this part works well for smaller trees
// //compute the symbolic formula in advance      
// // do this only once and give the functions to my_f_symbolic
//   cout << "Initializing the tree structure...\n";
//   rootTree(intree,NULL);
//   unprob_symb(intree);
//   initialize_tree(intree,leaves);
//   cout << "done.\nComputing probabilities along the tree...\n";
// //   cout << "We did get past the initialize tree\n";
// //   check_probs(tree,anzahl,para->rhoS);
//   comp_pkfhs(intree,leaves,x);
//   cout << "done.\nEstimation of theta and rho...\n";
// //    cout << "We computed the probabillties within the tree\n";
// //   check_probs(tree,anzahl,para->rhoS)
//   treegfs_symbolic_fast(intree,leaves,theogfs_symb_fast,x);
//   for(i=0; i<leaves; i++){
//     paraS->symbolicgfs[i] = theogfs_symb_fast[i];
// //     printf("\n\n\n\n");
// //     cout << theogfs_symb_fast[i];
// // // // //         
        
        
// // // // // //   this part works faster for larger trees
// //compute the pkfhs and the gfs along the tree during estimation
// this is faster than using the symbolic structure where many formulas appear many times
// my_f_numeric will need the tree to compute the pkfhs along the tree
  cout << "Initializing the tree structure...\n";
  rootTree(intree,NULL);
  unprob(intree);
  initialize_tree_numeric(intree,leaves);
  cout << "done.\nTesting to computing probabilities along the tree...\n";
//   cout << "We did get past the initialize tree\n";
//   check_probs(tree,anzahl,para->rhoS);
  comp_pkfhs_numeric(intree,leaves,2.8);
  cout << "done.\nEstimation of theta and rho...\n";
//    cout << "We computed the probabillties within the tree\n";
//   check_probs(tree,anzahl,para->rhoS)
//   treegfs_symbolic_fast(intree,leaves,theogfs_symb_fast,x);
//   for(i=0; i<leaves; i++){
//     paraS->symbolicgfs[i] = theogfs_symb_fast[i];
// //     printf("\n\n\n\n");
// //     cout << theogfs_symb_fast[i];
//   }
  // // // // //    
  
//  cout << "We successfully called treegfs_symbolic_fast\n";        
      estimate_numeric(&theta_hat,&rho_hat,para);
    }
  }
  //////////////////* END estimate theta and rho*//////////////////////////////////////////


  printf("parameters for neutral evolution:\t theta = %f \t rho = %f\n\n\n" , theta_hat , rho_hat ); 

// return 0;
  
  
  ///////////////* compute chisquare like statistic  *///////////////////////////

  float data_chi_square;
  

  // compute the expected GFS
  theoGfs_core(estimated_gfs_theo,leaves,theta_hat,rho_hat,0.);

//   /*compute probabilities for estimated theta and rho*/  //i can now do that much faster (see below)
//   rootTree(intree,NULL);
//   unprob(intree);
//   computeprobsall(intree,leaves,rho_hat);
//   rootTree(intree,NULL);
//   // compute the GFS given the tree
//   treegfs(intree,leaves,est_gfs_giventree_theo,rho_hat);  //this is not fast
  rootTree(intree,NULL);
//   treegfs_symbolic_fast(intree,leaves,theogfs_symb_fast,x); // this is fast but can only very slowly be evaluated
  unprob(intree);
  initialize_tree_numeric(intree,leaves);
  comp_pkfhs_numeric(intree,leaves,rho_hat);
  treegfs_numeric_fast(intree,leaves,est_gfs_giventree_theo,rho_hat);
//   int inc;
//   for (inc = 0; inc < leaves; inc++){
//      est_gfs_giventree_theo[inc] = to_double(ex_to<numeric>(theogfs_symb_fast[inc].subs(x == rho_hat).evalf()));
//   }
  
  

  if (details_flag == 1){
    printf("Estimated GFS:\n");
    printgfs(estimated_gfs_theo,leaves,1.,1.);
    printf("Input GFS:\n");
    printgfs(input_gfs,leaves,1.,1.);
    printf("Estimated GFS given the tree:\n");
//     printgfs(est_gfs_giventree_theo,leaves,theta_hat,rho_hat);
    printgfs(est_gfs_giventree_theo,leaves,theta_hat,rho_hat);
  }




  if (estimate_flag == 1){
    core_hat = input_gfs[leaves-1] - est_gfs_giventree_theo[leaves-1]*theta_hat/rho_hat;
    printf("estimated core size = %.2f\n",core_hat);
  }
  else{
    if (core_hat < input_gfs[leaves-1] - est_gfs_giventree_theo[leaves-1]*theta_hat/rho_hat - 5. ){
      if (details_flag == 1) printf("Note: core is set to %.0f a higher value would fit better.\n" , core_hat);
    }
    if (core_hat > input_gfs[leaves-1] - est_gfs_giventree_theo[leaves-1]*theta_hat/rho_hat + 5. ){
      if (details_flag == 1) printf("Note: core is set to %.0f a lower value would fit better.\n" , core_hat);
    }
  }



  data_chi_square = comp_chissquare(input_gfs, estimated_gfs_theo, leaves);
  if (details_flag == 1) printf("\n\nChi^2:\t%.2f\n\n", data_chi_square);



  // this statistic is a variant of the one above:
  // it takes into account how many genes in total appeared in the data and rescales the theo_gfs to this number
  // the other statistic, just compares the data gfs with the expected gfs for theta_hat and rho_hat

  // data_chi_square = comp_chissquare2(input_gfs, estimated_gfs_theo, leaves);
  // printf("alternative Chi^2:\t%.2f\n", data_chi_square);


  ///////////////*END compute chisquare like statistic  *///////////////////////////

  char outputname[] = "sims/chis.sim";


  if(notest_flag == 0){

    printf("Testing the hypothesis of neutral evolution:...this may take some time\n");
    if (details_flag == 1) printf("%d runs will be done:\n", runs);

    pvalue = compquality_without_estimating(leaves, theta_hat, rho_hat, runs,  outputname, data_chi_square,details_flag);

    if   (pvalue <0.05) printf("\n\nThe hypothesis of neutral evolution is not very probable. p-value is %.5f\n\n", pvalue);
    else printf("\n\nThe hypothesis of neutral evolution is not rejected as the p-value is %.3f\n\n", pvalue);

  }




}



//compute pairdiffs here in advance as intree might be scaled to another height later (during the sampling bias test result). result is printed at the end.
 float pairdiffs;
 pairdiffs = differentgenesinapair_giventhetree(intree, leaves, theta_hat, rho_hat)/2.0;











if (skipall_flag == 0){

  ////////////sampling bias test/////////////////////////////////////////77777
  if( samplingbias_flag == 1){

  float samplingtreeheight;
    
  if (details_flag == 1) printf("Re-estimating theta and rho for sampling bias\n");


    /* estimate the real tree height and scale the tree to that size*/
    if(scale_flag == 1){
      rootTree(intree,NULL);
      scalingfactor = estimate_treeheight_short(intree,leaves,treeparts,1);
      estimatedheight = scalingfactor*treeheight(intree,leaves);
      samplingtreeheight = estimatedheight;
      if (details_flag == 1){
	printf("Faktor is: %.2f\n", scalingfactor);
	printf("Estimated Height is: %.2f\n", estimatedheight);
      }
      rootTree(intree,NULL);
      scalerootedtree(intree,scalingfactor);
      rootTree(intree,NULL);
      if (details_flag == 1) printf("New tree height is now: %.2f\n",treeheight(intree,leaves));
    }

    
    /* estimate theta_samplingbias and rho_samplingbias */ 

    float theta_hat_samplingbias , rho_hat_samplingbias, core_hat_samplingbias;

    theta_hat_samplingbias = 0.;
    rho_hat_samplingbias = 0.;
    core_hat_samplingbias = 0.;


    // if notest_flag == 1 , estimate_flag == 1  and anyway samplingbias == 1   --> newly estimate parameters for sampling bias
    // 
    if (notest_flag == 1){
      if (estimate_flag == 1){
	printf("Estimating theta and rho under sampling bias...this may take some time\n");
	estimate(&theta_hat_samplingbias,&rho_hat_samplingbias,para);
      }
      if (estimate_flag == 0){
	printf("Setting theta ,rho and core for samplingbias to user input:\n");
	theta_hat_samplingbias = theta_input;
	rho_hat_samplingbias = rho_input;
	core_hat_samplingbias = core_input;
      }
    }



    // if notest_flag == 0  --> get the estimates for sampling bias from the estimates for neutrality
    if (notest_flag == 0){
      // rescale the parameters of theta_hat, rho_hat, core_hat to the sampling bias
      rho_hat_samplingbias = rho_hat*neutraltreeheight/samplingtreeheight;
      theta_hat_samplingbias = theta_hat*neutraltreeheight/samplingtreeheight;
    }


    printf("parameters under samplingbias:\ttheta = %f \t rho = %f\n\n\n" , theta_hat_samplingbias , rho_hat_samplingbias ); 


    // compute chisquare under samplingbias
    if (details_flag == 1) printf("computing chisquare under samplingbias\n");


    // compute the gfs with sampling bias
	  float *gfs_theo;
	  gfs_theo = (float *)malloc(leaves*sizeof(float));
	  float *fnotlost_theo;
	  fnotlost_theo = (float *)malloc(leaves*sizeof(float));
	  float *flost_theo;
	  flost_theo = (float *)malloc(leaves*sizeof(float));  
	  float *gfs_bias_theo;
	  gfs_bias_theo = (float *)malloc(leaves*sizeof(float));  
	  float *gfs_giventree_bias_theo;
	  gfs_giventree_bias_theo = (float *)malloc(leaves*sizeof(float)); 
      
	  // computing the function f_notlost
	  comp_fnotlost (leaves, fnotlost_theo, rho_hat_samplingbias, 10000 ); 
	  
	  // computing the function f_lost
	  comp_flost(leaves, flost_theo, fnotlost_theo);
	  
	  // compute the theoretical gfs without sampling bias for theta_samplingbais and rho_samplingbias)
	  theoGfs_core(gfs_theo, leaves, theta_hat_samplingbias, rho_hat_samplingbias, 0.);
	  
	  // compute the gfs with samplingbias based on the theoretical gfs without sampling bias from above
	  comp_gfs_bias(leaves, gfs_bias_theo, gfs_theo, flost_theo, fnotlost_theo, theta_hat_samplingbias , rho_hat_samplingbias, 0.);
	  
    // computed the gfs with sampling bias

    float data_chi_square_samplingbias;
    data_chi_square_samplingbias = comp_chissquare(input_gfs, gfs_bias_theo, leaves);
    if (details_flag == 1) printf("\nSampling bias:\tChi^2:\t%.2f\n\n", data_chi_square_samplingbias);

    /*compute probabilities for estimated/given theta and rho with samplingbias*/
    rootTree(intree,NULL);
    unprob(intree);
    computeprobsall(intree,leaves,rho_hat_samplingbias);
    rootTree(intree,NULL);

    // compute the GFS given the tree with samplingbias
    treegfs(intree,leaves,gfs_giventree_bias_theo,rho_hat_samplingbias);


    if (details_flag == 1){  
      printf("Estimated GFS with samplingbias:\n");
      printgfs(gfs_bias_theo, leaves,1.,1.);

      printf("Input GFS:\n");
      printgfs(input_gfs, leaves,1.,1.);

      printf("Estimated GFS with samplingbias given the tree:\n");
      printgfs(gfs_giventree_bias_theo,leaves,theta_hat_samplingbias,rho_hat_samplingbias);

    }
    // printf("This should be the same as:\n");
    // printgfs(est_gfs_giventree_theo,leaves,theta_hat,rho_hat);

    
    if (estimate_flag == 1){
      core_hat_samplingbias = input_gfs[leaves-1] - gfs_giventree_bias_theo[leaves-1]*theta_hat_samplingbias/rho_hat_samplingbias;
      printf("estimated core size under sampling bias = %.2f\n",core_hat_samplingbias);
    }
    else{
      if (core_hat_samplingbias < input_gfs[leaves-1] - gfs_giventree_bias_theo[leaves-1]*theta_hat_samplingbias/rho_hat_samplingbias - 5.){
	if (details_flag == 1) printf("Note: sampling core is set to %.0f a higher value would fit better.\n" , core_hat_samplingbias);
      }
      if (core_hat_samplingbias > input_gfs[leaves-1] - gfs_giventree_bias_theo[leaves-1]*theta_hat_samplingbias/rho_hat_samplingbias + 5.){
	if (details_flag == 1) printf("Note: sampling core is set to %.0f a lower value would fit better.\n" , core_hat_samplingbias);
      }
    }

  // Testing for neutral evolution + sampling bias
    printf("Testing the hypothesis of sampling bias:...this may take some time\n");
    if (details_flag == 1) printf("%d runs will be done:\n", runs);
    
    char samplingoutputname[] = "sims/chis_samplingbias.sim";
    pvalue_samplingbias = compquality_without_estimating_samplingbias(leaves, theta_hat_samplingbias, rho_hat_samplingbias, runs,  samplingoutputname, data_chi_square_samplingbias, details_flag);

    if   (pvalue_samplingbias <0.05) printf("\n\nThe hypothesis of neutral evolution + sampling bias is not very probable. p-value is %.5f\n\n", pvalue_samplingbias);
    else printf("\n\nThe hypothesis of neutral evolution + sampling bias is not rejected as the p-value is %.3f\n\n", pvalue_samplingbias);
  }

  ///////////////// END Sampling Bias Test ////////////////////////////////////////////
}





// printing some characteristics


printf("\n\n-----------------------------------------------------------------------\n");
printf("Some characteristics of a typical population with parameters:\n");
printf("theta = %.2f \t rho = %.2f \t core = %.2f\n\n\n", theta_hat, rho_hat, core_hat);

printf("The average number of genes per individual:\t\t\t\t\t%.2f\n", theta_hat/rho_hat + core_hat);
printf("For each fixed tree this number is still the same: \t\t\t\t%.2f\n\n", theta_hat/rho_hat + core_hat);

// printf("Expected number of different genes in a pair: %.2f \n", theta_hat/(rho_hat + 1.) );
// printf("Expected number of different genes in a pair taken from the tree: %.2f \n", pairdiffs);

printf("The average number of genes in two individuals:\t\t\t\t\t%.2f\n", theta_hat/rho_hat + theta_hat/(rho_hat+1.) + core_hat);  //== theo_curr_pansize_sum(theta_hat, rho_hat, 2 , core_hat)
printf("Expected number of genes within two individuals taken from the given tree:\t\t%.2f\n\n", theta_hat/rho_hat + 0.5*pairdiffs + core_hat);

printf("The average number of genes in %d individuals:\t\t\t\t\t%.2f\n",leaves, theo_curr_pansize_sum(theta_hat, rho_hat, leaves , core_hat) );

float fixedtreepangenomesize = core_hat;
for(i=0; i<leaves; i++){
fixedtreepangenomesize += est_gfs_giventree_theo[i]*theta_hat/rho_hat;
}
printf("The average number of genes in %d individuals for the given tree:\t\t%.2f\n\n", leaves, fixedtreepangenomesize );


printf("The average number of genes in 1000 individuals:\t\t\t\t%.2f\n\n", theo_curr_pansize_sum(theta_hat, rho_hat, 1000 , core_hat) );

float pansize = theo_curr_pansize_sum(theta_hat, rho_hat, 10000 , core_hat);
printf("The average number of genes in 10000 individuals:\t\t\t\t%.2f\n", pansize );
//  theo_curr_antisuprasize(float theta, float rho, uint64_t n, uint64_t border)
float antisuprasize =  theo_curr_antisuprasize(theta_hat,rho_hat,10000, 100);
printf("The average number of genes present in more than 100 of 10000 individuals:\t%.2f\n\n", pansize - antisuprasize);


printf("The average number of new genes found in the %d'ths individual:\t\t\t%.2f\n\n", leaves+1, theta_hat/(leaves+rho_hat) );

float mysuprasize = theo_suprasize(theta_hat,rho_hat,0.01)+core_hat;
printf("The persistant pangenome size (average number of genes present in more than 1 percent of the population) is:\t%f\n\n", mysuprasize);


printf("If you used panicmage please cite: Baumdicker, F., Hess, W. R., & Pfaffelhuber, P. (2012). The infinitely many genes model for the distributed genome of bacteria. Genome Biology and Evolution, 4(4), 443–456.");


// compute the size of the pangenome
if(pansize_flag == 1){

uint64_t effective_popsize;
if (millgenstoMRCA_input < 1.0){ 
  effective_popsize =  ( ( millgenstoMRCA_input*1000000.)/ neutraltreeheight );
}
else{
effective_popsize =  (uint64_t) ( ( (uint64_t) millgenstoMRCA_input*1000000.)/ neutraltreeheight );
}
printf("The effective population size is given by %" PRIu64 "\n", effective_popsize);
pansize = theo_curr_pansize_sum(theta_hat, rho_hat, effective_popsize , core_hat);
printf("The average number of genes in the population (pansize) is:\t%.2f\n", pansize );
// antisuprasize = theo_curr_antisuprasize(theta_hat,rho_hat,effective_popsize, (uint64_t) (effective_popsize/100) );
// printf("The average number of genes present in more than 1 percent of the population is:\t%f\n\n", pansize - antisuprasize);
float theta_gen, rho_gen;
int theta_power = 0, rho_power = 0;
theta_gen = theta_hat/(2.*effective_popsize);
while (theta_gen < 1.){
  theta_gen = theta_gen * 10.;
  theta_power++;
}
rho_gen = rho_hat/(2.*effective_popsize);
while (rho_gen < 1.){
  rho_gen = rho_gen * 10.;
  rho_power++;
}
printf("per generation rate of gene gain:\t\t %f * 10^-%d\n", theta_gen, theta_power );
printf("per generation rate of gene loss for each gene:\t %f * 10^-%d\n", rho_gen, rho_power );
}
else{
printf("NOTE: no number of generations up to the most recent commen ancestor (MRCA) given.\nPansize and per generation rates are only computable if the paramter -g is given\n");
printf("To compute the pansize and per generation rates without reestimating parameters use the options:\n -a -t %.5f -r %.5f -c %.2f -g [float]\n where [float] is the number of million generations up to the MRCA \n\n", theta_hat, rho_hat, core_hat);
}



// printing some results to file

FILE *RESULT_output;
RESULT_output = fopen("panicmage_estimatedparameters.txt","w");
fprintf(RESULT_output, "%.10f\t%.10f\t%.2f\n", theta_hat, rho_hat, core_hat);
fclose(RESULT_output);

if (notest_flag == 0){
  RESULT_output = fopen("panicmage_pvalue.txt","w");
  fprintf(RESULT_output, "%.2f\n", pvalue);
  fclose(RESULT_output);
}


  
return 0;
  
  
  
  
  
  
  
  
  
  
  
  
  
}