/***** treestatistics.h ***********************************************
 * Description: statistics and estimation procedure
 * Author: Franz Baumdicker, baumdicker@stochastik.uni-freiburg.de ; Peter Pfaffelhuber, pp@stochastik.uni-freiburg.de
 * File created 2010.
 *****************************************************************/

#include "treestructure.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



extern const gsl_rng_type * T;
extern gsl_rng * r;

void comp_fnotlost (int leaves, float *fnotlost_k, float rho, unsigned long int largeN );
void comp_flost (int leaves, float *flost_k,float *fnotlost_k );
void comp_gfs_bias (int leaves, float *gfs_bias_theo, float *gfs_theo, float *flost_k, float *fnotlost_k, float theta, float rho, float core);

void initialize_tree(Node *node, int leaves);
void comp_pkfhs(Node *node, int leaves, const ex & rhoS);
void check_probs(Node *node, int leaves, const ex & rhoS);
void treegfs_symbolic_fast(Node *tree, int leaves, ex *gfs_k_symb, const ex & rhoS);

void initialize_tree_numeric(Node *node, int leaves);
void comp_pkfhs_numeric(Node *node, int leaves, float myrho);
void treegfs_numeric_fast(Node *tree, int leaves, float *gfs_k_numeric, float myrho);






// distance function according to least square
double lsqu(float *gfs, float *theogfs, int anzahl,float t,float r){
  double sum = 0.0;
  int i;
  for(i=0;i<anzahl;++i)
    sum += std::pow((gfs[i]-theogfs[i]*t/r),2);
  return sum;
}

// distance function according to maximum likelihood
double llh(float *gfs, float *theogfs, int anzahl, float t, float r){
  double sum = 0.0;
  int i;
  for(i=0; i<anzahl; i++)
    sum += theogfs[i]*t/r - gfs[i] * std::log(theogfs[i]*t/r);
  return sum;
}



// theoretical means for the gene frequency spectrum
void theoGfs_core (float *gfs, int nsam, float theta, float rho, float core) {
  int k,l;
  float l2;
  for (k = 1; k <= nsam; k++){
    gfs[k-1] = theta/((float) k );
    for( l = nsam; l >= nsam-k+1; l--){
     gfs[k-1] *= l;
    }
    for( l2 = ((float)nsam)-1.+rho ; l2 > ((float)nsam)-k+rho-0.5; l2 -= 1.){
      gfs[k-1] *= 1./l2;
    }
  }
  gfs[nsam-1] = gfs[nsam-1]+core;  
  // R code to do the same
  //   for(k in 1:nsam) {
  //     gfs[k]<-theta/k * prod(nsam:(nsam-k+1))/(prod((nsam-1+rho):(nsam-k+rho)))
  //   }
  //     gfs[nsam] <- gfs[nsam]+core
}





// compute the expected treeheight (tree is initially scaled to this size
float expectedtreeheight(int leaves){
  float result = 0.;
  int i;
  for (i = 2; i <= leaves; i++){
    result += 2./( ((float)i)* (((float)i) -1.) );
  } 
  return result;
}


// compute the expected number of diffent genes in a randomly choosen pair for a given tree
float differentgenesinapair_giventhetree(Node *tree, int leaves, float theta, float rho){
  int i,j, count;
  float sum = 0.0;
  count = 0;
  for(i=0; i < leaves-1; i++){
    rootTree(tree+i,NULL);
//     printf("ID: %d\n", tree[i].id);
    for(j=i+1; j < leaves; j++){
      count++;
      sum += 2.*theta/rho*( 1.-exp(-rho/2.*(tree[j].disttoroot)) );
    }
  }
//   printf("Number of Pairs: %d\n", count);
  return sum/( (float) count);
}




/*------------function which is minimized-----------------------------------*/
double my_f (const gsl_vector *v, void *params)
{
  extern int g_includecoreflag;
  cout << "This is the OLD function. I should not be here... STOPPING:\n";
  return 0;
  double rho, theta, sum = 0.;
  Params *para = (Params *)params;
       
  int anzahl = para->anzahl;
  Node *tree = para->tree;
  
  rho = gsl_vector_get(v, 0); 
  theta = gsl_vector_get(v, 1);
  
  rootTree(tree,NULL);
  unprob(tree);
  computeprobsall(tree,anzahl,rho);
  rootTree(tree,NULL);
  float *theogfs;
  theogfs = (float *)malloc(anzahl*sizeof(float));
  treegfs(tree,anzahl,theogfs,rho);
  
  float *datagfs;
  datagfs = (float *)malloc(anzahl*sizeof(float));
  
  int i;
  for (i = 0; i < anzahl; i++) {
    datagfs[i] = para->datagfs[i];
  }
 
  //return lsqu(datagfs,theogfs,anzahl-1,theta,rho);
  //return llh(datagfs, theogfs, anzahl-1,theta,rho);
  if(g_includecoreflag == 1){
    return llh(datagfs, theogfs, anzahl,theta,rho);
  }else{
    return llh(datagfs, theogfs, anzahl-1,theta,rho);
  }
}







/*estimate the parameters theta and rho*/

// Params are:
// Params *para;
// para = malloc(sizeof(Params) + number * sizeof(double));
// para->anzahl = number;
// para->tree = tree;
// 
// /* put the gfs data in the array */
// for(i=0; i<leaves; i++){
//   para->datagfs[i] = randoms[i];
// }

void estimate(float *theta_hat, float *rho_hat, Params *para){   
           
       const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
       gsl_multimin_fminimizer *s = NULL;
       gsl_vector *ss, *x;
       gsl_multimin_function minex_func;
     
       size_t iter = 0;
       int status;
       double size;
       
       /*optional set the starting points to a first simple estimate*/
       
       int i, n;
       float c1, c2, firsttheta, firstrho;
       c1 = 0.;
       c2 = 0.;
       n = para->anzahl;
       for (i = 0; i < n; i++){
       c1 += (i+1.)/n * para->datagfs[i];
       c2 += (i+1.)*(n-i-1.)/(n*(n-1.)) * para->datagfs[i];
       }
       firstrho = c2/(c1-c2);
       firsttheta = c2 + firstrho*c2;
       if (firstrho > 20.) firstrho = 20.;
       if (firstrho < 0.1) firstrho = 0.1;
       if (firsttheta > 10000.) firsttheta = 10000.;
       if (firsttheta < 50.) firsttheta = 50;
       
       
       /* Starting point */
       x = gsl_vector_alloc (2);
       
        gsl_vector_set (x, 0, 1.);  // fixed starting point rho  
        gsl_vector_set (x, 1, 1000.0); // fixed starting point theta
       
//        gsl_vector_set (x, 0, firstrho);  // estimated starting point rho
//        gsl_vector_set (x, 1, firsttheta); // estimated starting point theta

       //        printf("firsttheta = %.0f\tfirstrho = %.2f\n",firsttheta, firstrho);
       
     
       /* Set initial step sizes to 1 */
       ss = gsl_vector_alloc (2);
       gsl_vector_set_all (ss, 1.0);
     
       /* Initialize method and iterate */
       minex_func.n = 2;
       minex_func.f = my_f;
       minex_func.params = para;
     
       s = gsl_multimin_fminimizer_alloc (T, 2);
       gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
     
       do
         {
           iter++;
           status = gsl_multimin_fminimizer_iterate(s);
           
           if (status) 
             break;
     
           size = gsl_multimin_fminimizer_size (s);
           status = gsl_multimin_test_size (size, 1e-1);
     
           if (status == GSL_SUCCESS)
             {
//                printf ("converged to minimum at\n");
             }
     
//         printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
//                    iter,
//                    gsl_vector_get (s->x, 0), 
//                    gsl_vector_get (s->x, 1), 
//                    s->fval, size);
         }
       while(status == GSL_CONTINUE && iter < 100);
       
//        printf("\n");
       if (iter != 100) {
	 theta_hat[0] = (float) gsl_vector_get (s->x, 1);
	 rho_hat[0]   = (float) gsl_vector_get (s->x, 0);
       }
       else {
	 printf("WARNING: minimizer did not converge: \t size: %.2f \t\t theta: %.2f \t\t rho: %.2f\n", size, (float) gsl_vector_get (s->x, 1),  (float) gsl_vector_get (s->x, 0) );
// 	 theta_hat[0] = 0;
// 	 rho_hat[0] = 0;
	 theta_hat[0] = (float) gsl_vector_get (s->x, 1);
	 rho_hat[0]   = (float) gsl_vector_get (s->x, 0);
       }
       gsl_vector_free(x);
       gsl_vector_free(ss);
       gsl_multimin_fminimizer_free (s);
}






/* compute the chi square like statistic for the gene frequency spectrum*/
float comp_chissquare(float *datagfs, float *gfs_theo, int leaves){
  float chi_square;
  int jdex;
  chi_square = 0.;
  for(jdex = 0; jdex<leaves-1; jdex++){
    if (gfs_theo[jdex] == 0.) printf("WARNING: zero in denominator at j = %d \n", jdex);
    chi_square += std::pow( (datagfs[jdex] - gfs_theo[jdex]) , 2.)  / (gfs_theo[jdex]);
  }
  return chi_square;
}



/* compute the chi square like statistic for the gene frequency spectrum taking the number of genes in the dataset into account*/ // not used in IMaGe
float comp_chissquare2(float *datagfs, float *gfs_theo, int leaves){
  float sumpies = 0.;
  float *pies;
  pies = (float *)malloc((leaves-1)*sizeof(float));
  int j;
  int numgenes = 0;
  for(j=0; j<leaves-1; j++){
    numgenes += ((int) datagfs[j]);
  }
  if (numgenes == 0) printf("WARNING: number of genes is zero");
  for(j = 0; j<leaves-1; j++){
      sumpies += gfs_theo[j];
  }
  for(j = 0; j<leaves-1; j++){
      pies[j] = ((float) numgenes)*gfs_theo[j]/sumpies;
  }
  return comp_chissquare(datagfs, pies, leaves);
}




//generate genes for a given tree and given parameters
void creategenenumbers(int leaves, float *gfs_mean , float *randoms, float theta, float rho){
  int i;
//   const gsl_rng_type * T;  global variable is used
//   gsl_rng * r;  global variable is used
  /* create a generator chosen by the 
          environment variable GSL_RNG_TYPE */
  for(i = 0; i<leaves; i++){
    randoms[i] = (float) gsl_ran_poisson(r, gfs_mean[i]*theta/rho );
  }
//   gsl_rng_free (r);
}




void simonetree(float *gfs_sim,int leaves,float rho_hat){
  
	int n;
	n = 2*leaves-1; 
        symbol localx("localrho");
        Node *list[leaves];
	Node * tree;
        tree = new Node[n];

     maketree(tree,list,leaves);   
     
     
  rootTree(tree,NULL);
  initialize_tree(tree,leaves);
//   cout << "initialized ok\n";
//   check_probs(tree,anzahl,para->rhoS);
  comp_pkfhs(tree,leaves,localx);
//    cout << "We computed the probabillties within the tree\n";
//   check_probs(tree,anzahl,para->rhoS)
  ex *theogfs_fast_local;
  theogfs_fast_local = new ex[leaves];
//   cout << "so far ok\n";
  treegfs_symbolic_fast(tree,leaves,theogfs_fast_local,localx);
//    cout << "treegfs_symbolic_fast is ok\n";
  int i;
  for (i = 0; i < leaves; i++){
     gfs_sim[i] = to_double(ex_to<numeric>(theogfs_fast_local[i].subs(localx == rho_hat).evalf()));
  }
  
  delete [] tree;
  delete [] theogfs_fast_local;
  return;
}






void simonetree_numeric(float *gfs_sim,int leaves,float rho_hat){ 
	int n;
	n = 2*leaves-1; 
        Node *list[leaves];
	Node * tree;
        tree = new Node[n]; 

     maketree(tree,list,leaves);   

     
     
  rootTree(tree,NULL);
//   cout << "rooted ok\n";
  initialize_tree_numeric(tree,leaves);
//   cout << "We did get past the initialize tree\n";
  comp_pkfhs_numeric(tree,leaves,rho_hat);
//    cout << "We computed the probabillties within the tree\n";
  float *theogfs_fast_local;
  theogfs_fast_local = new float[leaves];
  treegfs_numeric_fast(tree,leaves,gfs_sim,rho_hat);
  
  delete [] tree;
  delete [] theogfs_fast_local;
  return;
}



float compquality_without_estimating(int leaves, float true_theta, float true_rho, int runs, char *filename, float data_chisquare, int details){

  //  printf("Compute Quality of estimators:\n");
  
  float pvalue = 0.;

  float *gfs_sim, *randoms;
  gfs_sim = new float[leaves];
  randoms = new float[leaves];

  float *gfs_theo, *pies;
  float sumpies,chi_square,chi_square2;
 
  gfs_theo = new float[leaves];
  pies = new float[leaves];

  FILE *CHI_output;
  CHI_output = fopen( filename, "w");

  // FILE *CHI_output2;
  // CHI_output2 = fopen( "empirical2.dat", "w");

  FILE *GFS_output;
  GFS_output = fopen("sims/GFS.sim","w");

  int steps = ( (int)  (((float) runs )/100.)  );
  if (runs < 100) steps = 10;
  if (runs < 10) steps = 1;
  if (details == 0) steps = runs + 1;
  int idex;
  
  printf("0%%|");
  for(idex = 0; idex< 100; idex++){         
     printf(".");
  }
  printf("|100%%\n0%%|");

  float explength = expectedtreeheight(leaves), reallength;
  //  printf("expected length of tree for %d individuals is: \t%f\n", leaves, explength);

  int j;
//   Node *list[leaves];

//   Params *para;
//   para = (struct Params * )malloc(sizeof(Params) + leaves * sizeof(double));

  float theta_hat, rho_hat;

      Node *tree;
//       tree = 0;
//       tree = new Node[n];
  int i;
  
  // simulate the distribution of the chi square like statistic
  for(idex = 0; idex< runs; idex++){         

//     simonetree(gfs_sim,leaves,true_rho);   // this works as well but is probably slower than the numeric computation
    simonetree_numeric(gfs_sim,leaves,true_rho);    
      
     
//      cout << "so far okidoki\n";

     /* create random gfs for this tree*/
     creategenenumbers(leaves,gfs_sim,randoms,true_theta,true_rho);
     
      //      ESTIMATION IS NOT DONE
      //    (this is the old and thus slow algo, be carefull if reactivated)
      // 	this part estimates new parameters for each random tree
      //      /* estimate theta and rho based on the current tree*/
      // 	free(para);
      //      
      // 	/* scale tree as this info is never known */
      // 	rootTree(tree,NULL);
      // 	reallength = findmaxdepth(tree,leaves)/2.;
      // // 	printf("simulated height: \t%f\n", reallength);
      // // 	printf("Factor: \t%f\n", explength/reallength);
      // 	scalerootedtree(tree,explength/reallength);
      // // 	printf("New tree height is now: %.2f\n", treeheight(tree,leaves));
      // 	rootTree(tree,NULL);
      // 	//printf("scaling done\n");
      // 
      // 	// Params *para;
      // 	para = malloc(sizeof(Params) + leaves * sizeof(double));
      // 	para->anzahl = leaves;
      // 	para->tree = tree;
      // 	/* put the random numbers from above in the array */
      // 	for(j=0; j<leaves; j++){
      // 	  para->datagfs[j] = randoms[j];
      // 	}
      // 	
      //      estimate(&theta_hat,&rho_hat,para);
      //      //printf("estimating done\n");
      //      //printf("theta = %.2f \t rho = %.2f\n" , theta_hat , rho_hat );

      // 	WE SET THETA AND RHO TO THE TRUE PARAMETERS
      theta_hat = true_theta;
      rho_hat = true_rho;
     
      
     //printf("computing the GFS with these parameters:\n");
     theoGfs_core(gfs_theo,leaves,theta_hat,rho_hat,0);
     //printf("\n");
     //printgfs(gfs_theo,leaves,1,1);
     //printf("\n");
     
             
      //      chi_square2  = comp_chissquare2(randoms, gfs_theo, leaves);
     chi_square = comp_chissquare(randoms, gfs_theo, leaves);
      //      printf("%.2f \tvs.\t %.2f\n", chi_square, chi_square2 );
     
     
     fprintf(CHI_output, "%.2f\t", chi_square);
     fflush(CHI_output);

      //      fprintf(CHI_output2, "%.2f\t", chi_square2);
      //      fflush(CHI_output2);
     
     // print the simulated GFS in GFS.sims
     for(j=0; j<leaves-1; j++){
        fprintf(GFS_output, "%.1f\t" ,randoms[j]);
     }
     fprintf(GFS_output, "%.1f\n" , randoms[leaves-1]);
     fflush(GFS_output);
     
     if (chi_square > data_chisquare) pvalue += 1.;
     
     
//      if (!((idex+1) % steps) )  printf("%d\n", idex+1);
     if (!((idex+1) % steps) ){
        printf(".");
        fflush(stdout);
     }
        
//      rootTree(tree,NULL);
//      emptytree(tree);
  }

  printf("|100%%\n");
  fflush(stdout);

  fclose(CHI_output);
  // fclose(CHI_output2);
  fclose(GFS_output);

  pvalue = pvalue/((float)runs);
  return pvalue;
}














///////////////////////same with sampling bias
float compquality_without_estimating_samplingbias(int leaves, float true_theta, float true_rho, int runs, char *filename, float data_chisquare, int details){

  //  printf("Compute Quality of estimators:\n");

  // printf("WARNING: note that one should estimate other parameters under sampling bias\n"); this is now done by IMaGe

  float pvalue = 0.;

  int n;
  n = 2*leaves-1; 

  Node *tree;
  tree = (Node *)malloc(n*sizeof(Node));

  float *gfs_sim, *randoms;
  gfs_sim = (float *)malloc(leaves*sizeof(float));
  randoms = (float *)malloc(leaves*sizeof(float));  

  float *gfs_theo, *pies;
  float sumpies,chi_square, chi_square2;

  gfs_theo = (float *)malloc(leaves*sizeof(float));  
  pies = (float *)malloc((leaves-1)*sizeof(float));  
  
  FILE *CHI_output;
  CHI_output = fopen( filename, "w");

  //   FILE *CHI_output2;
  //   CHI_output2 = fopen( "empirical2.dat", "w");

  FILE *GFS_output;
  GFS_output = fopen("sims/GFS_samplingbias.sim","w");

  int steps = ( (int)  (((float) runs )/100.)  );
  if (runs < 100) steps = 10;
  if (runs < 10) steps = 1;
  if (details == 0) steps = runs + 1;
  int idex;

  float explength = expectedtreeheight(leaves), reallength;
  //  printf("expected length of tree for %d individuals is: \t%f\n", leaves, explength);

  int j;
  Node *list[leaves];

  Params *para;
  para = (struct Params * )malloc(sizeof(Params) + leaves * sizeof(double));

  float theta_hat, rho_hat;

  theta_hat = true_theta;
  rho_hat = true_rho;

  // compute the gfs with sampling bias     
//     if (details == 1) printf("computing the GFS with the given parameters:\n");
    float *fnotlost_theo;
    fnotlost_theo = (float *)malloc(leaves*sizeof(float));

    float *flost_theo;
    flost_theo = (float *)malloc(leaves*sizeof(float));  
    
    float *gfs_bias_theo;
    gfs_bias_theo = (float *)malloc(leaves*sizeof(float));  

    //comp_fnotlost (int leaves, float *fnotlost_k, float rho, int smalln, unsigned long int largeN )
    comp_fnotlost (leaves, fnotlost_theo, rho_hat, 10000 ); 
    comp_flost(leaves, flost_theo, fnotlost_theo);
    //theoGfs_core(float *gfs, int nsam, float theta, float rho, float core)
    theoGfs_core(gfs_theo, leaves, theta_hat, rho_hat, 0.);

    //comp_gfs_bias (int leaves, float *gfs_bias_theo, float *gfs_theo, float *flost_k, float *fnotlost_k, float theta, float rho, float core)
    comp_gfs_bias(leaves, gfs_bias_theo, gfs_theo, flost_theo, fnotlost_theo, theta_hat , rho_hat, 0.);
//     if (details == 1) printf("done\n");
  // computed the gfs with sampling bias


  // simulate the distribution of the chi square like statistic with sampling bias
  for(idex = 0; idex< runs; idex++){
      /* create a random coalescent tree*/
      free(tree);
      tree = (Node *)malloc(n*sizeof(Node));
      maketree(tree,list,leaves);
      
      enlargebranches(tree,leaves, generate_random_enlargementsize(leaves, 10000) );
      
      /* compute parameters of the poisson distribution*/
//       TODO this should as well be replaced by the fast algorithm
      rootTree(tree,NULL);
      unprob(tree);
      computeprobsall(tree,leaves,true_rho);
      rootTree(tree,NULL);
      treegfs(tree,leaves,gfs_sim,true_rho);

      /* create random gfs for this tree*/
      creategenenumbers(leaves,gfs_sim,randoms,true_theta,true_rho);

  //  printgfs(randoms,leaves,1.,1.);
      
      
  //    ESTIMATION IS NOT DONE
  // 	THEREFORE WE SET THETA AND RHO TO THE TRUE PARAMETERS
  // 	theta_hat = true_theta;  already done before the for loop
  // 	rho_hat = true_rho;
  
  //       printf("%.2f\n", chi_square);


      chi_square  = comp_chissquare(randoms, gfs_bias_theo, leaves);
      chi_square2 = comp_chissquare2(randoms, gfs_bias_theo, leaves); 
  //  printf("Vergleich:\t%.2f \tvs.\t %.2f\n", chi_square, chi_square2 );
      
      
      fprintf(CHI_output, "%.2f\t", chi_square);
      fflush(CHI_output);
      
      // print the simulated GFS in filename
      for(j=0; j<leaves-1; j++){
	  fprintf(GFS_output, "%.1f\t" ,randoms[j]);
      }
      fprintf(GFS_output, "%.1f\n" , randoms[leaves-1]);
      fflush(GFS_output);
      
      if (chi_square > data_chisquare) pvalue += 1.;
      
      if (!(idex % steps) )  printf("%d\n", idex);
  }

  fclose(CHI_output);
  //   fclose(CHI_output2);
  fclose(GFS_output);  

  pvalue = pvalue/((float)runs);
  return pvalue;
}
























///////////////// estimate treeheigth ////////////////////////////////////////
float disttoleafs(Node *tree, Node *innernode, int n){
  int j;
  float min = 1000.;
  for(j=0; j<n; j++){
    rootTree(tree+j,NULL);
    if (innernode->disttoroot < min){
      min = innernode->disttoroot;
    }
  }
  return min;
}



int compare(const void *i1, const void *i2) {
  return *(int *)i1 - *(int *)i2; /* sort from low to high */
}






// estimates the height of a tree for given (unscaled) lengths 

float estimate_treeheight_short(Node *tree, int n, float *treeparts, int samplingbias){

  //get the tree lengths from tree
  rootTree(tree,NULL);

  /* compute treeparts */
    int j;
    for (j=0; j<n-2; j++){
      treeparts[j] = disttoleafs(tree,tree+n+j,n);
    }
    treeparts[n-2] = findmaxdepth(tree, n)/2.;

    // sort treeparts
    qsort(treeparts, n-1, sizeof(float), compare);

    // take differences
    for (j=n-2; j>0; j--){
      treeparts[j] = treeparts[j]-treeparts[j-1];
    }
  /* treeparts computed */

  // // define treeparts for n = 5  // just for testing
  // treeparts[0] = 0.1*log(2); // *log(2) is the median *1 is mean
  // treeparts[1] = 0.16666666666*log(2);
  // treeparts[2] = 0.33333333333*log(2);
  // treeparts[3] = 1.*log(2);

  float sum1 = 0, height=0.;
  float t, lambda;

  for (j=0; j<n-1;j++){
  t = n-j;
  lambda = (t*(t-1.))/2.;  // lambda is \binom{t}{2}
  sum1 +=  lambda * treeparts[j] ;
  }
  height = 1./( sum1/(n-1.) ) ;  // this is the scaling factor


  // skip the last length T_leaves (j == 0) if we test for sampling bias
  float height2;
  sum1 = 0;
  for (j=1; j<n-1;j++){
  t = n-j;
  lambda = (t*(t-1.))/2.;  // lambda is \binom{t}{2}
  sum1 +=  lambda * treeparts[j] ;
  }

  // compute the last part of the sum which is different as this is now the enlarged branch
  // for this estimator T_leaves is assumed to be exponentially distributed
  sum1 += 1./( expectedtreeheight(10000)-expectedtreeheight(n-1) ) * treeparts[0] ;
  height2 = 1./( sum1/(n-1.) ) ;  // this is the scaling factor with sampling bias

  if (samplingbias == 0) return  height;
  if (samplingbias == 1) return  height2;
}











