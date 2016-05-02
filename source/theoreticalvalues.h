/***** theoreticalvalues.h ***********************************************
 * Description: computation of theoretical values
 * Author: Franz Baumdicker, baumdicker@stochastik.uni-freiburg.de ;  Peter Pfaffelhuber, pp@stochastik.uni-freiburg.de
 * File created 2010.
 *****************************************************************/


#include<stdint.h>

// average genome size
float theo_averGenes(float theta, float rho, float core){
return theta/rho + core;
}

// pairwise differences
float theo_aver_pairdiff(float theta, float rho){
return theta/(1.+rho) ;
}


// gene frequency spectrum
// defined in treestatistics.h
// void theoGfs_core (float *gfs, int nsam, float theta, float rho, float core)


// number of k clades
// just done in R so far





///////////size of the pangenome computed by sum --> faster than summing the gfs
float theo_curr_pansize_sum(float theta, long double rho,uint64_t n, float core){
  long double sum = 0;
  uint64_t k;

  //   this is for printing the progressbar
  uint64_t total = n;
  uint64_t steps = ((uint64_t) (total/50));
  if (n < 100) {total = n; steps = n;}
  if (n > 10000){
  printf("Computing the Pangenome size...\tthis may take some time\n");
  printf("0%%");
  int zaehl;
  for (zaehl = 1; zaehl <= 50; zaehl++) printf(".");
  printf("100%%\n0%%");
  }
  
  // calculation
  for (k=0; k < n; k++){
    sum = sum + (long double) (1./((long double) ((rho)+k) ));
    if (!(k % steps) )  {
      if (n > 10000){
	printf("."); fflush(stdout);
      }
    }
  }
  if (n > 10000) {printf("100%%\n");}
  return theta*sum+core;
}





// size of the antisupragenome
// this computes the number of genes in frequency smaller than border/n  
// this is faster than computing the genes above frequency border/n
// substracting this number from pansize results in the supragenome above frequency border/n
float theo_curr_antisuprasize(float theta, float rho, uint64_t n, uint64_t border){
    long double v=0, sum1=0, sum2=0;
    uint64_t total = border;
    uint64_t steps = ((uint64_t) (border/50.));
    if (border < 100) {total = border; steps = border;}
    if (border > 10000){
      printf("Computing the Anti-Supragenome size...\tthis will take less time\n");
      printf("0%%");
      int zaehl;
      for (zaehl = 1; zaehl <= 50; zaehl++) printf(".");
      printf("100%%\n0%%");
    }
    uint64_t k;
    for (k=1; k < border+1.; k++){
      sum1 = sum1 + std::log((long double) (n-k+1.));
      sum2 = sum2 + std::log((long double) (n-k+rho));
      v = v + ( std::exp(sum1-sum2) )/k; 
      if (!(k % steps) )  {
	if (n > 10000){
	  printf("."); fflush(stdout);
	}
      }
    }
    if (n > 10000) printf("100%%\n");

    return theta*v;
}




// size of the supragenome
// this computes the number of genes in frequency greater than border/n  
// this is faster than computing the genes below frequency border/n and(!) more accurate
float integrand(float x,float rho){
    return 1.0/(x*std::pow((1.0-x),(1-rho)));
}

float funcintegrate(float (*f)(float x,float rho), float a, float b, int n, float rho){
int i;
double x;
double k = (b - a) / n;
double s = 0.5 * (f(a,rho) + f(b,rho));
    for (i = 1; i < n; i++) {
         x = a + k * i;
         s = s + f(x,rho);
    }
return s * k;
}

float theo_suprasize(float theta, float rho, float border){
    long double v=0;
    v = funcintegrate(integrand,border,1.0,1000000,rho); 
    return theta*v;
}










