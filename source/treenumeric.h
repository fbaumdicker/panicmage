
#include "treestructure.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>






void initialize_tree_numeric(Node *node, int leaves){  // node should be root to get whole tree
  if (node != NULL){
//     cout << "in initialize_tree\n";
      node->compdone = 0;
      int counter;
      int isleave = 0;
      // set all expressions to -1
      for (counter = 0; counter < 500; counter++){
	node->pkfh1_numeric[counter] = -1.0; node->pkfh2_numeric[counter] = -1.0;
      }
//       cout << "so far ok\n";
      // set all expressions which are at most used to 0
      for (counter = 0; counter <= leaves; counter++){
	node->pkfh1_numeric[counter] = 0.0; node->pkfh2_numeric[counter] = 0.0;
      }
      if (node->child1 || node->child2){
	
      }else{ // this is a leave
	isleave = 1;
	node->compdone = 1;
// 	cout << "Here as well:\t" << node->id << "\n";
      }
            if (node->child1){
	initialize_tree_numeric(node->child1,leaves);
      }else{
	node->pkfh1_numeric[1] = 1.0;
// 	cout << node->id << "\n";
// 	cout << "The Computation has been done?:\t" << node->compdone << "\n";
      }
      if (node->child2){
	initialize_tree_numeric(node->child2,leaves);
      }else{
	node->pkfh2_numeric[1] = 1.0;
	if (isleave == 1){node->pkfh2_numeric[1] = 0.0; node->pkfh2_numeric[0] = 1.0; };
// 	cout << "";
// 	cout << node->disttoroot << "\n";
// 	cout << "THE Computation has been done?:\t" << node->compdone << "\n";
      }
    }
}




void comp_pkfhs_numeric(Node *node, int leaves, float myrho){
  int k,kprime;
  float summe;
  float length;
//   cout << "starting pkfhs \n";
  if (node->child1){
    if (node->child1->compdone == 0){
      comp_pkfhs_numeric(node->child1,leaves,myrho);
    }
  }
  if (node->child2){
    if (node->child2->compdone == 0){
      comp_pkfhs_numeric(node->child2,leaves,myrho);
    }
  }
  // now everything should be ready
  if(node->child1){
//       cout << "got to the one to compute\n";
    length = node->child1->disttoroot - node->disttoroot;
    for (k = 0; k<=leaves; k++){
      //compute sum
      summe = 0;
      for (kprime = 0; kprime <= k; kprime++){
	summe += node->child1->pkfh1_numeric[kprime] * node->child1->pkfh2_numeric[k-kprime];
      }
//       cout << "so far rhoS not used\n";
      node->pkfh1_numeric[k] = std::exp(-myrho/2*length) * summe;
//       cout << "now rhoS is used\n";
    }
    node->pkfh1_numeric[0] += 1-std::exp(-myrho/2*length);
  }
  // same for second child if it exists (the root has only one child)
  if(node->child2){
    length = node->child2->disttoroot - node->disttoroot;
    for (k = 0; k<=leaves; k++){
      //compute sum
      summe = 0;
      for (kprime = 0; kprime <= k; kprime++){
	summe += node->child2->pkfh1_numeric[kprime] * node->child2->pkfh2_numeric[k-kprime];
      }
      node->pkfh2_numeric[k] = std::exp(-myrho/2*length) * summe;
    }
    node->pkfh2_numeric[0] += 1-std::exp(-myrho/2*length);
  }
  
  
}


void addtogfs_numeric(Node *node, int leaves, float *zeiger, float myrho, int k){
    // get length to the parent to decide how many new genes one would expect at this node
    float length = -1.0;
    float expectednewatthisnode;
    float keepprob = 0.0;
    int subk;
    if(node->parent){
      length = node->disttoroot - node->parent->disttoroot;
//       expectednewatthisnode = 1/rhoS * (1 - exp(-rhoS/2*length) ); 
      expectednewatthisnode = 1.0/1.0 * (1.0 - std::exp(-myrho/2*length) ); 
    }else{
      // in this case the length is given by infinity
//       expectednewatthisnode = 1/rhoS; 
      expectednewatthisnode = 1.0/1.0; 
    }
    // compute the expected number of genes in k out of leaves individuals 
    for(subk = 0; subk<= k; subk++){
      keepprob += node->pkfh1_numeric[subk] * node->pkfh2_numeric[k-subk];
    }
    *zeiger += expectednewatthisnode * keepprob;
    if(node->child1){
      addtogfs_numeric(node->child1, leaves, zeiger,myrho,k);
    }
    if(node->child2){
      addtogfs_numeric(node->child2, leaves, zeiger,myrho,k);
    }
}


void treegfs_numeric_fast(Node *tree, int leaves, float *gfs_k_numeric, float myrho){
  int k;
  float *zeiger;
  for(k = 1; k <= leaves; k++){
//     printf("%d aus %d\n", k, leaves);
    zeiger = &gfs_k_numeric[k-1];
    gfs_k_numeric[k-1] = 0;
    addtogfs_numeric(tree,leaves,zeiger,myrho,k);
  }
}



/*------------function which is minimized-----------------------------------*/
double my_f_numeric (const gsl_vector *v, void *params)
{
  extern int g_includecoreflag;
//   cout << "This is the numeric estimation...\n";
  double rho, theta, sum = 0.;
  Params *para = (Params *)params;
       
  int anzahl = para->anzahl;
  Node *tree = para->tree;
  
  rho = gsl_vector_get(v, 0); 
  theta = gsl_vector_get(v, 1);
  
  rootTree(tree,NULL);
  unprob(tree);
  initialize_tree_numeric(tree,anzahl);
  comp_pkfhs_numeric(tree,anzahl,rho);
  
  
//   rootTree(tree,NULL);
  float *theogfs;
  theogfs = (float *)malloc(anzahl*sizeof(float));
//   treegfs(tree,anzahl,theogfs,rho);
  treegfs_numeric_fast(tree,anzahl,theogfs,rho);
  
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









/*estimate the parameters theta and rho using the numeric precomputation
 of pkfhs (probabilities to keep from here in k individuals) */

// Params are:
// Params *para;
// para = malloc(sizeof(Params) + number * sizeof(double));
// para->anzahl = number;
// para->tree = tree;
// }

void estimate_numeric(float *theta_hat, float *rho_hat, Params *para, int customstart_flag,float customstarttheta, float customstartrho){   
           
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
       
       
       if (customstart_flag == 1){
        firsttheta = customstarttheta;
        firstrho = customstartrho;
        printf("Using custom start values for numeric parameter estimation:\n");
        printf("theta start value = %.0f\t rho start value = %.2f\n",firsttheta, firstrho);
        }
       
       /* Starting point */
       x = gsl_vector_alloc (2);
       
//         gsl_vector_set (x, 0, 1.);  // fixed starting point rho  
//         gsl_vector_set (x, 1, 2000.0); // fixed starting point theta
       
       gsl_vector_set (x, 0, firstrho);  // estimated starting point rho
       gsl_vector_set (x, 1, firsttheta); // estimated starting point theta

//        printf("firsttheta = %.0f\tfirstrho = %.2f\n",firsttheta, firstrho);
       
     
       /* Set initial step sizes to 1 */
       ss = gsl_vector_alloc (2);
       gsl_vector_set_all (ss, 1.0);
     
       /* Initialize method and iterate */
       minex_func.n = 2;
       minex_func.f = my_f_numeric;
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
               printf ("converged to minimum at\n");
             }
     
        printf ("%5zu %10.3e %10.3e f() = %7.3f size = %.3f\n", 
                   iter,
                   gsl_vector_get (s->x, 0), 
                   gsl_vector_get (s->x, 1), 
                   s->fval, size);
         }
       while(status == GSL_CONTINUE && iter < 100);
       
//        printf("\n");
       if (iter != 100) {
     theta_hat[0] = (float) gsl_vector_get (s->x, 1);
     rho_hat[0]   = (float) gsl_vector_get (s->x, 0);
       }
       else {
     printf("WARNING: minimizer did not converge: \t size: %.2f \t\t theta: %.2f \t\t rho: %.2f\n", size, (float) gsl_vector_get (s->x, 1),  (float) gsl_vector_get (s->x, 0) );
//   theta_hat[0] = 0;
//   rho_hat[0] = 0;
     theta_hat[0] = (float) gsl_vector_get (s->x, 1);
     rho_hat[0]   = (float) gsl_vector_get (s->x, 0);
       }
       gsl_vector_free(x);
       gsl_vector_free(ss);
       gsl_multimin_fminimizer_free (s);
}









