
#include "treestructure.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





// marks the remaining part of the tree and computes the probability to loose a gene ...at all remaining leaves
ex losstoleaves_symbolic(Node *node){
  if ( !((node)->marker) ){
      node->marker = 2;
      if (node->parent->marker == 0) {
	return losstoleaves_symbolic(node->parent);
      }
      else if(node->parent->marker == 1){
        if (node->parent->neighbor1 == node) return node->parent->symblossprob1;
        else if (node->parent->neighbor2 == node) return node->parent->symblossprob2;
        else if (node->parent->neighbor3 == node) return node->parent->symblossprob3;
      }
      else return 1;
  }
  else return 1;
}

// ...at all remaining leaves
ex totalloss_symbolic(Node *tree, int leaves){
  int i;
  ex prob = 1;
  for(i = 0; i < leaves; i++){
    prob *= losstoleaves_symbolic(tree+i);
  }
  return prob;
}





////////////////////gene frequency spectrum
// add the expected number genes present in k of n indiviudal for all possible combinations
//there is a much faster way to compute the gfs by treegfs_symbolic_fast
void koutofn_symbolic(int feld[],int n,int k,int pos,int val, Node *tree, ex *zeiger, const ex & rhoS){
  int i;
  if(pos==k){
      unmarkTree(tree,0);
      float sublength;
      ex keepprobS;
      ex lossprobS;
      int subn = k;
      int *subtree;
      subtree = feld;
      sublength = marksubtree(subtree,subn,tree);
      keepprobS = exp(-rhoS/2*sublength);
      lossprobS = totalloss_symbolic(tree,n);
      *zeiger += keepprobS*lossprobS;
      //printf("keepprob:\t %.2f \t lossprob: \t %.5f \t product: \t %.5f \n", keepprob, lossprob, keepprob*lossprob);
      return; // print(feld,k);
  }
  for(i=val;i<n;++i){
      feld[pos]=i;
      koutofn_symbolic(feld,n,k,pos+1,i+1,tree,zeiger,rhoS);
  }
}




// compute the gene frequency spectrum for the given tree
//there is a much faster way to compute the gfs by treegfs_symbolic_fast
void treegfs_symbolic_slow(Node *tree, int leaves, ex *gfs_k_symb, const ex & rhoS){
  int k;
  int feld[leaves-1];
  ex *zeiger;
  for(k = 1; k <= leaves; k++){
    printf("%d aus %d\n", k, leaves);
    zeiger = &gfs_k_symb[k-1];
    gfs_k_symb[k-1] = 0;
    koutofn_symbolic(feld,leaves,k,0,0,tree,zeiger,rhoS);
  }
}

////////////////////////////////////////////////////

void unprob_symb(Node *node){  // node should be root to unprob whole tree
  if (node != NULL){
      node->symblossprob1 = (ex) -1; node->symblossprob2 = (ex) -1; node->symblossprob3 = (ex) -1;
      if (node->child1) unprob_symb(node->child1);
      if (node->child2) unprob_symb(node->child2);
    }
}

// computes the probabilities to loose/keep a gene in one direction ...
ex computeProbs_symbolic(Node *node, Node *child, const ex & rhoS){
  ex factor;
  if (child == NULL){
    if (node->symblossprob1 == -1 || node->symblossprob2 == -1 || node->symblossprob3 == -1 ){ 
      if (node->neighbor1 == child) node->symblossprob1 = 0;
      if (node->neighbor2 == child) node->symblossprob2 = 0;
      if (node->neighbor3 == child) node->symblossprob3 = 0;
    }
    return 0;
  }
  else{
      if (node->neighbor1 == child){
	if (node->symblossprob1 != -1) return node->symblossprob1;
	else{
	    factor = 1;
	    factor = factor * computeProbs_symbolic(child,child->child1,rhoS);
	    factor = factor * computeProbs_symbolic(child,child->child2,rhoS);
	    node->symblossprob1 = ( 1-exp(-rhoS/2*(node->length1)) ) + exp(-rhoS/2*(node->length1)) * factor;
	    return node->symblossprob1;
	}
      }
      if (node->neighbor2 == child){
	if (node->symblossprob2 != -1) return node->symblossprob2;
	else{
	    factor = 1;
	    factor = factor * computeProbs_symbolic(child,child->child1,rhoS);
	    factor = factor * computeProbs_symbolic(child,child->child2,rhoS);
	    node->symblossprob2 = ( 1-exp(-rhoS/2*(node->length2)) ) + exp(-rhoS/2*(node->length2)) * factor;
	    return node->symblossprob2;
	  
	}
      } 
      if (node->neighbor3 == child){
	if (node->symblossprob3 != -1) return node->symblossprob3;
	else{
	    factor = 1;
	    factor = factor * computeProbs_symbolic(child,child->child1,rhoS);
	    factor = factor * computeProbs_symbolic(child,child->child2,rhoS);
	    node->symblossprob3 = ( 1-exp(-rhoS/2*(node->length3)) ) + exp(-rhoS/2*(node->length3)) * factor;
	    return node->symblossprob3;
	  
	}
      }
  }
  
}


// ...for all possibile directions
void computeprobsall_symbolic(Node *tree, int leaves, const ex & rhoS){
  int i;
  for(i = 0; i < leaves; i++){
    printf("%d\n",i);
    rootTree(tree+i,NULL);
    printf("rooted the tree\n");
    computeProbs_symbolic(tree+i,(tree+i)->child1,rhoS);
    printf("child1\n");
    computeProbs_symbolic(tree+i,(tree+i)->child2,rhoS);
    printf("child2\n");
  }
} 






// new computation

void initialize_tree(Node *node, int leaves){  // node should be root to get whole tree
  if (node != NULL){
//     cout << "in initialize_tree\n";
      node->compdone = 0;
      int counter;
      int isleave = 0;
      // set all expressions to -1
      for (counter = 0; counter < 500; counter++){
	node->pkfh1[counter] = (ex) -1; node->pkfh2[counter] = (ex) -1; node->pkfh3[counter] = (ex) -1;
      }
//       cout << "so far ok\n";
      // set all expressions which are at most used to 0
      for (counter = 0; counter <= leaves; counter++){
	node->pkfh1[counter] = (ex) 0; node->pkfh2[counter] = (ex) 0; node->pkfh3[counter] = (ex) 0;
      }
      if (node->child1 || node->child2){
	
      }else{ // this is a leave
	isleave = 1;
	node->compdone = 1;
// 	cout << "Here as well:\t" << node->id << "\n";
      }
            if (node->child1){
	initialize_tree(node->child1,leaves);
      }else{
	node->pkfh1[1] = (ex) 1;
// 	cout << node->id << "\n";
// 	cout << "The Computation has been done?:\t" << node->compdone << "\n";
      }
      if (node->child2){
	initialize_tree(node->child2,leaves);
      }else{
	node->pkfh2[1] = (ex) 1;
	if (isleave == 1){node->pkfh2[1] = (ex) 0; node->pkfh2[0] = (ex) 1; };
// 	cout << "";
// 	cout << node->disttoroot << "\n";
// 	cout << "THE Computation has been done?:\t" << node->compdone << "\n";
      }
    }
}



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





void comp_pkfhs(Node *node, int leaves, const ex & rhoS){
  int k,kprime;
  ex summe;
  float length;
//   cout << "starting pkfhs \n";
  if (node->child1){
    if (node->child1->compdone == 0){
      comp_pkfhs(node->child1,leaves,rhoS);
    }
  }
  if (node->child2){
    if (node->child2->compdone == 0){
      comp_pkfhs(node->child2,leaves,rhoS);
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
	summe += node->child1->pkfh1[kprime] * node->child1->pkfh2[k-kprime];
      }
//       cout << "so far rhoS not used\n";
      node->pkfh1[k] = exp(-rhoS/2*length) * summe;
//       cout << "now rhoS is used\n";
    }
    node->pkfh1[0] += 1-exp(-rhoS/2*length);
  }
  // same for second child if it exists (the root has only one child)
  if(node->child2){
    length = node->child2->disttoroot - node->disttoroot;
    for (k = 0; k<=leaves; k++){
      //compute sum
      summe = 0;
      for (kprime = 0; kprime <= k; kprime++){
	summe += node->child2->pkfh1[kprime] * node->child2->pkfh2[k-kprime];
      }
      node->pkfh2[k] = exp(-rhoS/2*length) * summe;
    }
    node->pkfh2[0] += 1-exp(-rhoS/2*length);
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



void check_probs(Node *node, int leaves, const ex & rhoS){
  ex summe1 = 0;
  ex summe2 = 0;
  float testerrho = 0;
  int k;
  for(k=0;k<=leaves;k++){
    summe1 += node->pkfh1[k];
  }
  for(k=0;k<=leaves;k++){
    summe2 += node->pkfh2[k];
  }
  cout << "Knoten:" << node->id << "\n" << summe1.subs(rhoS == testerrho).evalf() << "\n" << summe2.subs(rhoS == testerrho).evalf() << "\n";
  cout << "Knoten:" << node->id << "\n";
  for(k=0;k<=leaves;k++){
   cout << node->pkfh1[k].subs(rhoS == testerrho).evalf()  << "\t";
  }
  cout << "---------------------------------------------------------------------\n";
  for(k=0;k<=leaves;k++){
   cout << node->pkfh2[k].subs(rhoS == testerrho).evalf()  << "\t";
  }
  cout << endl;
  
  if(node->child1){cout << "Es gibt Kind 1\n"; }
  else{cout << "1: NICHT!!\n";}
  if(node->child2){cout << "Es gibt Kind 2\n"; }
  else{cout << "2: NICHT!!\n";}
  if(node->child1){
    check_probs(node->child1,leaves,rhoS);
  }
  if(node->child2){
    check_probs(node->child2,leaves,rhoS);
  }
}


void addtogfs(Node *node, int leaves, ex *zeiger, const ex & rhoS, int k){
    // get length to the parent to decide how many new genes one would expect at this node
    float length = -1.0;
    ex expectednewatthisnode;
    ex keepprob = 0;
    int subk;
    if(node->parent){
      length = node->disttoroot - node->parent->disttoroot;
//       expectednewatthisnode = 1/rhoS * (1 - exp(-rhoS/2*length) ); 
      expectednewatthisnode = 1/1 * (1 - exp(-rhoS/2*length) ); 
    }else{
      // in this case the length is given by infinity
//       expectednewatthisnode = 1/rhoS; 
      expectednewatthisnode = 1/1; 
    }
    // compute the expected number of genes in k out of leaves individuals 
    for(subk = 0; subk<= k; subk++){
      keepprob += node->pkfh1[subk] * node->pkfh2[k-subk];
    }
    *zeiger += expectednewatthisnode * keepprob;
    if(node->child1){
      addtogfs(node->child1, leaves, zeiger,rhoS,k);
    }
    if(node->child2){
      addtogfs(node->child2, leaves, zeiger,rhoS,k);
    }
}



void treegfs_symbolic_fast(Node *tree, int leaves, ex *gfs_k_symb, const ex & rhoS){
  int k;
  ex *zeiger;
//   printf("called treegfs_symbolic_fast\n");
  for(k = 1; k <= leaves; k++){
//     printf("%d aus %d\n", k, leaves);
    zeiger = &gfs_k_symb[k-1];
    gfs_k_symb[k-1] = 0;
    addtogfs(tree,leaves,zeiger,rhoS,k);
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
double my_f_symbolic (const gsl_vector *v, void *params)
{
  extern int g_includecoreflag;
  double rho, theta, sum = 0.;
  Params_symbolic *para = (Params_symbolic *)params;
       
  int anzahl = para->anzahl;
//   Node *tree = para->tree;
  
  rho = gsl_vector_get(v, 0); 
  theta = gsl_vector_get(v, 1);
  
//   cout << "igot this rho:" << rho << "\n";
  if (rho == 0){
    rho = 0.00000000000000001;
//     cout << "but i use this rho:" << rho << "\n";
  }
//   cout << "igot this theta:" << theta << "\n";
  

  ex * mysymbolicgfs;
  mysymbolicgfs = new ex[anzahl];
  
  int inc;
  for (inc = 0; inc < anzahl; inc++) {
    mysymbolicgfs[inc] = para->symbolicgfs[inc];
  }
  
//   printgfs_thin_symb(mysymbolicgfs,anzahl);
    
  
  float *theogfs_float;
  theogfs_float = new float[anzahl];
  for (inc = 0; inc < anzahl; inc++){
//      printf("%d\n",inc);
     theogfs_float[inc] = to_double(ex_to<numeric>(mysymbolicgfs[inc].subs(para->rhoS == rho).evalf()));
  }
  
  //old
//   float *theogfs;
//   theogfs = (float *)malloc(anzahl*sizeof(float));
//   treegfs(tree,anzahl,theogfs,rho);
  
  float *datagfs;
  datagfs = (float *)malloc(anzahl*sizeof(float));
  
  int i;
  for (i = 0; i < anzahl; i++) {
    datagfs[i] = para->datagfs[i];
  }
  
//   printgfs(datagfs,11,1.0,1.0);
//   printgfs(theogfs_float,11,1.0,1.0);
 
  //return lsqu(datagfs,theogfs,anzahl-1,theta,rho);
//   cout << "This is the llh:" <<  llh(datagfs, theogfs_float, anzahl-1,theta,rho) << "\n";
  if(g_includecoreflag == 1){
    return llh(datagfs, theogfs_float, anzahl,theta,rho);
  }else{
    return llh(datagfs, theogfs_float, anzahl-1,theta,rho);
  }
}





// // try to compute the gfs in a clever way:
// rootTree(intree,NULL);
// initialize_tree(intree,leaves);
// check_probs(intree,leaves,x);
// 
// printf("_____________________________\n");
// comp_pkfhs(intree,leaves,x);
// check_probs(intree,leaves,x);
// 
//   ex *theogfs_symb_fast;
//   theogfs_symb_fast = new ex[leaves];
//   treegfs_symbolic_fast(intree,leaves,theogfs_symb_fast,x);
//   for (inc = 0; inc < leaves; inc++){
//     theogfs_symb_fast[inc] =  theogfs_symb_fast[inc].subs(x == 0.5).evalf();
//   }
//   printgfs_thin_symb(theogfs_symb_fast,leaves);
//   printf("_--------------------------------------\n");
//   printgfs_thin_symb(theogfs_symb,leaves);
//   printf("_--------------------------------------\n");
//   printgfs(theogfs_numeric,leaves,1.0,0.5);










/*estimate the parameters theta and rho using the symbolic expression*/

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

void estimate_symbolic(float *theta_hat, float *rho_hat, Params_symbolic *paraS){   
           
       const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
       gsl_multimin_fminimizer *s = NULL;
       gsl_vector *ss, *x;
       gsl_multimin_function minex_func;
     
       size_t iter = 0;
       int status;
       double size;
       
//        cout << "so far fine\n";
       
       /*optional set the starting points to a first simple estimate*/
       
       int i, n;
       float c1, c2, firsttheta, firstrho;
       c1 = 0.;
       c2 = 0.;
       n = paraS->anzahl;
       for (i = 0; i < n; i++){
       c1 += (i+1.)/n * paraS->datagfs[i];
       c2 += (i+1.)*(n-i-1.)/(n*(n-1.)) * paraS->datagfs[i];
       }
       firstrho = c2/(c1-c2);
       firsttheta = c2 + firstrho*c2;
       if (firstrho > 20.) firstrho = 20.;
       if (firstrho < 0.1) firstrho = 0.1;
       if (firsttheta > 10000.) firsttheta = 10000.;
       if (firsttheta < 5.) firsttheta = 5.;
       
       
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
       minex_func.f = my_f_symbolic;
       minex_func.params = paraS;
     
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
               printf("converged to minimum at\n");
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













