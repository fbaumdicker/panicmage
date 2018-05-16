/***** treefunctions.h ***********************************************
 * Description: functions operating on the given newick tree
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


// ersetze 
// genrand_real1()
//durch
// gsl_ran_flat(r,0.,1.)


// roots the tree at node, initialize with preroot NULL
void rootTree(Node *node,Node *preroot) {
 if(node != NULL){
   node->parent = preroot;
   if(node->neighbor1 == preroot){
     node->child1 = node->neighbor2;
     node->child2 = node->neighbor3;
     if (preroot != NULL) node->disttoroot = preroot->disttoroot + node->length1;
     else node->disttoroot = 0;
     rootTree(node->neighbor2, node);
     rootTree(node->neighbor3, node);
   }
   else if(node->neighbor2 == preroot){
     node->child1 = node->neighbor1;
     node->child2 = node->neighbor3;
     if (preroot != NULL) node->disttoroot = preroot->disttoroot + node->length2;
     else node->disttoroot = 0;
     rootTree(node->neighbor1, node);
     rootTree(node->neighbor3, node);
   }
   else if(node->neighbor3 == preroot){
     node->child1 = node->neighbor1;
     node->child2 = node->neighbor2;
     if (preroot != NULL) node->disttoroot = preroot->disttoroot + node->length3;
     else node->disttoroot = 0;
     rootTree(node->neighbor1, node);
     rootTree(node->neighbor2, node);
   }
   else { printf("WARNING: not able to root tree in node %d\n", node->id); }
 }
}


// finds the current root of the tree (node)
Node *findroot(Node *node){
  if (node != NULL){
    if (node->parent){
      findroot(node->parent);
    }
    else return node;
  }
}


// sets all lengths to same value
void setlengths(Node *node , float length){  // node should be root to get whole tree
if (node != NULL){
    if (node->neighbor1) node->length1 = length;
    if (node->neighbor2) node->length2 = length;
    if (node->neighbor3) node->length3 = length;
    if (node->child1) setlengths(node->child1, length);
    if (node->child2) setlengths(node->child2, length);
  }
}



//////////markings
void unmark(Node *node, int opt){
  if (node != NULL){
    if (opt == 0)  node->marker = 0;
    else if (opt == 1 && node->marker == 1) node->marker = 0;
    else if (opt == 2 && node->marker == 2) node->marker = 0;
    if (node->child1) unmark(node->child1,opt);
    if (node->child2) unmark(node->child2,opt);
  }
}

void unmarkTree(Node *node, int opt){
  Node *var;
  var = findroot(node);
  unmark(var,opt);
}

// marks the path up to the root or the next marked node and returns the length up to there
float marktoroot(Node *node){
  float subtreelength;
  if (node != NULL){
    node->marker = 1;
    if (node->parent) {
      if (node->parent->marker == 0){
	subtreelength = marktoroot(node->parent);
	subtreelength += node->disttoroot - node->parent->disttoroot;
	return subtreelength;
      }
      else{
       return node->disttoroot - node->parent->disttoroot;
      }
      
    }
    else return node->disttoroot;
  }
}

//marks a complete subtree and returns its length
float marksubtree (int *subtree, int subn, Node *tree){
  float subtreelength = 0;
  rootTree(tree+subtree[0],NULL);
  unmark(tree+subtree[0],0);
  int i;
  for(i = 0; i < subn; i++){
    subtreelength += marktoroot(tree+subtree[i]);
  }
  return subtreelength;
}


////////////////////////////////


/////////////////probabilities

void unprob(Node *node){  // node should be root to unprob whole tree
  if (node != NULL){
      node->lossprob1 = -1; node->lossprob2 = -1; node->lossprob3 = -1;
      if (node->child1) unprob(node->child1);
      if (node->child2) unprob(node->child2);
    }
}

// computes the probabilities to loose/keep a gene in one direction ...
float computeProbs(Node *node, Node *child, float rho){
  float factor;
  if (child == NULL){
    if (node->lossprob1 < 0 || node->lossprob2 < 0 || node->lossprob3 < 0){
      if (node->neighbor1 == child) node->lossprob1 = 0;
      if (node->neighbor2 == child) node->lossprob2 = 0;
      if (node->neighbor3 == child) node->lossprob3 = 0;
    }
    return 0;
  }
  else{
      if (node->neighbor1 == child){
	if (node->lossprob1 >= 0) return node->lossprob1;
	else{
	    factor = 1.0;
	    factor = factor * computeProbs(child,child->child1,rho);
	    factor = factor * computeProbs(child,child->child2,rho);
	    node->lossprob1 = ( 1-std::exp(-rho/2*(node->length1)) ) + std::exp(-rho/2*(node->length1)) * factor;
	    return node->lossprob1;
	}
      }
      if (node->neighbor2 == child){
	if (node->lossprob2 >= 0) return node->lossprob2;
	else{
	    factor = 1.0;
	    factor = factor * computeProbs(child,child->child1,rho);
	    factor = factor * computeProbs(child,child->child2,rho);
	    node->lossprob2 = ( 1-std::exp(-rho/2*(node->length2)) ) + std::exp(-rho/2*(node->length2)) * factor;
	    return node->lossprob2;
	  
	}
      } 
      if (node->neighbor3 == child){
	if (node->lossprob3 >= 0) return node->lossprob3;
	else{
	    factor = 1.0;
	    factor = factor * computeProbs(child,child->child1,rho);
	    factor = factor * computeProbs(child,child->child2,rho);
	    node->lossprob3 = ( 1-std::exp(-rho/2*(node->length3)) ) + std::exp(-rho/2*(node->length3)) * factor;
	    return node->lossprob3;
	  
	}
      }
  }
  
}


// ...for all possibile directions
void computeprobsall(Node *tree, int leaves, float rho){
  int i;
  for(i = 0; i < leaves; i++){
    rootTree(tree+i,NULL);
    computeProbs(tree+i,(tree+i)->child1,rho);
    computeProbs(tree+i,(tree+i)->child2,rho);
  }
} 


// marks the remaining part of the tree and computes the probability to loose a gene ...at all remaining leaves
float losstoleaves(Node *node){
  if ( !((node)->marker) ){
      node->marker = 2;
      if (node->parent->marker == 0) {
	return losstoleaves(node->parent);
      }
      else if(node->parent->marker == 1){
        if (node->parent->neighbor1 == node) return node->parent->lossprob1;
        else if (node->parent->neighbor2 == node) return node->parent->lossprob2;
        else if (node->parent->neighbor3 == node) return node->parent->lossprob3;
      }
      else return 1;
  }
  else return 1;
}

// ...at all remaining leaves
float totalloss(Node *tree, int leaves){
  int i;
  float prob = 1;
  for(i = 0; i < leaves; i++){
    prob *= losstoleaves(tree+i);
  }
  return prob;
}


//////////////////////////////////////////////

////////////////////gene frequency spectrum
// add the expected number genes present in k of n indiviudal for all possible combinations
void koutofn(int feld[],int n,int k,int pos,int val, Node *tree, float *zeiger, float rho){
  int i;
  if(pos==k){
      unmarkTree(tree,0);
      float sublength,lossprob,keepprob;
      int subn = k;
      int *subtree;
      subtree = feld;
      sublength = marksubtree(subtree,subn,tree);
      keepprob = std::exp(-rho/2*sublength);
      lossprob = totalloss(tree,n);
      *zeiger += keepprob*lossprob;
      //printf("keepprob:\t %.2f \t lossprob: \t %.5f \t product: \t %.5f \n", keepprob, lossprob, keepprob*lossprob);
      return; // print(feld,k);
  }
  for(i=val;i<n;++i){
      feld[pos]=i;
      koutofn(feld,n,k,pos+1,i+1,tree,zeiger,rho);
  }
}



// compute the gene frequency spectrum for the given tree
void treegfs(Node *tree, int leaves, float *gfs_k, float rho){
  int k;
  int feld[leaves-1];
  float *zeiger;
  for(k = 1; k <= leaves; k++){
    //printf("%d aus %d\n", k, leaves);
    zeiger = &gfs_k[k-1];
    gfs_k[k-1] = 0;
    koutofn(feld,leaves,k,0,0,tree,zeiger,rho);
  }
}

////////////////////////////////////////////////////




// returns the time to the root for a node
float depthofnode(Node *node, float sofar){
  if (node != NULL){
    if (node->parent != NULL){
      if (node->parent == node->neighbor1){
	sofar += node->length1;
      }
      else if (node->parent == node->neighbor2){
	sofar += node->length2;
      }
      else if (node->parent == node->neighbor3){
	sofar += node->length3;
      }
      else{ printf("WARNING: parent is not a neighbor\n"); }
      return depthofnode(node->parent,sofar);
    }
    else{
      return sofar;
    }
  }
}



// scales all lengths of a tree by factor
void scalerootedtree(Node *node , float factor){  // node should be root to get whole tree
  if (node != NULL){
    node->disttoroot *= factor;
    if (node->neighbor1) node->length1 *= factor;
    if (node->neighbor2) node->length2 *= factor;
    if (node->neighbor3) node->length3 *= factor;
    if (node->child1) scalerootedtree(node->child1, factor);
    if (node->child2) scalerootedtree(node->child2, factor);
  }
}



// finds the longest distance in the tree
float findmaxdepth(Node *tree, int leaves){
  float mem = 0.;
  int i, max_id;
  for(i = 0; i<leaves; i++){
    if (mem <= depthofnode(tree+i,0.)){
    mem = depthofnode(tree+i,0.);
    max_id = (tree+i)->id;
    }  
  }
  //printf("maximum depth at Node %d:\t%f\n", max_id, mem );
  return mem;
}

// determines the treeheight of a tree
float treeheight(Node *tree, int leaves){
rootTree(tree,NULL);
return findmaxdepth(tree,leaves)/2.;
}



// removes the root of a tree in order to get an unrooted tree (later the root will always be one of the leaves)
void removeroot(Node *root){
  root->child1->length3 += root->child2->length3; 
  root->child2->length3 = root->child1->length3;
  root->child1->parent = root->child2;
  root->child1->neighbor3 = root->child2;
  root->child2->parent = root->child1;
  root->child2->neighbor3 = root->child1;
}




// generates a random tree according to the Kingman's coalescent
void maketree(Node *tree, Node **list, int leaves){ 
  /* tree should be leave number one 
     note that disttoroot is used as distance to the leave in this function*/
  int i,j;
  double t;

  Node *root;
  int n = leaves;
  root = tree+2*n-2;
  root->parent = NULL;

  /* set all values zero */
  for(i=0; i<2*n-1; i++){
    tree[i].lossprob1 = 0.;
    tree[i].length1 = 0.;
    tree[i].lossprob2 = 0.;
    tree[i].length2 = 0.;
    tree[i].lossprob3 = 0.;
    tree[i].length3 = 0.;
    tree[i].marker = 0;
    tree[i].disttoroot = 0.;
  }

  /* set leaves children NULL  */
  for(i=0; i<leaves; i++){
    tree[i].child1 = NULL;
    tree[i].child2 = NULL;
  }

  /* initialize list */
    for(i=0;i<n;i++){
      //printf("%d\n", i);
      list[i] = tree + i;
    }

  /* set ids of all nodes*/
  for(i=0; i<2*n-1; i++){
    tree[i].id = i;
  }

  /* generate distance to root of internal nodes */
    t = 0.0;
    for(i=n; i>1; i--){
      t += -2.0*log(1.0-gsl_ran_flat(r,0.,1.))/((double)(i*(i-1)));
      //t += 1.;
      tree[2*n-i].disttoroot = t;
    }

  /* generate tree topology */
    for(i=n;i>1;i--){
      j = i*gsl_ran_flat(r,0.,1.);
      list[j]->parent = tree + 2 * n - i;
      tree[2*n-i].child1 = list[j];
      list[j] = list[i-1];
      j = (i-1)*gsl_ran_flat(r,0.,1.);
      list[j]->parent = tree + 2 * n - i;
      tree[2*n-i].child2 = list[j];
      list[j] = tree + 2 * n - i;
    }

  /* set neighbors*/
  for(i=0; i<2*n-1; i++){
    tree[i].neighbor1 = tree[i].child1;
    tree[i].neighbor2 = tree[i].child2;
    tree[i].neighbor3 = tree[i].parent;
  }

  /* set lengths */
  for(i=0; i<2*n-1; i++){
    if (tree[i].parent) tree[i].length3 = (tree+i)->parent->disttoroot - (tree+i)->disttoroot;
    if (tree[i].child1) tree[i].length1 = (tree+i)->child1->length3;
    if (tree[i].child2) tree[i].length2 = (tree+i)->child2->length3;
  }

  /* remove the root*/
  removeroot(root);
  // that is:
  // root->child1->length3 += root->child2->length3; 
  // root->child2->length3 = root->child1->length3;
  // root->child1->parent = root->child2;
  // root->child1->neighbor3 = root->child2;
  // root->child2->parent = root->child1;
  // root->child2->neighbor3 = root->child1;

  rootTree(tree,NULL);
}


// empty tree this function deletes any information which is given in the tree.
// just used for debugging
void emptytree(Node *node){  // node should be root to unprob whole tree
  int inc;
  if (node != NULL){
      node->lossprob1 = -1; node->lossprob2 = -1; node->lossprob3 = -1;
      if (node->child1) unprob(node->child1);
      if (node->child2) unprob(node->child2);
      node->id = 0;
//   struct Node *parent;
  node->parent = NULL;
  node->child1 = NULL;
  node->child2 = NULL;
  node->neighbor1 = NULL;
  node->neighbor2 = NULL;
  node->neighbor3 = NULL;
  
  for(inc = 0;inc < 500; inc++){
    node->pkfh1_numeric[inc] = 0;
  }
  node->lossprob1 = 0;
  node->length1 = 0;
  for(inc = 0;inc < 500; inc++){
    node->pkfh2_numeric[inc] = 0;
  }
  node->lossprob2 = 0;
  node->length2 = 0;
  node->lossprob3 = 0;
  node->length3 = 0;
  node->marker = 0;
  node->disttoroot = 0;
  node->compdone = 0;
    }
}


// enlarges the branches in order to test for sampling bias
void enlargebranches(Node *tree, int leaves, float addsize){
  int i; 
  rootTree(tree+1,NULL);
  for(i = 0; i < leaves; i++){
    if (tree[i].neighbor1) tree[i].length1 += addsize;
    if (tree[i].neighbor2) tree[i].length2 += addsize;
    if (tree[i].neighbor3) tree[i].length3 += addsize;
    
    if ((tree+i)->parent->neighbor1 == (tree+i)) (tree+i)->parent->length1 += addsize;
    if ((tree+i)->parent->neighbor2 == (tree+i)) (tree+i)->parent->length2 += addsize;
    if ((tree+i)->parent->neighbor3 == (tree+i)) (tree+i)->parent->length3 += addsize;
    if (i == 0) rootTree(tree,NULL);
  }
  rootTree(tree,NULL);
}


// generates the additional time the branches have due to sampling bias
float generate_random_enlargementsize(int leaves, int enlargenumber){
  /* generate distance to root of internal nodes */
    int i;
    float t = 0.0;
    for(i=leaves+1; i < enlargenumber+1; i++){
      t += -2.0*log(1.0-gsl_ran_flat(r,0.,1.))/((double)(i*(i-1)));
    }
    return t;
}


