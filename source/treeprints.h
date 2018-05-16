/***** treeprints.h ***********************************************
 * Description: functions to print trees etc. to the standard interface
 * Author: Franz Baumdicker, baumdicker@stochastik.uni-freiburg.de ; Peter Pfaffelhuber, pp@stochastik.uni-freiburg.de
 * File created 2010.
 *****************************************************************/

#include "treestructure.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
// using namespace GiNaC;


void printprobs(Node *node){
  float pp,c1p,c2p;
  int one = 0, two = 0, three=0;
  if (node->neighbor1 == node->parent) {pp = node->lossprob1 ; one = 1;}
  else if (node->neighbor2 == node->parent)  {pp = node->lossprob2 ; two = 1;}
  else if (node->neighbor3 == node->parent)  {pp = node->lossprob3 ; three = 1;}
  else {printf("WARNING: parent is not a neighbor\n"); }
  
  if (node->neighbor1 == node->child1)  {c1p = node->lossprob1 ; one = 1;}
  else if (node->neighbor2 == node->child1)  {c1p = node->lossprob2 ; two = 1;}
  else if (node->neighbor3 == node->child1) {c1p = node->lossprob3 ; three = 1;}
  else {printf("WARNING: child1 is not a neighbor\n"); }
  
  if (node->neighbor1 == node->child2) {c2p = node->lossprob1 ; one = 1;}
  else if (node->neighbor2 == node->child2) {c2p = node->lossprob2 ; two = 1;}
  else if (node->neighbor3 == node->child2) {c2p = node->lossprob3 ; three = 1;}
  else {printf("WARNING: child2 is not a neighbor\n"); }
  
  if (one+two+three != 3 && (node->child1) && (node->child2) ) printf("WARNING: relatives != neighbors\n");
  
  printf("id: \t %d \t dist: %.2f \t lossp: \t %.2f  \t lossp: \t %.2f  \t lossp: \t %.2f \n", node->id, node->disttoroot , pp, c1p, c2p);
}



void printlengths(Node *node){
  float pl,c1l,c2l;
  int one = 0, two = 0, three=0;
  if (node->neighbor1 == node->parent) {pl = node->length1 ; one = 1;}
  else if (node->neighbor2 == node->parent)  {pl = node->length2 ; two = 1;}
  else if (node->neighbor3 == node->parent)  {pl = node->length3 ; three = 1;}
  else {printf("WARNING: parent is not a neighbor\n"); }
  
  if (node->neighbor1 == node->child1)  {c1l = node->length1 ; one = 1;}
  else if (node->neighbor2 == node->child1)  {c1l = node->length2 ; two = 1;}
  else if (node->neighbor3 == node->child1) {c1l = node->length3 ; three = 1;}
  else {printf("WARNING: child1 is not a neighbor\n"); }
  
  if (node->neighbor1 == node->child2) {c2l = node->length1 ; one = 1;}
  else if (node->neighbor2 == node->child2) {c2l = node->length2 ; two = 1;}
  else if (node->neighbor3 == node->child2) {c2l = node->length3 ; three = 1;}
  else {printf("WARNING: child2 is not a neighbor\n"); }
  
  if (one+two+three != 3 && (node->child1) && (node->child2) ) printf("WARNING: relatives != neighbors\n");
  
  printf("id: \t %d \t dist: %.2f \t length: \t %.2f  \t length: \t %.2f  \t length: \t %.2f \n", node->id, node->disttoroot , pl, c1l, c2l);
}







void printrootedTree(Node *node) {
  int c1, c2, n1, n2, n3, p, i;
  i = c1 = c2 = n1 = n2 = n3 = p = -1;
  if(node) {
    i = node->id;
    if (node->parent) p = node->parent->id;
    if (node->child1) c1 = node->child1->id;
    if (node->child2) c2 = node->child2->id;
    
    printf("id: \t %d \t marked: %d \t parent: \t %d  \t child1: \t %d  \t child2: \t %d \n", i, node->marker , p, c1, c2);
    printprobs(node);
    printlengths(node);
//     printsymbprobs(node);
    printf("\n");
    printrootedTree(node->child1);
    printrootedTree(node->child2);
  }
}


void print(int feld[],int anzahl){
  int i;
     for(i=0;i<anzahl;++i)
         printf( "%d\t", feld[i]);
     printf("\n");
}



void printfloats(float feld[],int anzahl){
  int i;
     for(i=0;i<anzahl;++i)
         printf( "%.6f   ", feld[i]);
     printf("\n");
}



void printgfs(float feld[],int anzahl,float t, float r){
  int i;
     for(i=0;i<anzahl;++i)
         printf( "%4.1f\t", feld[i]*t/r);
     printf("\n");
}



