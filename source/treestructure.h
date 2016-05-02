/***** treestructure.h ***********************************************
 * Author: Franz Baumdicker, baumdicker@stochastik.uni-freiburg.de ; Peter Pfaffelhuber, pp@stochastik.uni-freiburg.de
 * File created 2010.
 *****************************************************************/


#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

typedef struct Testsymbol{ // TODO this Testsymbol can be removed 
  ex ausdruck;
} Testsymbol;

typedef struct Node{
  int id;
  struct Node *parent;
  struct Node *child1;
  struct Node *child2;
  struct Node *neighbor1;
  ex symblossprob1;
  ex pkfh1[500];
  float pkfh1_numeric[500];
  float lossprob1;
  float length1;
  struct Node *neighbor2;
  ex symblossprob2;
  ex pkfh2[500];
  float pkfh2_numeric[500];
  float lossprob2;
  float length2;
  struct Node *neighbor3;
  ex symblossprob3;
  ex pkfh3[500]; // actually this is not used
  float lossprob3;
  float length3;
  int marker;
  float disttoroot;
  int compdone;
} Node;


typedef struct Params{
  int anzahl;
  struct Node *tree;
  float datagfs[];   /* Incomplete array type */
} Params;

typedef struct Params_symbolic{
  int anzahl;
  float datagfs[500];
  ex  symbolicgfs[500];
  symbol  rhoS;
} Params_symbolic;

typedef struct treepart{
 int anzahl;
 float treeparts[]; /* Incomplete array type */
} treepart;




#endif

// void printrootedTree(Node *node);