// this file contains the function that relied on the symbolic version of this program and have not been in the file treesymbolic.h

// from treestatistics.h

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


// from treestructure.h

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


typedef struct Params_symbolic{
  int anzahl;
  float datagfs[500];
  ex  symbolicgfs[500];
  symbol  rhoS;
} Params_symbolic;




// from treefunctions.h

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
  
  node->symblossprob1 = 0;
  for(inc = 0;inc < 500; inc++){
    node->pkfh1[inc] = 0;
  }
  node->lossprob1 = 0;
  node->length1 = 0;
  node->symblossprob2 = 0;
  for(inc = 0;inc < 500; inc++){
    node->pkfh2[inc] = 0;
  }
  node->lossprob2 = 0;
  node->length2 = 0;
  node->symblossprob3 = 0;
  for(inc = 0;inc < 500; inc++){
    node->pkfh3[inc] = 0;
  }
  node->lossprob3 = 0;
  node->length3 = 0;
  node->marker = 0;
  node->disttoroot = 0;
  node->compdone = 0;
    }
}




// from treeprints.h


void printsymbprobs(Node *node){
  ex pp,c1p,c2p;
  int one = 0, two = 0, three=0;
  if (node->neighbor1 == node->parent) {pp = node->symblossprob1 ; one = 1;}
  else if (node->neighbor2 == node->parent)  {pp = node->symblossprob2 ; two = 1;}
  else if (node->neighbor3 == node->parent)  {pp = node->symblossprob3 ; three = 1;}
  else {printf("WARNING: parent is not a neighbor\n"); }
  
  if (node->neighbor1 == node->child1)  {c1p = node->symblossprob1 ; one = 1;}
  else if (node->neighbor2 == node->child1)  {c1p = node->symblossprob2 ; two = 1;}
  else if (node->neighbor3 == node->child1) {c1p = node->symblossprob3 ; three = 1;}
  else {printf("WARNING: child1 is not a neighbor\n"); }
  
  if (node->neighbor1 == node->child2) {c2p = node->symblossprob1 ; one = 1;}
  else if (node->neighbor2 == node->child2) {c2p = node->symblossprob2 ; two = 1;}
  else if (node->neighbor3 == node->child2) {c2p = node->symblossprob3 ; three = 1;}
  else {printf("WARNING: child2 is not a neighbor\n"); }
  
  if (one+two+three != 3 && (node->child1) && (node->child2) ) printf("WARNING: relatives != neighbors\n");
  
  //printf("id: \t %d \t dist: %.2f \t lossp: \t %.2f  \t lossp: \t %.2f  \t lossp: \t %.2f \n", node->id, node->disttoroot , pp, c1p, c2p);
  printf("id: \t %d \t dist: %.2f \t", node->id, node->disttoroot);
  cout << " lossp: \t " << pp.evalf() << "\t lossp: \t " << c1p.evalf()  << " \t lossp: \t " << c2p.evalf()  << endl;
}


void printgfs_thin_symb(ex feld[],int anzahl){
  int i;
     for(i=0;i<anzahl;++i)
//          printf( "%4.1f\t", feld[i]*t/r);
       cout << feld[i] << " \t\t\t ";
     printf("\n");
}




// from panicmage.c

// in main function
    // this defines the global variable for the paremeter rho (the gene loss rate) 
//     symbol x("rhoS");
