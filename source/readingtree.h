/***** readingtree.h ***********************************************
 * Description: reading the newick tree and transforming this tree to custom data format
 * Author: Franz Baumdicker, baumdicker@stochastik.uni-freiburg.de ; Peter Pfaffelhuber, pp@stochastik.uni-freiburg.de
 * File created 2010.
 *****************************************************************/

#include "treestructure.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>






int split_childs_info(char *s1,char *childs, char *idlength){
  int leaveornot;
  char *ptr;
  
  ptr = strrchr(s1, ')');

  if (ptr == NULL){
    childs[0] = '\0';
    strcpy(idlength,s1);
  }
  else{
  *ptr = '\0';
  strcpy(childs,s1+1);
//   printf("%s\n" , childs);
  strcpy(idlength, ptr+1);
  }  
  if ( (childs[0] == '\0') || ( childs[0] == '(' && childs[1] == ')' && childs[2] == '\0' ) ){
//     printf("empty string or '()' found => node is leave\n");
    leaveornot = 1;
  }else{
    leaveornot = 0;
  }
  return leaveornot;
}





int assignnode(int leaveornot, int n){
    static int innerid;
    static int leaveid;
    if (leaveornot == 0){ // an inner node
      innerid++;
      return n-1-(innerid-1);
    }
    if (leaveornot == 1){ // a leave
      leaveid++;
      return leaveid-1;
    }
}
 
 
  
float readidlength(char *idlength_input, int place, int parent, Node *tree, int sec_child){
      //static int leaveid;
      //static int innerid;
      char *ptr;
      
      char idlength[50000], str1[50000], str2[50000];
      strcpy(idlength, idlength_input);
      // start reading idlength
      char trennzeichen[] = ":";
      int counter = 0;
      int ID = -50;
      float length = -2.;
      ptr = strtok(idlength, trennzeichen);
      if(ptr == NULL) {  // thus it is an empty string
	// search ID
	ID = place;
	length = 0.;
      }
      if(ptr != NULL){
// 	printf("%s\n", ptr);
	// warning if ptr is not of type float not included
	strcpy(str1,ptr);
	counter++;
	ptr = strtok(NULL, trennzeichen);
      }
      if(ptr != NULL){
// 	printf("%s\n", ptr);
	// warning if ptr is not of type float not included
	strcpy(str2,ptr);
	counter++;
	ptr = strtok(NULL, trennzeichen);
      }
      if( ptr != NULL){
// 	printf("%s\n", ptr);
	counter++;
	printf("WARNING: too many ':' in newick file additional information is ignored\n");
      }
//       printf("Counter: %d\n", counter);
      if (counter == 1){
	if( idlength_input[0] == ':'){
	  length = atof(str1);
	  // a ID has to be set to some integer
	  ID = place;
	}
	else{
	ID = atoi(str1);
	length = 0.;
	}
      }
      else if (counter >= 2){
	ID = atoi(str1);
	length = atof(str2);
      }
    
//     printf("ID: %d \t Length: %.3f to parent %d \t\t of node %d \n", ID, length, parent, place);
//     printf("btw sec_child is: %d\n", sec_child);
    (tree+place)->id = ID;
    (tree+place)->length3 = length;
    if (parent == -1){ // then place is root
      (tree+place)->neighbor3 = NULL;
      (tree+place)->parent = NULL;
    }
    else{  
      (tree+place)->neighbor3 = tree+parent;
      (tree+place)->parent = tree+parent;
    }
    if (parent >= 0 && sec_child == 0){
      (tree+parent)->child1 = tree+place;
      (tree+parent)->length1 = length;
      (tree+parent)->neighbor1 = tree+place;
    }
    if (parent >= 0 && sec_child == 1){
      (tree+parent)->child2 = tree+place;
      (tree+parent)->length2 = length;
      (tree+parent)->neighbor2 = tree+place;
    }   
    return length;
}
  
  
  
  
  
int splitchilds(char *childs, char *child1, char *child2){
  if (childs[0] == '\0') {  
  }
  // Start split childs
  else{    
    int dex = 0, count = 0, marker = 0, numchilds = 1;
    char *ch2ptr;
    char save;
    int pos;
    
    while (childs[dex] != '\0'){
//        printf("%d\n",dex);
       if (childs[dex] == '(') count++;
       if (childs[dex] == ')') count--;
       if (childs[dex] == ',' && count == 0 && numchilds >= 2) { 
	 printf("ERROR: node with more than three neighbours found\n");
	 abort;
       }
       if (childs[dex] == ',' && count == 0){
// 	 printf("second child found at %d \n", dex);
	 save = childs[dex];
	 pos = dex;
	 childs[dex] = '\0';
	 numchilds++;
	 marker = 1;
	 ch2ptr = childs+dex+1;
       }
       dex++;
    }
    strcpy(child1, childs);
    if (numchilds == 2) strcpy(child2, ch2ptr);
    childs[pos] = save;
    return numchilds;
  }
// End split childs
}




void parsenewick(char *s1 , int parent , int n, Node *tree, int sec_child){
int leaveornot = 0;
char childs[50000], idlength[50000], child1[50000], child2[50000];
int place,numchilds;
float nodelength;

leaveornot = split_childs_info(s1,childs,idlength);
place = assignnode(leaveornot,n);

if (sec_child == -1){  // then place is the root
  (tree+place)->parent = NULL;
  (tree+place)->neighbor3 = NULL;
  (tree+place)->length3 = 0.;
}
nodelength = readidlength(idlength,place,parent,tree,sec_child);
numchilds = splitchilds(childs,child1,child2);


// printf("parent:\t %s\n", childs);
// printf("child1:\t %s\n" , child1);
// if (numchilds == 2) printf("child2:\t %s\n\n" , child2);
// else printf("child2: \t xxxxx\n\n");

  if (leaveornot == 0){
    parsenewick(child1,place,n,tree,0);
    if (numchilds == 2) parsenewick(child2,place,n,tree,1);
    else{
      (tree+place)->child2 = NULL;
      (tree+place)->neighbor2 = NULL;
      (tree+place)->length2 = 0.;
    }
  }

  if (leaveornot == 1 && sec_child == 0){
    (tree+parent)->child1 = tree+place;
    (tree+parent)->neighbor1 = tree+place;
    (tree+parent)->length1 = nodelength;
  }
  if (leaveornot == 1 && sec_child == 1){
    (tree+parent)->child2 = tree+place;
    (tree+parent)->neighbor2 = tree+place;
    (tree+parent)->length2 = nodelength;
  }

}
