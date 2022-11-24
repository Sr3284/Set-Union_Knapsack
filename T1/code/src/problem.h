#ifndef __PROBLEM__
#define __PROBLEM__
#include<stdio.h>
#include "scip/scip.h"

/** structure for each item */
typedef struct{
  int label; 
  int value;
  int weight;
}itemType;

typedef struct{
   int n;
   int m;
   int C;
   int *weight;
   int **R;
   int nR;
   itemType *item; /**< data for each item in 0..n-1 */
} instanceT;

void freeInstance(instanceT* I);
void createInstance(instanceT** I, int n, int m, int C);
void printInstance(instanceT* I);
// load instance from a file
int loadInstance(char* filename, instanceT** I);
// load instance problem into SCIP
int loadProblem(SCIP* scip, char* probname, instanceT* in);
#endif
