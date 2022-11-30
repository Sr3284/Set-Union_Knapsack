/**@file   problem.c
 * @brief  This file contains routines specific for the problem and the functions loadInstance(), freeInstance, 
 * printInstance, and loadProblem must be implemented 
 *
 **/ 
#include<stdio.h>
#include<math.h>
#include "scip/scip.h"
#include "problem.h"
#include "probdata_mochila.h"

void freeInstance(instanceT** I)
{
  int i;
  
  if(*I){
     free((*I)->item);
     for(i=0;i<(*I)->n;i++)
        free((*I)->R[i]);
     free((*I)->R);
    free(*I);
    *I = NULL;
  }
}
void createInstance(instanceT** I, const char* instanceName, int n, int m, int C)
{
  int i;
  
  *I = (instanceT*) malloc(sizeof(instanceT));
  (*I)->item = (itemType*) malloc(sizeof(itemType)*n);
  (*I)->weight = (int*) malloc(sizeof(int)*m);
  (*I)->R = (int**) malloc(sizeof(int*)*n);
  for(i=0;i<n;i++){
    (*I)->R[i] = (int*) malloc(sizeof(int)*m);
  }
  (*I)->n = n;
  (*I)->m = m;
  (*I)->C = C;
  (*I)->instanceName = instanceName;
}
void printInstance(instanceT* I)
{
  int i,j;
  printf("\nInstance with n=%d items, m=%d elements", I->n, I->m);
  for(i=0;i<I->m;i++){
     printf("\nElement %d weight=%d", i+1, I->weight[i]);
  }
  printf("\nItems= \n");
  for(i=0;i<I->n;i++){
     printf("%d value=%d\nR={", I->item[i].label, I->item[i].value);
     for(j=0;j<I->m;j++){
       if(I->R[i][j]){
         printf("%d ", j);
       }
     }
     printf("}\n");
  }
}
int loadInstance(const char* filename, instanceT** I)
{
  FILE* fin;
  int n, m, i, C, j, nR;
  fin = fopen(filename, "r");
  if(!fin){
    printf("\nProblem to open file %s\n", filename);
    return 0;
  }

  fscanf(fin,"%d %d %d\n", &n, &m, &C);
  createInstance(I, filename, n, m, C);
  for(i=0; i<n; i++){
     fscanf(fin, "%d\n", &((*I)->item[i].value));
     //     printf("item: %d valor=%d\n", i, ((*I)->item[i].value));
     (*I)->item[i].label=i;
     (*I)->item[i].weight=0;
  }
  for(i=0; i<m; i++){
     fscanf(fin, "%d\n", &((*I)->weight[i]));
     //     printf("Element: %d peso=%d\n", i, ((*I)->weight[i]));
  }
  nR=0;
  for(i=0; i<n; i++){
     //    printf("Set of item %d={", i);
    for(j=0; j<m; j++){
      fscanf(fin, "%d", &((*I)->R[i][j]));
      if((*I)->R[i][j]){
         nR++;
         //         printf("%d ", j);
         (*I)->item[i].weight+=(*I)->weight[j];
      }
    }
    //    printf("}. Peso do item=%d.\n", (*I)->item[i].weight);
  }
  (*I)->nR = nR;
#ifdef DEBUG_INSTANCE
  printInstance(*I);
#endif
  fclose(fin);
  return 1;
}
// load instance problem into SCIP
int loadProblem(SCIP* scip, instanceT* I, int relaxed, int* fixed)
{
  SCIP_RETCODE ret_code;

  ret_code = SCIPprobdataCreate(scip, I->instanceName, I, relaxed, fixed);
  if(ret_code!=SCIP_OKAY)
    return 0;
  return 1;
}
