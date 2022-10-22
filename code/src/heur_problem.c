/**@file   heur_problem.c
 * @brief  This file contains auxiliar routines used by primal heuristics (rounding and diving) and also contains 
 *         specific routines for the problem 
 *
 **/ 
#include <assert.h>

#include "probdata_mochila.h"
#include "parameters_mochila.h"
#include "heur_problem.h"

int randomIntegerB (int low, int high)
{
  int k;
  double d;

  d = (double) rand () / ((double) RAND_MAX + 1);
  k = d * (high - low + 1);
  return low + k;
}
/**
 * @brief split lambda variables in three groups according to its LP value:
         varlist[0..n1-1]     = whose LP value is equal to 1.0;
         varlist[n1..n0-1]    = fractional values; and 
         varlist[n0..nvars-1] = whose LP value is equal to 0.0
         nfrac = total of vars in the second group 
 *
 * @param scip problem
 * @param pvars vector of splitted vars
 * @param pn1 total of variables in the first group, whose LP values are equal to 1.0 
 * @param pnfrac total of variables in the second group
 * @param pn0 index in pvar for the first variable in the third group
 * @param pnlpcands total of compact variables that have fractionary LP value 
 * @return int 1 if LP is avaliable and variables are splitted, 0 otherwise.
 */
int getLPsolution(SCIP* scip, SCIP_VAR** pvars, int *pn1, int *pnfrac, int *pn0, int *pnlpcands)
{
   int i, nvars, n0, n1, nfrac, nlpcands;
   SCIP_PROBDATA* probdata;
   SCIP_VAR **vars, *var;
   SCIP_Real solval;

   nlpcands = 0;

   if(pvars==NULL)
      return 0;

   /* continua somente se LP concluido */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return 0;
   /* recupera os dados do problema */
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);

   vars = SCIPprobdataGetVars(probdata);
   nvars = SCIPprobdataGetNVars(probdata);
   n1 = 0;
   n0 = nvars;
   nfrac = 0;
   for(i=0;i<nvars;i++){
      var = vars[i];
      solval = SCIPgetVarSol(scip,var);
      if(solval> 1 - EPSILON){
         pvars[n1+nfrac]=pvars[n1];
         pvars[n1++]=var;
      }
      else if(solval < EPSILON){
         pvars[--n0]=var;
      }
      else{
         pvars[n1+nfrac]=var;
         nfrac++;
      }
   }
   *pn1=n1;
   *pn0=nvars-n0;
   *pnfrac = nfrac;
   if(pnlpcands!=NULL)
      (*pnlpcands) = nlpcands;
   return 1;
}
/**
 * @brief Print in the standard output the list of variables in a specific order as obtained by getLPsolution:
         varlist[0..n1-1]     = whose LP value is equal to 1.0;
         varlist[n1..n0-1]    = fractional values; and 
         varlist[n0..nvars-1] = whose LP value is equal to 0.0
         nfrac = total of vars in the second group 
 *
 * @param scip problem
 * @param vars vector of splitted vars
 * @param n1 total of variables in the first group, whose LP values are equal to 1.0 
 * @param nfrac total of variables in the second group
 * @param n0 index in vars for the first variable in the third group
 */
void printLPvars(SCIP* scip, SCIP_VAR** pvars, int n1, int nfrac, int n0)
{
   int i;
   SCIP_Real solval;
   
   printf("\nTotal of vars with value 1.0 = %d\n", n1);
   for(i=0;i<n1;i++){
      solval = SCIPgetVarSol(scip,pvars[i]);
      printf("%s (%lf)\n", SCIPvarGetName(pvars[i]), solval);
   }
   printf("\nTotal of vars comwith value fractional = %d\n", nfrac);
   for(i=n1;i<n1+nfrac;i++){
      solval = SCIPgetVarSol(scip,pvars[i]);
      printf("%s (%lf)\n", SCIPvarGetName(pvars[i]), solval);
   }
   printf("\nTotal of vars with valur 0.0 = %d\n", n0);
#ifdef DEBUG
   for(i=n1+nfrac;i<n1+nfrac+n0;i++){
      solval = SCIPgetVarSol(scip,pvars[i]);
      printf("%s (%lf)\n", SCIPvarGetName(pvars[i]), solval);
   }
#endif
}
/* TODO: it should contain codes that depend on the problem */
/**
 * @brief Update the cost and the vector of itens covered
 *
 * @param var the variable that is used to update the current solution (given by covered and cost)
 * @param n the size of itens, which refers to the part of consids's var that has information concerning itens
 * @param covered the vector of itens already covered (and that must be updated) 
 * @param nCovered the number of itens in covered
 * @param cost the current cost of the solution, that must be updated
 * @return always return 1 (not used yet). 
 */
int updateSolution(SCIP_VAR* var, instanceT *I, int* covered, int *nCovered, int *cost)
{
   int a, i;
// update item covered by the current solution
   const char* name;

   name = SCIPvarGetName(var);
   a = 0;
   for(i=4;name[i]!='\0';i++){
      a = a*10 + name[i]-48;
   }
   covered[a]=1;
   (*nCovered)++;
   *cost += I->item[a].weight;
   return 1;
}
// TODO: create solution based on vars in solution (specific for the problem)
SCIP_Real createSolution(SCIP* scip, SCIP_SOL* sol, SCIP_VAR** solution, int nSolution, int *infeasible, int *covered)
{
   int s, i;
   SCIP_VAR *var;
   SCIP_PROBDATA* probdata;
   SCIP_Real value;
   int weight; 
   instanceT* I;

  /* recover the data problem */
  probdata=SCIPgetProbData(scip);
  assert(probdata != NULL);

  I = SCIPprobdataGetInstance(probdata);

  value = 0;
  weight = 0;
  for(i=0;i<I->n;i++)
  {
     if(covered[i]){
        value += I->item[i].value;
        weight += I->item[i].weight;
     }
  }
  for(s=0;s<nSolution;s++){
    var = solution[s];
    SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );

#ifdef DEBUG_PRIMAL
    printf("\nAdd %s in sol.", SCIPvarGetName(var));
#endif
  } // for each var at solution

  *infeasible = 0;
  if(weight > I->C){
     *infeasible = 1;
  }

  if(*infeasible){
    printf("\nSolution has became infeasible!!\n");
  }
  return value;
}

/* TODO: it depends on the problem. */
int isFeasibleColumn(SCIP* scip, SCIP_VAR** solution, int nInSolution, int* covered, SCIP_VAR* var)
{
  SCIP_PROBDATA* probdata;
  int i, feasible, weight;
  instanceT* I;
  int a;
// update item covered by the current solution
  const char* name;

  name = SCIPvarGetName(var);
  a = 0;
  for(i=4;name[i]!='\0';i++){
     a = a*10 + name[i]-48;
  }
  
  probdata=SCIPgetProbData(scip);
  assert(probdata != NULL);
  I = SCIPprobdataGetInstance(probdata);
  
  feasible = 1;
  /* recover problem's data */
  probdata=SCIPgetProbData(scip);
  assert(probdata != NULL);
  // TODO: we are considering that there is one constraint for each compact variable in PMR
  // TODO: should consider specific requirements for the problem
  weight = 0;
  for(i=0;i<I->n && weight <= I->C; i++){
     if(covered[i]){
        weight+=I->item[i].weight;
     }
  }
  if(weight + I->item[a].weight > I->C){ // Here, we just do not allow non-disjoint columns and non-relaxed columns
    feasible = 0;
  }
  return feasible;
}
/* TODO: can be improved, based on specific properties of the problem */
SCIP_RETCODE selectCand(SCIP* scip, SCIP_VAR** solution, int nInSolution, int cost, SCIP_VAR** pvar, SCIP_VAR** varlist, int n1, int nfrac, int* covered)
{
// only select actived var in scip_cp and whose rounding up is valid for the problem
   int c; // v, m, nvars
   SCIP_VAR* var, *bestSol;//, **vars;
   SCIP_PROBDATA* probdata;
   SCIP_Real solval, bestSolVal;
   int found;
   instanceT* I;
   
   /* recupera os dados do problema */
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);
   I = SCIPprobdataGetInstance(probdata);
   
   found = 0;
   bestSol = NULL;
   if(!found){
      // select the best fractionary variable 
      bestSolVal = 0;
      for (c = n1; c < n1 + nfrac; c++)
      {
         var = varlist[c];
         solval = SCIPgetVarSol(scip,var);    
         if(isFeasibleColumn(scip,solution, nInSolution, covered, var)){
#ifdef DEBUG_ROUNDING
            printf("\nVar %s = %lf", SCIPvarGetName(var), SCIPgetVarSol(scip,var));
            printf("... feasible");
#endif
            if(bestSolVal<solval){
               bestSolVal = solval;
               bestSol = var;
               found = 1;
            }
         }
      }
   }
   // if there is no fractionary variable, select the first one not used feasible variable.
   for (c = n1 + nfrac; !found && c < I->n; c++)
   {
      var = varlist[c];
      if(isFeasibleColumn(scip,solution, nInSolution, covered, var)){
#ifdef DEBUG_ROUNDING
         printf("\nVar %s = %lf", SCIPvarGetName(var), SCIPgetVarSol(scip,var));
         printf("... feasible");
#endif
         found = 1;
         bestSol = var;
      }
   }
   *pvar = bestSol;
#ifdef DEBUG_ROUNDING
   if(found)
      printf("\nSelect var %s = %lf", SCIPvarGetName(bestSol), bestSolVal);
#endif
   return found;
}

SCIP_RETCODE SCIPtrySolMine(SCIP* scip, SCIP_SOL* sol, SCIP_Bool printreason, SCIP_Bool checkbounds, SCIP_Bool checkintegrality, SCIP_Bool checklprows, SCIP_Bool *stored)
{
#if (defined SCIP_VERSION_MAJOR)  ||  (SCIP_VERSION_MAJOR==6)
   return SCIPtrySol(scip,sol, TRUE,printreason,checkbounds,checkintegrality,checklprows,stored);
#else
   return SCIPtrySol(scip,sol,printreason,checkbounds,checkintegrality,checklprows,stored);
#endif
}

/* TODO: it depends on the problem */
int isCompleteSolution(SCIP_VAR** solution, int nInSolution, int maxInSolution, int* covered, int nCovered, int n)
{
   return 1;//nInSolution>=maxInSolution || nCovered==n;
}
