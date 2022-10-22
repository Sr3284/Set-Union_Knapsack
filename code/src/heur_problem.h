#ifndef __SCIP_HEUR_PROBLEM_H__
#define __SCIP_HEUR_PROBLEM_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
   int c;
   SCIP_Real solval;
} candidatoT;

int randomIntegerB (int low, int high);
int getLPsolution(SCIP* scip, SCIP_VAR** pvars, int *pn1, int *pnfrac, int *pn0, int *pnlpcands);
void printLPvars(SCIP* scip, SCIP_VAR** pvars, int n1, int nfrac, int n0);
int updateSolution(SCIP_VAR* var, instanceT* I, int* covered, int *nCovered, int *cost);
SCIP_Real createSolution(SCIP* scip, SCIP_SOL* sol, SCIP_VAR** solution, int nSolution, int *infeasible, int *covered);
int isFeasibleColumn(SCIP* scip, SCIP_VAR** solution, int nInSolution, int* covered, SCIP_VAR* var);
SCIP_RETCODE selectCand(SCIP* scip, SCIP_VAR** solution, int nInSolution, int cost, SCIP_VAR** pvar, SCIP_VAR** varlist, int n1, int nfrac, int* covered);
SCIP_RETCODE SCIPtrySolMine(SCIP* scip, SCIP_SOL* sol, SCIP_Bool printreason, SCIP_Bool checkbounds, SCIP_Bool checkintegrality, SCIP_Bool checklprows, SCIP_Bool *stored);
int isCompleteSolution(SCIP_VAR** solution, int nInSolution, int maxInSolution, int* covered, int nCovered, int n);
#ifdef __cplusplus
}
#endif

#endif
