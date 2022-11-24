/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_myrounding.c
 * @brief  rounding primal heuristic
 * @author Edna Hoshino (based on template provided by Tobias Achterberg)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "probdata_mochila.h"
#include "parameters_mochila.h"
#include "heur_myrounding.h"
#include "heur_problem.h"

//#define DEBUG_ROUNDING 1
/* configuracao da heuristica */
#define HEUR_NAME             "myrounding"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         'r'
#define HEUR_PRIORITY         3 /**< heuristics of high priorities are called first */
#define HEUR_FREQ             1 /**< heuristic call frequency. 1 = in all levels of the B&B tree */
#define HEUR_FREQOFS          0 /**< starts of level 0 (root node) */
#define HEUR_MAXDEPTH         10 /**< maximal level to be called. -1 = no limit */
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE /**< when the heuristic should be called? SCIP_HEURTIMING_DURINGLPLOOP or SCIP_HEURTIMING_AFTERNODE */
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#ifdef DEBUG
   #define PRINTF(...) printf(__VA_ARGS__)
#else
   #define PRINTF(...) 
#endif

/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
/*struct SCIP_HeurData
{
};
*/

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRounding)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRounding)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRounding)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRounding)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolRounding)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolRounding)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/**
 * @brief Core of the rounding heuristic: it builds one solution for the problem by rounding procedure.
 *
 * @param scip problem
 * @param sol pointer to the solution structure where the solution wil be saved
 * @param heur pointer to the rounding heuristic handle (to contabilize statistics)
 * @return int 1 if solutions is found, 0 otherwise.
 */
int rounding(SCIP* scip, SCIP_SOL** sol, SCIP_HEUR* heur)
{
   int found, infeasible, nInSolution;
   unsigned int stored;
   int nvars, n1, nfrac, n0;
   int *covered, n, custo, nCovered, nCands;
   SCIP_VAR *var, **solution, **varlist;//, **subvars;
   //  SCIP* scip_cp;
   SCIP_Real valor, bestUb;
   SCIP_PROBDATA* probdata;
   int i, error;
   instanceT* I;
   
   found = 0;
   infeasible = 0;
   
#ifdef DEBUG_ROUNDING
   printf("\n============== New rounding heur at node: %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif

   /* recover the problem data*/
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);

   nvars = SCIPprobdataGetNVars(probdata);
   I = SCIPprobdataGetInstance(probdata);
   n = I->n;
    
   varlist = (SCIP_VAR**)malloc(sizeof(SCIP_VAR*)*(nvars));
   solution = (SCIP_VAR**) malloc(sizeof(SCIP_VAR*)*n);
   covered = (int*) calloc(n,sizeof(int));
   nInSolution = 0;
   nCovered = 0;
   custo = 0;

   // get LP solution, by spliting variables in three groups according to its LP value
   getLPsolution(scip, varlist, &n1, &nfrac, &n0, NULL);
#ifdef DEBUG_ROUNDING
   printLPvars(scip, varlist, n1, nfrac, n0);
#endif
   // first, select all variables with LP value equal to 1.0
   for(i=0;i<n1;i++){
      var = varlist[i];
      solution[nInSolution++]=var;
      // update vertex covered by the current solution
      infeasible = !updateSolution(var, I, covered, &nCovered, &custo);
#ifdef DEBUG_ROUNDING
      printf("\nSelected var= %s. TotalItems=%d weight=%d infeasible=%d", SCIPvarGetName(var), nInSolution, custo, infeasible);
#endif
   }
   // complete solution using variables with fractional values
   error = 0;
   found = 0;//isCompleteSolution(solution, nInSolution, I->n, covered, nCovered, n);
   for(i=0;i<nfrac && !found && !error;i++){
      var = varlist[n1+i];
      nCands = selectCand(scip, solution, nInSolution, custo, &var, varlist, n1, nfrac, covered);
      // only select actived var in scip and whose rounding up is valid for the problem
      if(nCands>0){
         // include selected var in the solution
         solution[nInSolution++]=var;
         // update vertex covered by the current solution
         infeasible = !updateSolution(var, I, covered, &nCovered, &custo);
#ifdef DEBUG_ROUNDING
         printf("\n\nSelected var= %s. TotalItems=%d cost=%d infeasible=%d\n", SCIPvarGetName(var), nInSolution, custo, infeasible);
#endif
      }
      // check if the solution is complete
      found = isCompleteSolution(solution, nInSolution, I->n, covered, nCovered, n);
      if(!nCands){
#ifdef DEBUG_ROUNDING
         printf("Problem has no candidate!\n");
#endif
         found = nInSolution>0;
         error = !found;
      }
   }
   if(found){
      /* create SCIP solution structure sol */
      SCIP_CALL( SCIPcreateSol(scip, sol, heur) );
      // save found solution in sol
      valor = createSolution(scip, *sol, solution, nInSolution, &infeasible, covered);
      bestUb = SCIPgetPrimalbound(scip);
#ifdef DEBUG_ROUNDING
      printf("\nFound solution...\n");
      //      SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
      printf("\ninfeasible=%d value = %lf > bestUb = %lf? %d\n\n", infeasible, valor, bestUb, valor > bestUb + EPSILON);
#endif
      if(!infeasible && valor > bestUb + EPSILON){
#ifdef DEBUG_ROUNDING
         printf("\nBest solution found...\n");
         SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif
         
         /* check if the solution is feasible */
         SCIP_CALL( SCIPtrySolMine(scip, *sol, TRUE, TRUE, FALSE, TRUE, &stored) );
         if( stored )
         {
#ifdef DEBUG_PRIMAL
            printf("\nSolution is feasible and was saved! Total of items = %d", nInSolution);
            SCIPdebugMessage("found feasible rounding solution:\n");
            SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
            //       *result = SCIP_FOUNDSOL;
         }
         else{
            found = 0;
#ifdef DEBUG_ROUNDING
            printf("\nCould not found\n. BestUb=%lf", bestUb);
#endif
         }
      }
   }
#ifdef DEBUG_ROUNDING
   getchar();
#endif
   free(varlist);
   free(solution);
   free(covered);
   return found;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRounding)
{  /*lint --e{715}*/
   SCIP_SOL*             sol;                /**< solution to round */
   int nlpcands;

   assert(result != NULL);
   //   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* continue only if the LP is finished */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* continue only if the LP value is less than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;


   /* check if there exists integer variables with fractionary values in the LP */
   SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &nlpcands, NULL, NULL) );
   //Fractional implicit integer variables are stored at the positions *nlpcands to *nlpcands + *nfrac - 1
  
   /* stop if the LP solution is already integer   */
   if ( nlpcands == 0 )
     return SCIP_OKAY;

   /* solve rounding */
   if(rounding(scip, &sol, heur)){
     *result = SCIP_FOUNDSOL;
   }
   else{
     *result = SCIP_DIDNOTFIND;
#ifdef DEBUG_PRIMAL
     printf("\nRounding could not find feasible solution!");      
#endif
   }
   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the rounding_crtp primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMyRounding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create rounding primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, param.heur_freq, param.heur_freqofs,
         param.heur_maxdepth, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyRounding, heurFreeRounding, heurInitRounding, heurExitRounding, heurInitsolRounding, heurExitsolRounding, heurExecRounding,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, param.heur_round_freq, param.heur_round_freqofs,
         param.heur_round_maxdepth, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRounding, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRounding) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRounding) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRounding) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRounding) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolRounding) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolRounding) );
#endif

   /* add rounding primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
