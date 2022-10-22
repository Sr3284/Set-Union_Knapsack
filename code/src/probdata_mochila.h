/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_mochila.h
 * @brief  Problem data for mochila problem
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_MOCHILA__
#define __SCIP_PROBDATA_MOCHILA__

#include "scip/scip.h"
#include "problem.h"

/* constants */

/* macros */
#define EPSILON 0.000001
#ifdef DEBUG
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...) 
#endif

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the mochila, all variables which are created initially
 */
struct SCIP_ProbData
{
   const char*           probname;           /**< problem name */
   SCIP_VAR**            vars;               /**< array of variables */
   SCIP_CONS**           conss;              /**< all constraints */
   int                   nvars;              /**< total of vars */
   int                   ncons;              /**< number of constraints */
   instanceT*            I;                  /**< instance of knapsack */
};

/** sets up the problem data */
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   instanceT*            I                   /**< instance of K-coloring */
   );

/** adds given variable to the problem data */
extern
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var                 /**< variables to add */
   );

/** returns Probname of the instance */
extern
const char* SCIPprobdataGetProbname(
   SCIP_PROBDATA*        probdata            /**< problem data */
			      );
/** returns array of all variables ordered in the way they got generated */
extern
SCIP_VAR** SCIPprobdataGetVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of variables */
extern
int SCIPprobdataGetNVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of set partitioning constrains */
extern
SCIP_CONS** SCIPprobdataGetConss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of set partitioning constrains */
extern
int SCIPprobdataGetNcons(
   SCIP_PROBDATA*        probdata            /**< problem data */
			 );

/** returns instance I */
extern
instanceT* SCIPprobdataGetInstance(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );
#endif
