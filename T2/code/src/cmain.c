/* cmain.c
Heuristicas baseadas em modelos matematicos para resolver o problema da mochila
usando a biblioteca SCIP para resolver o PL ou o MIP (mixed integer programming)
By: Edna A. Hoshino 

Problema da mochila: 
dados um conjunto de itens I={1,2,...,n} e uma mochila K
em que cada item i tem um peso pi e um valor vi
e a mochila tem uma capacidade C,
encontrar um subconjunto dos itens que podem ser transportados na mochila de modo
a maximizar o valor transportado. Ou seja, a soma do valor dos itens escolhidos
deve ser maxima e a soma dos pesos dos itens escolhidos nao pode
ultrapassar a capacidade da mochila.

Modelo de PLI para o problema da mochila (KP):

max v1x1 + v2x2 + ... vnxn
s.a.:
   p1x1 + p2x2 + ... pnxn <= C
   
variaveis: xi binarias = 1 se o item i eh transportado na mochila

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "probdata_mochila.h"
#include "problem.h"

#define EPSILON 0.000001

#ifdef DEBUG
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...)
#endif

typedef struct{
   // global settings
   int time_limit; /* limit of execution time (in sec). Default = 1800 (-1: unlimited) */
   int display_freq; /* frequency to display information about B&B enumeration. Default = 50 (-1: never) */
   int nodes_limit; /* limit of nodes to B&B procedure. Default = -1: unlimited (1: onlyrootnode) */
   char* parameter_stamp;
} parametersT;

double getDualBound(instanceT I, double *tempo_gasto);
int RandomInteger(int low, int high);
int comparador(const void *valor1, const void *valor2);
double heuristica_LNS(instanceT I, double *x, double z, int tempo, double porcDestroy);
double heuristica(instanceT I, int tipo, double *x, int tempo, float tam, int cLNS, double porcDestroy, double* z0, double* tempoLNS);
int heuristica_Gulosa (instanceT I, double *x, double *aux, int tempo);
void heuristica_Aleatoria (instanceT I, double *x, double *aux);
int heuristica_relax_fix(instanceT I, double* x, float percentualTam, int tempo, double *aux);
// altera cabeçalho da otimiza_PLI para usar a variavel info (necessária para a callback)
double otimiza_PLI(instanceT I, int tipo, double *x, parametersT param, double* finalDualBound, int *nodes, int *ativos);
// cabeçalho da callback
int validarSolucao(instanceT I, double *x, double z, int tipo);
void salva_solucao(char *name, int tipo, double z, double *x, instanceT* I, int validacao, int tempo, int cLNS, double porcDestroy, float tam);
void printSol(SCIP* scip, char* outputname, parametersT param);
SCIP_RETCODE configScip(SCIP** pscip, parametersT param);
void  buildEdges(instanceT* I, instanceT* II);

/** 
 * creates a SCIP instance with default plugins, and set SCIP parameters 
 */
SCIP_RETCODE configScip(
  SCIP** pscip,
  parametersT param
  )
{
   SCIP* scip = NULL;
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) ); 
   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   /* for column generation, disable restarts */
   //   SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) ); 
   /* disable presolving */
   //   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) ); // turn off
   /* turn off all separation algorithms */
   //   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) ); // turn off
   /* disable heuristics */
   //   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );  // turn off
   /* for column generation, usualy we prefer branching using pscost instead of relcost  */
   //   SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", 1000000) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", param.display_freq) );
   /* set time limit */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", param.time_limit) );
   // for only root, use 1
   SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", param.nodes_limit) );
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );   
   *pscip = scip;
   return SCIP_OKAY;
}

// salva a solucao em arquivo
void salva_solucao(char *name, int tipo, double z, double *x, instanceT* I, int validacao, int tempo, int cLNS, double porcDestroy, float tam)
{
  int i, itens;
  char filename[500];
  char *instance_filename;
  FILE *fout;
  char diretorio[100] = "";
  
  instance_filename = strrchr(name, '/');
  instance_filename++;
  strcat(diretorio, instance_filename);

  sprintf(filename, "%s-%d-%d-%d-%.0lf-(%f).sol", diretorio, tipo, tempo, cLNS, porcDestroy, tam);
  if (!(fout = fopen(filename, "w")))
  {
    printf("\nProblema para salvar solucao em %s\n", filename);
    return;
  }

  itens = 0;
  for (i = 0; i < I->n; i++)
  {
    if (x[i] > EPSILON)
    {
      itens++;
    }
  }
  fprintf(fout, "%.0lf %d\n", z, itens);
  for (i = 0; i < I->n; i++)
  {
    if (x[i] > EPSILON)
    {
      fprintf(fout, "%d (v=%d p=%d)\n", i + 1, I->item[i].value, I->item[i].weight);
    }
  }
  
  if (validacao != 1)
  {
    switch (validacao)
    {
      case -1:
        fprintf(fout, "Solução incorreta! A soma dos pesos dos itens selecionados excede a capacidade da mochila!\n");
        break;
      case -2:
        fprintf(fout, "Solução incorreta! A soma dos valores dos itens selecionados está diferente do valor encontrado pela heurística!\n");
        break;
    }
  }
  
  fclose(fout);
}


/* sorteia um numero aleatorio entre [low,high] */
int RandomInteger(int low, int high)
{
  int k;
  double d;
  //  srand(time(NULL));

  d = (double)rand() / ((double)RAND_MAX + 1);
  k = d * (high - low + 1);
  return low + k;
}


// TODO: Get the best solution found and write the solution in a file. It depends on the problem!
void printSol(SCIP* scip, char* outputname, parametersT param)
{
   SCIP_PROBDATA* probdata;
   SCIP_SOL* bestSolution;
   SCIP_VAR** vars;
   SCIP_Real solval;
   FILE *file;
   int v, nvars;
   instanceT* I;
   char filename[SCIP_MAXSTRLEN];
   struct tm * ct;
   const time_t t = time(NULL);

   assert(scip != NULL);
   bestSolution = SCIPgetBestSol(scip);
   if( bestSolution == NULL )
     return;
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   I = SCIPprobdataGetInstance(probdata);
   nvars = SCIPprobdataGetNVars(probdata);
   vars = SCIPprobdataGetVars(probdata);

   (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s.sol", outputname);
   file = fopen(filename, "w");
   if(!file)
     {
       printf("\nProblem to create solution file: %s", filename);
       return;
     } 
   fprintf(file, "\nValue: %lf\nItems: ", SCIPsolGetOrigObj(bestSolution));

   for( v=0; v< nvars; v++ )
     {
       solval = SCIPgetSolVal(scip, bestSolution, vars[v]);
       if( solval > EPSILON )
	 {
	   fprintf(file, "%d ", I->item[v].label);
	 }
     }

   fprintf(file, "\n");
   //
   fprintf(file, "Parameters settings file=%s\n", param.parameter_stamp);
   fprintf(file, "Instance file=%s\n", SCIPgetProbName(scip));
   ct = localtime(&t);
   fprintf(file, "Date=%d-%.2d-%.2d\nTime=%.2d:%.2d:%.2d\n", ct->tm_year+1900, ct->tm_mon, ct->tm_mday, ct->tm_hour, ct->tm_min, ct->tm_sec);
   fclose(file);
}
/* resolve o problema de PLI usando o SCIP
   parameter info = guarda informações de execução do B&B
*/
double otimiza_PLI(instanceT I, int tipo, double *x, parametersT param, double* finalDualBound, int *nodes, int *ativos)
{
  SCIP* scip;
  SCIP_PROBDATA* probdata;
  double z, valor;
  int i;
#ifdef DEBUG
  int status;
  FILE* fsol;
#endif
  SCIP_VAR** vars;
  SCIP_SOL* bestSolution;
  //  clock_t start, end;
  
  // create scip and set scip configurations
  configScip(&scip, param);

  // carga do lp
  // load problem into scip
  if(!loadProblem(scip,&I,tipo==1?1:0, NULL)){
    printf("\nProblem to load instance problem\n");
    return -SCIPinfinity(scip);
  }

  // solve scip problem
  //  start=clock();
  SCIP_CALL( SCIPsolve(scip) );
  //  end = clock();

#ifdef DEBUG
  status = SCIPgetStatus(scip);
  printf("\nstatus=%d", status);
#endif
  
  probdata = SCIPgetProbData(scip);
  assert(probdata != NULL);

   // Recupera solucao
  bestSolution = SCIPgetBestSol(scip);
  *nodes = SCIPgetNTotalNodes(scip);
  *ativos = SCIPgetNNodesLeft(scip);
  *finalDualBound = SCIPgetDualbound(scip);
  if(bestSolution==NULL){
     PRINTF("\nNo solution found\n");
     return -SCIPinfinity(scip);
  }
  z = SCIPgetPrimalbound(scip);
  vars = SCIPprobdataGetVars(probdata);

  for (i = 0; i < I.n; i++)
  {
    valor = SCIPgetSolVal(scip, bestSolution, vars[i]);
    if (valor > EPSILON){
      PRINTF("x%d = %.2lf\n", I.item[i].label, valor);
    }
    x[i] = valor;
  }

#ifdef DEBUG
  // Grava solucao e PL
  PRINTF("\n---LP gravado em mochila.lp e solucao em mochila.sol");
  // print problem
  SCIP_CALL( SCIPwriteOrigProblem(scip, "mochila.lp", "lp", FALSE) );
  fsol = fopen("mochila.sol", "w");
  if(fsol!=NULL){
    SCIP_CALL( SCIPprintBestSol(scip, fsol, FALSE) );
    fclose(fsol);
  }
#endif
  // Destroi problema
  SCIP_CALL( SCIPfree(&scip) ); 
  return z;
}


/*
 tam: percentual do total de itens. Valor entre (0,1).
 */
int heuristica_relax_fix(instanceT I, double* x, float percentualTam, int tempo, double *aux)
{
  SCIP* scip;
  SCIP_PROBDATA* probdata;
  parametersT param;
#ifdef DEBUG
  int status;
  FILE* fsol;
#endif
  SCIP_VAR** vars;
  SCIP_SOL* bestSolution;
  //  clock_t start, end;
  unsigned int infeasible;
  double valor, zheur, z;
  int i, k, *particao, K, parte, frac, *fixed, tam;
#ifdef DEBUG
  char destaque='*', branco=' ';
#endif
  itemType *cand;
  int nCand, capacRes;
  int *selecionado;
#ifdef TESTE
  int ii;
#endif
  
  // aloca marcador de itens selecionados
  selecionado = (int*) calloc(I.n, sizeof(int));
  // aloca candidatos
  cand = (itemType*)malloc(sizeof(itemType)*I.n);
  fixed = (int*)calloc(I.n, sizeof(int)); // fixed[i]=0, if item i is not fixed, fixed[i]=1 if item i is fixed in 1.0, fixed[i]=-1 if item i is fixed in 0.

  // inicializa lista de candidatos
  nCand = 0;
  for(i=0;i<I.n;i++){
    cand[nCand].label = I.item[i].label;
    cand[nCand++].value = I.item[i].value;
  }

  qsort(cand, nCand, sizeof(itemType), comparador);

  // define tamanho de cada parte da particao
  tam = ceill(percentualTam * I.n);
  K = ceill(1/percentualTam);
  // aloca vetor das particoes
  particao = (int*)malloc(sizeof(int)*(I.n));
  for(parte=0;parte<K;parte++){
     for(i=0;i<tam && nCand>0;i++){
        k = RandomInteger(0,nCand-1); // sorteia um candidato
        particao[cand[k].label] = parte;
        // remove candidato
        cand[k] = cand[--nCand];
     }
  }

#ifdef DEBUG
  for(i=0;i<I.n;i++){
    printf("\ni=%d parte=%d", i, particao[i]);
  }
#endif

  // create scip and set scip configurations
  param.time_limit = tempo;
  param.display_freq = -1;
  param.nodes_limit = -1;
  configScip(&scip, param);
  
  frac = 1;
  capacRes = I.C;
  for(parte=0;parte<K && frac;parte++){ // itera para cada parte da particao, mas pode parar antes se a solucao ja eh inteira.
    PRINTF("\n=====parte=%d K=%d\n", parte, K);
    // carga do lp
    // load problem into scip
    if(!loadProblem(scip,&I, 1, fixed)){ // relaxation
      printf("\nProblem to load instance problem\n");
      return 0;
    }
    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);
    vars = SCIPprobdataGetVars(probdata);

    for(i=0;i<I.n;i++){
      if(particao[i]==parte){ // torna as variaveis da parte atual como binarias
        //          PRINTF("\nvar x%d binario", i+1);
         SCIPchgVarType(scip, vars[i], SCIP_VARTYPE_BINARY, &infeasible);
      }
    }
    // Executa Solver de PL
    //    start=clock();
    SCIP_CALL( SCIPsolve(scip) );
    //    end = clock();
  
#ifdef DEBUG
    status = SCIPgetStatus(scip);
    printf("\nstatus=%d", status);
#endif
     
    // Recupera solucao
    bestSolution = SCIPgetBestSol(scip);
    z = SCIPgetPrimalbound(scip);//SCIPsolGetOrigObj(bestSolution);

#ifdef DEBUG
    // Grava solucao e PL
    PRINTF("\n---LP gravado em relax_fix.lp e solucao em mochila.sol");
    // print problem
    SCIP_CALL( SCIPwriteOrigProblem(scip, "relaxFix.lp", "lp", FALSE) );
    fsol = fopen("mochila.sol", "w");
    if(fsol!=NULL){
       SCIP_CALL( SCIPprintBestSol(scip, fsol, FALSE) );
       fclose(fsol);
    }
#endif

    PRINTF("\n*****\n\t z=%lf\n", z);
    frac = 0;
    for (i = 0; i < I.n; i++)
    {
       valor = SCIPgetSolVal(scip, bestSolution, vars[i]); // recupera o valor da variavel xi
       if(valor>EPSILON){
          PRINTF("x%2d = %.2lf (parte = %2d)%c%c\n", i+1, valor, particao[i], particao[i]<=parte?destaque:branco, particao[i]==parte?destaque:branco);
          if(valor < 1.0 - EPSILON){
             frac = 1;
          }
       }
       if(particao[i]==parte){ // fixa variaveis da parte atual
          x[i]= valor>EPSILON? 1.0:0.0;
          PRINTF("x%2d = %.2lf (parte = %2d)%c%c\n", i+1, valor, particao[i], particao[i]<=parte?destaque:branco, particao[i]==parte?destaque:branco);
          if(x[i]>EPSILON){
            fixed[i]=1;
            //             SCIP_CALL( SCIPchgVarLb(scip, vars[i], x[i]) );
          }
          else{
            fixed[i]=-1;
            //             SCIP_CALL( SCIPchgVarUb(scip, vars[i], x[i]) );
          }
          if(valor>EPSILON){
             selecionado[i]++;
             capacRes -= I.item[i].weight;
          }
       }
    }
    //    getchar();
  } // next part

  zheur = 0;
  if(!frac){
    PRINTF("\nSol heur relax fix = %lf\n", z);
    zheur = z;
    if(parte<K){
      // falta copiar a solução das demais variaveis
       for (i = 0; i < I.n; i++)
       {
          valor = SCIPgetSolVal(scip, bestSolution, vars[i]); // recupera o valor da variavel xi
          if(valor>EPSILON && particao[i]>=parte){
             x[i]= valor>EPSILON? 1.0:0.0;
          }
       }
    }
#ifdef TESTE
#ifdef DEBUG
    printf("\ncapacRes=%d", capacRes);

    for(i=0;i<I.n;i++){
      if(!selecionado[i]){
        printf("\nitem %d peso=%d", i+1, I.item[i].weight);
      }
    }
#endif
    // repescagem
    for(i=0;i<I.n;i++){
      ii = cand[i].label - 1;
      if(!selecionado[ii]){
        if(capacRes>=I.item[ii].weight){
          selecionado[ii]=1;
          zheur += I.item[ii].value;
          x[ii] = 1;
          capacRes = capacRes - I.item[ii].weight;
          PRINTF("\n*Novo item=%d de peso=%d na mochila c/ peso residual=%d valor=%lf", ii+1, I.item[ii].weight, capacRes, zheur);
        }
      }
    }
#endif
  }
  
  for (i = 0; i < I.n; i++)
  {
    if (selecionado[i] > EPSILON)
    {
      aux[1] += 1.0;
    }
  }
  aux[0] = zheur;
  
  // Destroi problema
  free(particao);
  free(cand);
  free(selecionado);
  free(fixed);
  // clear problem
  //  SCIP_CALL( SCIPfree(&scip) );
  return 0;
}

double heuristica_LNS(instanceT I, double *x, double z, int tempo, double porcDestroy)
{
  SCIP* scip;
  SCIP_PROBDATA* probdata;
#ifdef DEBUG
  int status;
#endif
  SCIP_VAR** vars;
  SCIP_SOL* bestSolution;
  parametersT param;
  
  itemType *sair;
  int capacRes, nSair;
  double valor;
  int i, ii, saiu, totalSaiu, totalSelecionados;
  double zz, perda, vv;

  // aloca candidatos a sair da solução inicial
  sair = (itemType*)malloc(sizeof(itemType)*I.n);

  // inicializa capacidade residual da mochila
  capacRes = I.C;

  // Verifica quais itens foram selecionados e atualiza capacidade residual
  totalSelecionados = 0;
  for(i=0;i<I.n;i++){
    if(x[i] > EPSILON){
      capacRes -= I.item[i].weight;
      totalSelecionados++;
    }
  }

  // Decide quem sairá da solução com base no peso de cada item
  perda = 0;
  totalSaiu = 0;
  nSair = 0;

  for(i=0;i<I.n;i++){
    if(x[i]>EPSILON){
      sair[nSair].label = i;
      sair[nSair++].value = -(I.item[i].weight);
    }
  }
  
  qsort(sair, nSair, sizeof(itemType), comparador); // Ordena o vetor com base no peso dos itens selecionados
  saiu= nSair*(porcDestroy/100.0); // Destroi parte da solução
  PRINTF("\nporcDestroy=%lf saiu=%d nSair=%d\n", porcDestroy, saiu, nSair);
  for(i=0;i<saiu;i++){
    ii = sair[i].label;
    capacRes += I.item[ii].weight;
    perda += I.item[ii].value;
    PRINTF("\nRemove %d (peso=%d valor=%d) da mochila (capac residual=%d)", ii, I.item[ii].weight, I.item[ii].value, capacRes);
    x[ii]=0;
    totalSaiu++;
  }

#ifdef DEBUG
  printf("\nDestrui %.2lf%% (equivalente a um valor = %lf)\n", (100.0*totalSaiu/totalSelecionados), perda);
#endif

  // create scip and set scip configurations
  param.time_limit = tempo;
  param.display_freq = -1;
  param.nodes_limit = -1;
  configScip(&scip, param);
  
  // carga do lp
  // load problem into scip
  if(!loadProblem(scip,&I, 0, NULL)){
     printf("\nProblem to load instance problem\n");
     return -1;
  }
  probdata = SCIPgetProbData(scip);
  assert(probdata != NULL);
  // Recupera solucao
  vars = SCIPprobdataGetVars(probdata);
  // Fixa os itens após o destroy
  for (i = 0; i < I.n; i++)
  {
    if (x[i] > EPSILON)
       SCIP_CALL( SCIPchgVarLb(scip, vars[i], x[i]) );
  }  

  // print problem
  SCIP_CALL( SCIPwriteOrigProblem(scip, "lns.lp", "lp", FALSE) );
     
  // solve scip problem
  SCIP_CALL( SCIPsolve(scip) );
#ifdef DEBUG
  SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
#endif
  
#ifdef DEBUG
  status = SCIPgetStatus(scip);
  printf("\nstatus=%d", status);
#endif
     
  // Recupera solucao
  bestSolution = SCIPgetBestSol(scip);
  zz = SCIPgetPrimalbound(scip);//SCIPsolGetOrigObj(bestSolution);
 
  if (zz > z){
    vv=0;
     for (i = 0; i < I.n; i++)
     {
        valor = SCIPgetSolVal(scip, bestSolution, vars[i]);
        x[i] = valor;
        if(x[i]>EPSILON){
          vv+=I.item[i].value;
          PRINTF("\n%d (peso=%d valor=%d) da mochila valor=%lf", i+1, I.item[i].weight, I.item[i].value, vv);
        }
     }
  }
  else{
    // recupera a solucao anterior ao destroy
    for (i = 0; i < saiu; i++){
      x[sair[i].label] = 1;
    }
    zz = z;
  }

  // libera memoria
  free(sair);
  // Destroi problema
  SCIP_CALL( SCIPfree(&scip) ); 
  return zz;
}

/* heuristica a ser implementada */
double heuristica(instanceT I, int tipo, double *x, int tempo, float tam, int cLNS, double porcDestroy, double* z0, double *tempoLNS)
{
  double z = 0.0, z_lns;
  double aux[2]; // vetor que salva o valor e o total de itens da solução inicial
  clock_t start, end;
  
  aux[0] = 0.0;
  aux[1] = 0.0;
  
  // Solução inicial
  switch(tipo)
  {
    case 3:
       heuristica_Gulosa(I, x, aux, tempo);
      break;
    case 4:
       heuristica_Aleatoria(I, x, aux);
      break;
    case 5:
      heuristica_relax_fix(I, x, tam, tempo, aux);
      break;
  }
  
#ifdef DEBUG
  printf("\nValor da solução inicial: %.1lf\n", aux[0]);
  printf("Total de itens selecionados: %.0lf\n", aux[1]);
#endif
  
  z = aux[0];
  *z0 = z;
  
  start=clock();
  if (cLNS)
  {
     z_lns = heuristica_LNS(I, x, z, tempo, porcDestroy);
     if(z_lns>z){
       z = z_lns;
     }
  }
  end = clock();
  *tempoLNS=((double) end-start)/CLOCKS_PER_SEC;  
  return z;
}

int heuristica_Gulosa (instanceT I, double *x, double *aux, int tempo)
{
  instanceT II; // Instancia usada para escolher os itens
  int a, i, capacRes;
  int n;
  double valor;
  
  II.item = (itemType*)malloc(I.n * sizeof(itemType));
  
  for (i = 0; i < I.n; i++)
  {
    II.item[i].label = I.item[i].label;
    // Escolher o tipo de ordernação desejada
    II.item[i].value = -(I.item[i].value);
    II.item[i].weight = I.item[i].weight;
  }
  
  qsort(II.item, I.n, sizeof(itemType), comparador);
  
  // Constroi a solucao
  n = 0;
  valor = 0;
  capacRes = I.C;
  for (i = 0; i < I.n && capacRes > 0; i++)
  {
    a = II.item[i].label; // Selecionado o proximo item "a" da lista ordenada
    if(II.item[i].weight <= capacRes){
       x[a] = 1.0;
       capacRes -= II.item[i].weight;
       valor += -II.item[i].value;
       n++;
    }
    else{
       x[a]=0.0;
    }
  }

  aux[0]=valor;
  aux[1]=n;

  free(II.item);
  return 0;
}

void heuristica_Aleatoria (instanceT I, double *x, double *aux)
{
  int *L, nL;
  int i, a, r;
  int capacRes;
  
  L = (int *)malloc(I.n * sizeof(int));
  nL = I.n;
  
  for (i = 0; i < I.n; i++)
  {
    L[i] = I.item[i].label; // Cria a lista de itens
  }
  
  capacRes = I.C;
  while (nL > 0 && capacRes > 0)
  {
    r = RandomInteger(0, nL-1);
    a = L[r]; // Seleciona o item "a" de forma aleatória
    // Remove o item selecionado da lista
    L[r] = L[--nL];
    
    if (capacRes >= I.item[a].weight) // Tenta colocar o item na mochila
    {
       x[a] = 1.0;
       aux[0] += I.item[a].value;
       aux[1] += 1.0;
       capacRes -= I.item[a].weight;
    }
    else
    {
       x[a] = 0.0;
    }    
  }
  
  free(L);
}


// Função auxiliar de comparacao para o qsort
int comparador(const void *valor1, const void *valor2)
{
  if ((*(itemType *)valor1).value > (*(itemType *)valor2).value)
  {
    return 1;
  }
  else if ((*(itemType *)valor1).value == (*(itemType *)valor2).value)
  {
    return 0;
  }
  else
  {
    return -1;
  }
}

//Função que valida a solução encontrada pela Matheurística
int validarSolucao(instanceT I, double *x, double z, int tipo)
{
  double somaValor;
  int somaPeso;
  int *coberto;
  
  if (tipo != 1)
  {
    coberto = (int*)calloc(I.m, sizeof(int));
    somaPeso = 0;
    somaValor = 0.0;
    for (int i = 0; i < I.n; i++)
    {
      if (x[i] > EPSILON)
      {
        //        somaPeso += I.item[i].weight;
        // calcula o peso dos elementos cobertos pelo item i
        for(int j = 0; j < I.m; j++){
          if(I.R[i][j]>0 && !coberto[j]){
            somaPeso += I.weight[j];// soma peso dos elementos cobertos pela primeira vez
            coberto[j] = 1; // marca o elemento como coberto
          }
        }
        somaValor += I.item[i].value;
      }
    }
    free(coberto);
  
    if (somaPeso > I.C)
    {
      printf("Solução excedeu a capacidade por %d!\n", (somaPeso-I.C));
      return -1;
    }
    
    if (fabs(somaValor-z) > EPSILON)
    {
      if (somaValor > z + EPSILON)
      {
        printf("A soma dos valores dos itens excede o valor retornado pela heurística. \nSoma: %.1lf, z: %.1lf\n", somaValor, z);
      }
      
      if (somaValor < z - EPSILON)
      {
        printf("A soma dos valores dos itens está menor que o valor retornado pela heurística. \nSoma: %.1lf, z: %.1lf\n", somaValor, z);
      }
      return -2;
    }
  }
  return 1;
}
/* resolve o problema de PLI usando o SCIP
   parameter info = guarda informações de execução do B&B
*/
double getDualBound(instanceT I, double *tempo_gasto)
{
  SCIP* scip;
  parametersT param;
  double dualBound;
#ifdef DEBUG
  int status;
#endif
  clock_t start, end;
  
  // create scip and set scip configurations
  param.nodes_limit = 1;
  param.time_limit = 10;
  param.display_freq = -1;
  configScip(&scip, param);

  // carga do lp
  // load problem into scip
  if(!loadProblem(scip, &I, 1, NULL)){
    printf("\nProblem to load instance problem\n");
    return SCIPinfinity(scip);
  }

  // solve scip problem
  start=clock();
  SCIP_CALL( SCIPsolve(scip) );
  end = clock();

  *tempo_gasto = ((double) end-start)/CLOCKS_PER_SEC;

#ifdef DEBUG
  status = SCIPgetStatus(scip);
  printf("\nstatus=%d", status);
#endif
  dualBound = SCIPgetDualbound(scip);  

#ifdef DEBUG
  // Grava PL
  PRINTF("\n---LP gravado em mochila.lp");
  // print problem
  SCIP_CALL( SCIPwriteOrigProblem(scip, "mochila.lp", "lp", FALSE) );
#endif
  // Destroi problema
  //  SCIP_CALL( SCIPfree(&scip) ); 
  return dualBound;
}

/* programa principal */
int main(int argc, char **argv)
{
  double z, *x, porcDestroy, z0;
  clock_t antes, agora;
  int tipo, tempo, validacao;
  instanceT* I;
  int cLNS;
  float tam;
  char *instance_filename;
  parametersT param;
  int nodes, ativos;
  double dualBound, gap, tempo_lp, finalDualBound, tempoLNS;
  
  srand(time(NULL)); // Inicialização da semente da função rand()

  // checa linha de comando
  if (argc < 6)
  {
    printf("\nSintaxe: mochila <instancia.txt> <tipo> <tempo Solver (em segundos)> <com ou sem LNS> <porcentagem de Destroy no LNS> [<percentual_tam particao> in (0,1)]\n\t<tipo>: 1 = relaxacao linear, 2 = solucao inteira, 3 (heuristica - Gulosa), 4 - (heuristica - Aleatoria), 5 (heuristica - Relax&Fix)\n");
    exit(1);
  }

  // ler a entrada
  if (!loadInstance(argv[1], &I))
  {
    printf("\nProblema na carga da instância: %s", argv[1]);
    exit(1);
  }

  tipo = atoi(argv[2]);
  if (tipo < 1 || tipo > 5)
  {
    printf("Tipo invalido\nUse: tipo=1 (relaxacao linear), 2 (solucao inteira), 3 (heuristica - Gulosa), 4 - (heuristica - Aleatória), 5 (heuristica - Relax&Fix)\n");
    exit(1);
  }
  
  tempo = atoi(argv[3]);
  cLNS = atoi(argv[4]);
  porcDestroy = atof(argv[5]);

  // aloca memoria para a solucao
  x = (double *)calloc(I->n, sizeof(double));
  
  // inicializa dual bound
  dualBound = getDualBound(*I, &tempo_lp);
  gap = 0;
  nodes = 0;
  ativos = 0;

  tam = argc>6?atof(argv[6]):3;
  
  antes = clock();
  if (tipo < 3)
  {
    param.time_limit = tempo;
    param.display_freq = -1;
    param.nodes_limit=-1;
    z = otimiza_PLI(*I, tipo, x, param, &finalDualBound,&nodes, &ativos);
  }
  else
  {
    z = heuristica(*I, tipo, x, tempo, tam, cLNS, porcDestroy, &z0, &tempoLNS);
  }
  agora = clock();

  PRINTF("\nValor da solucao: %lf\tTempo gasto=%lf\n", z, ((double)agora - antes) / CLOCKS_PER_SEC);

  if(z>0){
     validacao = validarSolucao(*I, x, z, tipo);
  }else{
     validacao = -1;
  }
  instance_filename = strrchr(argv[1], '/');
  instance_filename++;

  if(tipo!=2){
     gap = 100*(dualBound - z)/dualBound;
  }
  else{
     gap = 100*(finalDualBound-z)/finalDualBound;
  }

  // mostra dados da callback
  printf("%s;%d;%d;%d;%lf;%f;%d;%d;%.0lf;%lf;%lf;%.3lf;%.3lf;%d;%d;%d;%lf;%lf;%lf\n", instance_filename, tipo,tempo, cLNS,porcDestroy, tam, I->n, I->m, tipo!=1?z:0, tipo==2?finalDualBound:dualBound, gap, dualBound, tempo_lp, nodes, ativos, validacao, ((double)agora - antes) / CLOCKS_PER_SEC, cLNS?z0:z, tempoLNS);

  #ifdef NDEBUG
  if(z>0 && tipo!=1){
     salva_solucao(argv[1], tipo, z, x, I, validacao, tempo, cLNS, porcDestroy, tam);
  }
  #endif

  // libera memoria alocada
  freeInstance(&I);
  free(x);
  return 0;
}

/* eof */
