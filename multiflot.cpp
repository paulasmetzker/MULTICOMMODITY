// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This is an implementation of the Volume algorithm for multicommodity
// capacited fixed charge network design problems. See the references
// F. Barahona, R. Anbil, "The Volume algorithm: producing primal solutions
// with a subgradient method," IBM report RC 21103, 1998.
// F. Barahona, F. Chudak, "Solving large scale uncapacitated facility
// location problems," IBM report RC 21515, 1999.
//#include <valarray>
#include <cstdio>
#include <iostream>
#include <set>
#include <list>
#include <vector>
#include <queue>
#include <cmath>
#ifndef WIN32
#include <sys/times.h>
#endif
#include "multiflot.hpp"

using namespace std;

void FLOT_read_data(const char* fname, FLOT& data);

int main(int argc, char* argv[]) {

   // read in problem specific parameters and initialize data structures
   FLOT_parms flot_par("multiflot.par"); 
   FLOT flot_data;
   // read data
   FLOT_read_data(flot_par.fdata.c_str(), flot_data);

   // create the VOL_problem from the parameter file
   VOL_problem volp("multiflot.par");
   volp.psize = flot_data.commd * flot_data.arc + flot_data.arc;
   
   volp.dsize = flot_data.arc;
   
   bool ifdual = false;
   if (flot_par.dualfile.length() > 0) {
     // read dual solution
      ifdual = true;
      VOL_dvector& dinit = volp.dsol;
      dinit.allocate(volp.dsize);
      // read from file
      FILE * file = fopen(flot_par.dualfile.c_str(), "r");
      if (!file) {
	 printf("Failure to open file: %s\n ", flot_par.dualfile.c_str());
	 abort();
      }
      const int dsize = volp.dsize;
      int idummy;
      for (int i = 0; i < dsize; ++i) {
	 fscanf(file, "%i%lf", &idummy, &dinit[i]);
      }
      fclose(file);
   }

   volp.dual_lb.allocate(volp.dsize);
   volp.dual_lb = -1e31;
   volp.dual_ub.allocate(volp.dsize);
   volp.dual_ub = 0;													/* problema Ax <= b *, entao UB(dual) = 0 */


#ifndef WIN32
   // start time measurement
   double t0;
   struct tms timearr; clock_t tres;
   tres = times(&timearr); 
   t0 = timearr.tms_utime; 
#endif

    //invoke volume algorithm
   if (volp.solve(flot_data, ifdual) < 0) {
      printf("solve failed...\n");
   } 
   cout << "resolveu"<< endl;
  /* else {
      // recompute the violation
     // const int n = flot_data.node;
      const int a = flot_data.arc;
      const int k = flot_data.commd;

      VOL_dvector v(volp.dsize);
      const VOL_dvector& psol = volp.psol;
      v = 1;
      for (int l = 1; l <= k; ++l){
		for (int m = 0; m < a; ++m) {
		  //calcula violação
		}
      }

      double vc = 0.0;
      for (int l = 1; l <= k; ++l)
		for (int m = 0; m < a; ++m)
         // MS Visual C++ has difficulty compiling 
         //   vc += std::abs(v[i]);
         // so just use standard fabs.
         // CORRIGIR CALCULO ! vc += fabs(v[i]);
      vc /= a;
      printf(" Average violation of final solution: %f\n", vc);

     if (flot_par.dual_savefile.length() > 0) {
		// save dual solution
		 FILE* file = fopen(flot_par.dual_savefile.c_str(), "w");
		 const VOL_dvector& u = volp.dsol;
		 int n = u.size();
		 for (int t = 0; t < n; ++t) {
			fprintf(file, "%8i  %f\n", t+1, u[t]);							// melhorar para saber qual o arco dessa var dual
		 }
		 fclose(file);
      }

      // run a couple more heuristics
      // double heur_val;
     // for (i = 0; i < ufl_par.h_iter; ++i) {
	 //heur_val = DBL_MAX;
	 //ufl_data.heuristics(volp, psol, heur_val);
      //} 
      // save integer solution
      //if (flot_par.int_savefile.length() > 0) {
	 //FILE* file = fopen(ufl_par.int_savefile.c_str(), "w");
	 //const VOL_ivector& x = ufl_data.ix;
	 //const int n = ufl_data.nloc;
	 //const int m = ufl_data.ncust;
	 //int i,j,k=n;
	 //fprintf(file, "Open locations\n");
	 //for (i = 0; i < n; ++i) {
	   //if ( x[i]==1 )
	  //  fprintf(file, "%8i\n", i+1);
	// }
	// fprintf(file, "Assignment of customers\n");
	 //for (i = 0; i < n; ++i) {
	   //for (j = 0; j < m; ++j) {
//	     if ( x[k]==1 ) 
	  //     fprintf(file, "customer %i  location %i\n", j+1, i+1);
	//     ++k;
//	   }
//	 }
//	 fclose(file);
   //   }
   }

//   printf(" Best integer solution value: %f\n", ufl_data.icost);
//   printf(" Lower bound: %f\n", volp.value);
 /*  
#ifndef WIN32
   // end time measurement
   tres = times(&timearr);
   double t = (timearr.tms_utime-t0)/100.;
   printf(" Total Time: %f secs\n", t);
#endif
*/
   return 0;
}



//############################################################################

void FLOT_read_data(const char* fname, FLOT& data)
{

   FILE * file = fopen(fname, "r");
   if (!file) {
      printf("Failure to open ufl datafile: %s\n ", fname);
      abort();
   }


    VOL_dvector &fcost = data.fcost;
	VOL_dvector &u = data.u;
	VOL_dvector &c_cost = data.c_cost;
	VOL_dvector &demand = data.demand;
	
	VOL_ivector &orig = data.orig;
	VOL_ivector &dest = data.dest;
	VOL_ivector &i = data.i;
	VOL_ivector &j = data.j;
	
   int& node = data.node;
   int& arc = data.arc;
   int& commd = data.commd;
   int len;

   char s[500];
   fgets(s, 500, file);
   len = strlen(s) - 1;
   if (s[len] == '\n')
      s[len] = 0;
   // read number of node , number of arcs and number of commodities
   sscanf(s, "%d%d%d", &node, &arc, &commd);
    
   i.allocate(arc);
   j.allocate(arc);
   c_cost.allocate(arc*commd);
   u.allocate(arc);
   fcost.allocate(arc);
   
   orig.allocate(commd);
   dest.allocate(commd);  
   demand.allocate(commd);
   
   int ii;
   int jj;
   double cost_v;
   double cap;
   double cost_f;
   
   double dmd;
   
   int l,k;
   // read location costs
   //para cada arco
   for (l = 0; l < arc; ++l) { 
     fgets(s, 500, file);
     len = strlen(s) - 1;
     if (s[len] == '\n')
	s[len] = 0;
     sscanf(s,"%d%d%lf%lf%lf",&ii, &jj, &cost_v, &cap, &cost_f); // read to each arc the pair origin-destination, the trasnportation
																//cost, the capacity and the designer cost
     i[l] = ii;
     j[l] = jj;
     for(int k = 0 ; k < commd ; k++)
          c_cost[ l + k * arc] = cost_v;
     u[l] = cap;
     fcost[l]=cost_f;
   }
   
   for (l = 0; l < commd; ++l) { 
     fgets(s, 500, file);
     len = strlen(s) - 1;
     if (s[len] == '\n')
	s[len] = 0;
     k = sscanf(s,"%d%d%lf",&ii, &jj, &dmd); // read to each commodity the demand
     if(k != 3) break;
     if(ii == -1)break;
     orig[l] = ii;
     dest[l] = jj;
     demand[l] = dmd;
    }

   fclose(file);

 }

//###### USER HOOKS		
int
FLOT::all_negative(const VOL_ivector &N)
{
	int rt = 0;
	for(int it = 0 ; it < N.size() ; it++)
		if (N[it] == 1){
			rt = -1;
			return rt;
		}
	return rt;
}												
																		/** como terei o segundo coeficiente lagrangiano */
// compute reduced costs
int
FLOT::compute_rc(const VOL_dvector& dual_rc, VOL_dvector& rc)
{
	rc = 0;
	if(dual_rc.size() != 0){
		
		int it = arc * commd;
		for (int a = 0; a < arc; a++){
			for (int k = 0 ; k < commd ; k++){
			   rc[a + k * arc] = c_cost[a + k * arc] + dual_rc[a] ;		
			} 
			rc[it] = fcost[a] - (u[a] * dual_rc[a]) ;
			it++;
		}
	}
	else{
		rc = 0;
		int it = arc * commd;
		for (int a = 0; a < arc; a++){
			for (int k = 0 ; k < commd ; k++){
			   rc[a + k * arc] = c_cost[a + k * arc] ;	
			} 
			rc[it] = fcost[a] ;
			it++;
		}
	}
	return 0;
}

int
FLOT::find_index_min_cost(const VOL_dvector& cost, const VOL_ivector& N){

	double min = 1e31;
		int id = -1;//int id = NULL;
		for (int tmp = 0 ; tmp < node ; tmp ++)
			if(N[tmp] == 1 && cost[tmp] < min){
				min= cost[tmp];
				id = tmp;
			}
		return id;
}

// IN: dual vector dual
// OUT: primal solution to the Lagrangian subproblem (x) 
//      optimal value of Lagrangian subproblem (lcost)
//      v = difference between the rhs and lhs when substituting
//                  x into the relaxed constraints (v)
//      objective value of x substituted into the original problem (pcost)
//      xrc
// return value: -1 (volume should quit) 0 ow
  /** Solve the subproblem for the subgradient step.
       @param dual (IN) the dual variables
       @param rc (IN) the reduced cost with respect to the dual values
       @param lcost (OUT) the lagrangean cost with respect to the dual values
       @param x (OUT) the primal result of solving the subproblem
       @param v (OUT) b-Ax for the relaxed constraints
       @param pcost (OUT) the primal objective value of <code>x</code>
   */
int
FLOT::solve_subproblem(const VOL_dvector& dual, const VOL_dvector& rc,
		      double& lcost, VOL_dvector& x, VOL_dvector& v, double& pcost)
{
	
	x = 0;
	x = 0;
	v = 0;
	lcost = 0;
	pcost = 0;
	//Dijkstra 
	for (int k = 0 ; k < commd ; k++){
		int e = -1 ;
		VOL_ivector N(node);	
		VOL_dvector cost(node);
		VOL_ivector pred(node);
		cost = 1e29;
		pred = -1;
		//pred = NULL;
		N = 1;
		int index = orig[k] ;
		cost[index - 1] = 0;
		while (all_negative(N) != 0 ){
			e = find_index_min_cost(cost, N);
			if (e == -1){
				cout << " RETORNO DE MINIMO INVALIDO " << endl;
				break;}
			N[e] = -1;
			for(int a = 0 ; a < arc ; a++){
				int f = (j[a]-1);
				if ( ( i[a] == (e+1) ) && ( N[f] == 1 ) ){
						if (cost[f] > ( cost[e] + rc[ a + k * arc] ) ){
							pred[f] = e;
							cost[f] = cost[e] + rc[a + k * arc];
							} // end if
						if( f == (dest[k] - 1) ){
							cout << " Chegou no destino " << endl;
							goto continuation ;
						}
						
				} // end if
			} // end for
		} // end while

		continuation: 
		cout << "solve to x" << endl;								
		for (int it = ( dest[k] -1) ; it != ( orig[k] - 1 ) ; it = pred[it]){
			for(int a = 0 ; a < arc ; a++){
				if(i[a] ==( pred[it] + 1 ) && j[a] == ( it + 1) )	// pred[it] et it sont values [0, node-1] et i[] j[] de 1 a node
					x[a + k * arc] = demand[k];
			}
		}
	} // end for commodity
	// end Dijkstra   daiquistra
	// atribuição dos y (decisao se usa ou não o arco)
	int it = 0;
	cout << "solve to y" << endl;
	for(int ind = arc * commd ; ind < (arc * commd) + arc ; ind ++){
		double sum = ( fcost[it] - (dual[it] * u[it]) ) ; 
		if (sum <= 0)
			x[ind] = 1; 
		else
			x[ind] = 0;
		it ++;
	} // end construction solution in y
	
	//calculando subgradiente para cada arc (em relação ao y)
	double sum = 0;
	it = arc * commd ;
	for (int a = 0 ; a < arc ; a++){
		for(int k = 0 ; k < commd ; k++){
			pcost += ( c_cost[a + k * arc] ) * x[a + k * arc] ;
			lcost += ( c_cost[a + k * arc] + dual[arc] ) * x[a + k * arc] ;
			sum += x[a + k * arc];
			v[a + k * arc] = x[a + k * arc] - ( demand[k] * x[it] ) ;
		}
		pcost += fcost[a] * x[it];
		lcost += ( fcost[a] - (dual[a] * u[a]) ) * x[it];
		v[it] = sum + u[a];
		it++;
	}
   return 0;
}

// IN:  fractional primal solution (x),
//      best feasible soln value so far (icost)
// OUT: integral primal solution (ix) if better than best so far
//      and primal value (icost)
// returns -1 if Volume should stop, 0/1 if feasible solution wasn't/was
//         found.
// We use randomized rounding. We look at x[i] as the probability
// of opening facility i.
int
FLOT::heuristics(const VOL_problem& p,
		const VOL_dvector& x, double& new_icost)
{
   return 0;
}
