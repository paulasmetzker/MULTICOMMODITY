// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef __MULTIFLOT_HPP__
#define __MULTIFLOT_HPP__
/////////// FAZER INDEXAÇÃO A * COMMODITY E SALVAR ASSIM
#include <cstdio>
#include <cfloat>
#include <string>
#include <cstring>


#include "VolVolume.hpp"

// parameters controlled by the user
class FLOT_parms {
public:
  std::string fdata; // file with the data
  std::string dualfile; // file with an initial dual solution
  std::string dual_savefile; // file to save final dual solution
  std::string int_savefile; // file to save primal integer solution
  
  //int h_iter; // number of times that the primal heuristic will be run
  // after termination of the volume algorithm*/
   
  FLOT_parms(const char* filename);
  ~FLOT_parms() {}
};



class FLOT : public VOL_user_hooks { 
public:
  // for all hooks: return value of -1 means that volume should quit
  // compute reduced costs
   int compute_rc(const VOL_dvector& u, VOL_dvector& rc);
  // solve relaxed problem
   int solve_subproblem(const VOL_dvector& cost, const VOL_dvector& rc,
			double& lcost, VOL_dvector& x, VOL_dvector&v,
			double& pcost);
  // primal heuristic
  // return DBL_MAX in heur_val if feas sol wasn't/was found 
   int heuristics(const VOL_problem& p, 
		  const VOL_dvector& x, double& heur_val);
		  
   int all_negative(const VOL_ivector &N);
   int find_index_min_cost(const VOL_dvector& cost, const VOL_ivector& N);

  // original data for multicommodity capacited fixed
public:
	// parameters
  int node, arc, commd; // number of installations , number of arcs,  number of commodity
  VOL_dvector fcost; // designer cost
  VOL_dvector u; // capacity in each arc
  /** o custo de transporte foi assumido como o mesmo para qualquer commodity */
  VOL_dvector c_cost;// transportation cost
  VOL_dvector demand; // demande of each commodity
 
/** for each commodity K we must have a pair of origin and destination. Their indexs are save in 2 vectors  */
  VOL_ivector orig; // origin of each commodity K
  VOL_ivector dest; // destination of each commodity K
  
  VOL_ivector i; // origin of each arc
  VOL_ivector j; // destination of each arc
  
  /** solution */														/** indexação de x = i * commd + j  */
  VOL_ivector ix;   // best integer feasible solution so far
  double      icost;  // value of best integer feasible solution 
public:
  FLOT() : icost(DBL_MAX) {}
  ~FLOT() {}  				// ANTES ERA VIRTUAL 
};

//#############################################################################
//########  Member functions          #########################################
//#############################################################################

//****** FLOT_parms
// reading parameters specific to facility location
FLOT_parms::FLOT_parms(const char *filename) :
   fdata("") 
 //  h_iter(10)
{
   char s[500];
   FILE * file = fopen(filename, "r"); 
   if (!file) { 
      printf("Failure to open ufl datafile: %s\n ", filename);
      abort();
   }
   
   while (fgets(s, 500, file)) {
      const int len = strlen(s) - 1;
      if (s[len] == '\n')
	 s[len] = 0;
      std::string ss;
      ss = s;
      
      if (ss.find("fdata") == 0) {
	 int j = ss.find("="); 
	 int j1 = ss.length() - j + 1; 
	 fdata = ss.substr(j+1, j1); 
	 
      } else if (ss.find("dualfile") == 0) {
	 int j = ss.find("=");
	 int j1 = ss.length() - j + 1;
	 dualfile = ss.substr(j+1, j1); 
	 
      } else if (ss.find("dual_savefile") == 0) {
	 int j = ss.find("=");
	 int j1 = ss.length() - j + 1;
	 dual_savefile = ss.substr(j+1, j1); 

      } else if (ss.find("int_savefile") == 0) {
	 int j = ss.find("=");
	 int j1 = ss.length() - j + 1;
	 int_savefile = ss.substr(j+1, j1); 
	 
      } /*else if (ss.find("h_iter") == 0) {
	 int i = ss.find("=");  
	 h_iter = atoi(&s[i+1]);
      }	*/
   }
   fclose(file);
}

#endif
