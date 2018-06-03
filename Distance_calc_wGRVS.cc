/*-----------------------------------------------------------------------------
*
* Distance_calc_wGRVS.cc
* Author: PJM
* email: paul@astro.lu.se
* Github: https://github.com/PaulMcMillan-Astro/
*
* Description: Code which calculates distances based on a Bayesian scheme, a la
* McMillan 2018
*
*-----------------------------------------------------------------------------*/



#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>


#include "Pi.h"
#include "Prior.h"
#include "Gauss-Legendre-header.h"

using namespace std;

// Value of a Gaussian, given value, mean, sigma
template <class T>
   T gaussian(T x, T x0, T sig) {
   T isig = 1./sig, tmp = (x-x0)*isig, tmp2 = -0.5*tmp*tmp;
   const T isrttpi = 1./ 2.506628274631000502415765284811;
   return isrttpi*isig*exp(tmp2);
 }

// Value of a Gaussian, given value, mean, sigma, unnormalised
// Only use when making comparisons!
 template <class T>
   T gaussian_unnorm(T x, T x0, T sig) {
   T isig = 1./sig, tmp = (x-x0)*isig, tmp2 = -0.5*tmp*tmp;
   return exp(tmp2);
 }

/*--               Give G_RVS, where input is G and GRP                     --*/
/*--    Formula given by Gaia Collaboration, Brown et al 2018 (eq 2 & 3)    --*/

 double get_GRVS(double G, double R) {
   double GRVS;
   double C = G-R;
   if(C < 1.4) {
     // extrapolation beyond blue end as constant (rarely important - Teff cut)
     if(C<0.1) C=0.1;
     GRVS = R + 0.042319-0.65124*C + 1.0215*C*C - 1.3947*C*C*C + 0.53768*C*C*C*C;
   } else {
     // extrapolation beyond red end - also constant
     if(C>1.7) C = 1.7;
     GRVS = R + 132.32 - 377.28*C + 402.32*C*C - 190.97*C*C*C + 34.026*C*C*C*C;
   }
   return GRVS;
 }


const double deg2rad  = Pi/180.;

/* ----- Corrections, taken from Lindegren et al 2018 (table 4)  -------------*/
const double gaia_parallax_zeropoint = -0.029;
const double gaia_correlated_parallax_error = 0.043;



/*

Main programme

*/
int main(int argc, char *argv[] ) {
// Prior chosen - can alter by hand
  GRVS_prior Prior("thin_DR2");

  ifstream from;
  ofstream to;
  if(argc>2) {
    from.open(argv[1]);
    to.open(argv[2]);
  } else {
    cerr << "Input and output file needed\n"
    << "Input lines: source_id l b parallax_obs parallax_obs_error Gmag RPmag\n"
    << "Space separated, and make sure that no values are missing!\n";
    return 1;
  }
  // input
  double l_in,b_in,parallax_obs_in, parallax_obs_error_in;
  double l,b,parallax_obs, parallax_obs_error, Gmag, RPmag,GRVS;
  // derived (integrated over) quantities
  double p,s,ls, s_pc;
  // integration parameters
  int nintegrate = 512;                  // Only possibility allowed
  double max_distance_considered = 100.; // in kpc
  double parmin, parmax,smin, smax,lsmin,lsmax, ls_centre, ls_hwidth;
  // Outputs from integration
  double prob_int, norm=0., s_av=0., s2_av=0.;
  double distance_out, distance_error_out;
  string source_id;

  to << "source_id,distance,distance_error\n";

  while(true) {
  // Read in data - l, b, parallax, parallax_error, G, G_RP
    from  >> source_id >> l_in >> b_in >> parallax_obs_in;
    if(from.eof()) break;
    from >> parallax_obs_error_in >> Gmag >> RPmag;

    // Correct for zero-point & systematic errors
    parallax_obs = parallax_obs_in - gaia_parallax_zeropoint;
    parallax_obs_error = sqrt(pow(gaia_correlated_parallax_error,2)
                              + pow(parallax_obs_error_in,2));

    // Convert to Gaia G_RVS magnitude
    GRVS = get_GRVS(Gmag,RPmag);
    norm=0.; s_av=0.; s2_av=0.;
    // convert l,b to radians
    l = l_in*deg2rad;
    b = b_in*deg2rad;

    Prior.prepare_sightline(l, b);

    // Set limits in sensible range given parallax
    parmin = parallax_obs-5*parallax_obs_error;
    parmax = parallax_obs+5*parallax_obs_error;
    smin = 1./parmax; // in kpc
    smax = 1./parmin; // in kpc

    // Make sure we're considering a sensible range
    if(smax<0.|| smax> max_distance_considered) smax = max_distance_considered;
    if(smin<0. || smin> max_distance_considered) smin = 0.1*max_distance_considered;

    // Integral in log distance - deals with different scales well.
    lsmin = log(smin); lsmax = log(smax);
    ls_centre = 0.5*(lsmin+lsmax);
    ls_hwidth = 0.5*(lsmax-lsmin);
    /* --------   Perform Gauss-Legendre integral          ---------*/
    for (int i=0;i!=nintegrate/2;i++){
      // Gauss-Legendre integration
      ls = x512[i]*ls_hwidth + ls_centre;
      s = exp(ls); // could speed things up with multiplication
      p = 1./s;
      // Prior is in pc (historical reasons)
      s_pc = 1000.*s;
      prob_int = w512[i]*Prior.prior(s_pc,GRVS)*s_pc*
       gaussian_unnorm<double>(p,parallax_obs,parallax_obs_error);
      norm  += prob_int;
      s_av  += s_pc*prob_int;
      s2_av += s_pc*s_pc*prob_int;
      //to << x512[i] << ' ' << s << ' ' << prob_int<< '\n'<< std::flush;

      // Integrate either side of centre
      ls = -x512[i]*ls_hwidth + ls_centre;
      s = exp(ls); // could speed things up with multiplication
      p = 1./s;
      s_pc = 1000.*s;
      prob_int = w512[i]*Prior.prior(s_pc,GRVS)*s_pc*
       gaussian_unnorm<double>(p,parallax_obs,parallax_obs_error);
      norm  += prob_int;
      s_av  += s_pc*prob_int;
      s2_av += s_pc*s_pc*prob_int;

    }
    distance_out = s_av/norm;
    distance_error_out = sqrt(s2_av/norm - distance_out*distance_out);
    to << source_id << ',' << distance_out << ',' << distance_error_out << '\n';
    /*to << source_id << ' ' <<  l_in << ' ' <<  b_in << ' '
      <<  parallax_obs_in << ' ' <<  parallax_obs_error_in << ' ' <<  Gmag << ' ' <<  RPmag
      << ' ' << distance_out << ' ' << distance_error_out << '\n'; */
  }
}
