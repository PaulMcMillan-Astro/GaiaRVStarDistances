/*-----------------------------------------------------------------------------
*
* Prior.h
* Author: PJM
* email: paul@astro.lu.se
* Github: https://github.com/PaulMcMillan-Astro/
*
* Description: Header file for priors used in Bayesian distance determination.
*    Specifically here: For distances derived from parallaxes given in Gaia's
*    DR2
*-----------------------------------------------------------------------------*/


#ifndef _PRIOR_BASIC_
#define _PRIOR_BASIC_ 1

#include <string>
using namespace std;


/*-----------------------------------------------------------------------------
*
* Class: Basic_prior
* Purpose: Gives probability p(s)ds as a function of distance s (in pc),
* given l & b (in radians) initially.
*
* Currently, a fixed density model. Any changes have to be done here.
*   possibly flexible in future.
*-----------------------------------------------------------------------------*/

class Basic_prior {
  const double
    R0, z0,                            // solar position
    iRd_thin, izd_thin,                // properties of thin disc
    norm_thick, iRd_thick, izd_thick,  // properties of thick disc
    norm_halo, rcut_halo, power_halo;  // properties of halo
  const double qb, gammab,betab,r0b,rcutb; // properties of bulge
  double rho0b;                            // '' - not const for ease of setup
  void setup(double s);                    // Used by prior(s)
  //void setup(valarray<double>);

  double R, z, r;                      // Position in galaxy
  double l, b;                         // Galactic coords (radians)
  double cosl, sinl, cosb, sinb;       // Saved for speed
  bool offgrid;                        // If out-of-range
 public:
  void prepare_sightline(double L, double B); // Setup new sightline (l,b)
  double prior(double s);                     // Return prior given s
  double prior_thin();                        // Thin disc component
  double prior_thick();                       // Thick disc component
  double prior_halo();                        // Halo component
  double prior_bulge();                       // Bulge component
  Basic_prior();                              // Constructor

  ~Basic_prior() {};                          // Destructor
};



/*-----------------------------------------------------------------------------
*
* Class: GRVS_prior
* Purpose: Gives probability p(s,G_{RVS})ds dG_{RVS}, given l & b
*   initially. Note that the density part is done by Basic_prior
*
*   P(M_{G_{RVS}}) is a log-linear fit (plus Gaussian red clump) to the
* pdf found from PARSEC isochrones. Default is as described in McMillan (2018)
*
* Currently: two options for photometric prior. Possibility to add more.
*-----------------------------------------------------------------------------*/

class GRVS_prior {
  Basic_prior Non_mag_prior;            // Handles density

  /*---------------- Parameters of the photometric prior  --------------------*/
  double Mmin, Mmax,  div0_photprior, div1_photprior,  div2_photprior ,
    m0_photprior, m1_photprior, m2_photprior, m3_photprior ,
    c0_photprior, c1_photprior , c2_photprior, c3_photprior,
    RC_norm, RC_cen, RC_sig;

  public:
  double M_GRVS;                        // abs magnitude - to be shared around
  void prepare_sightline(double, double);         // Setup new sightline (l,b)
  double prior(double s_in, double GRVS_in);      // Return prior given s, G_RVS
  double prior_mag();                             // Abs magnitude prior
  //void setup_flat();         // deprecated
  //void setup_thin();         // deprecated
  //void setup_thin_DR2_old(); // decrepated
  void setup_flat_DR2();       // Prior assuming flat metallicity & age dist.
  void setup_thin_DR2();       // Prior assuming thin-disc-like distribution
  GRVS_prior();                // Constructor - default prior
  GRVS_prior(string type);     // Constructor - choice of prior
  ~GRVS_prior() {};            // Desctructor
};



/*-----------------------------------------------------------------------------
*
* Class: Dusty_prior
* Purpose: Gives probability p(s,G_{RVS})ds dG_{RVS}, given l & b
*   initially. Note that the density part is done by Basic_prior
*
*   TBD.
*-----------------------------------------------------------------------------*/

//
// class Dusty_prior {
// public:
//   GRVS_prior DustFree_prior;
//   double h_R, h_z,  A_V0, R_w, gamma_w,
//     R_fl,gamma_fl;
//   const int ngrid_AV;
//   double *A_Vgrid, *sgrid;
//   double rho_dust(double);
//   void get_A_V_table();
//   double find_AV();
//   double M_GRVS;
//   double Mmin, Mmax,  div1_photprior,  div2_photprior ,  m1_photprior,
//     m2_photprior, m3_photprior , c1_photprior , c2_photprior, c3_photprior,
//     RC_norm, RC_cen, RC_sig;
//
//   void prepare_sightline(double, double,double);
//   double prior(double, double,double);
//   double prior_mag();
//   void setup_flat();
//   void setup_thin();
//   Dusty_prior();
//   Dusty_prior(string type);
//   ~Dusty_prior() {};
// };



#endif
