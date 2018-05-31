
#ifndef _PRIOR_BASIC_
#define _PRIOR_BASIC_ 1

#include <string>
using namespace std;

class Basic_prior {
  const double
    R0, z0, iRd_thin, izd_thin,  norm_thick,
    iRd_thick, izd_thick,
    norm_halo, rcut_halo, power_halo;
  const double qb, gammab,betab,r0b,rcutb;
  double rho0b; // Just makes life easier to have it here
  void setup(double);
  //void setup(valarray<double>);

  double R, z, r;
  double l, b;
  double cosl, sinl, cosb, sinb;
  bool offgrid;
 public:
  void prepare_sightline(double, double);
  double prior(double);
  double prior_thin();
  double prior_thick();
  double prior_halo();
  double prior_bulge();
  Basic_prior();

  ~Basic_prior() {};
};


class GRVS_prior {
public:
  Basic_prior Non_mag_prior;
  double M_GRVS;
  double Mmin, Mmax,  div0_photprior, div1_photprior,  div2_photprior ,
    m0_photprior, m1_photprior, m2_photprior, m3_photprior ,
    c0_photprior, c1_photprior , c2_photprior, c3_photprior,
    RC_norm, RC_cen, RC_sig;

  void prepare_sightline(double, double);
  double prior(double, double);
  double prior_mag();
  void setup_flat();
  void setup_thin();
  void setup_thin_DR2_old();
  void setup_flat_DR2();
  void setup_thin_DR2();
  GRVS_prior();
  GRVS_prior(string type);
  ~GRVS_prior() {};
};

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
