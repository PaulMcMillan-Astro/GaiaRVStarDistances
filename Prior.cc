
#include <cmath>
#include <iostream>
#include "Prior.h"

void Basic_prior::prepare_sightline(double L,double B) {
  l=L; b=B;
  cosl=cos(l); sinl=sin(l);
  cosb=cos(b); sinb=sin(b);
}

void Basic_prior::setup(double s) {
  // the below are slighty wrong because the Galactic coordinate
  // system isn't aligned with the plane of the disc.
  R =  sqrt(R0*R0 + s*s*cosb*cosb - 2.*R0*s*cosb*cosl);
  z = s*sinb+z0;
  r = sqrt(R*R+z*z);
  offgrid = ( s<0.)? true : false;
}


double Basic_prior::prior(double s) {
  setup(s);
  if(offgrid) return 0.;
  return (prior_thin()+prior_thick()+
	  prior_halo()+prior_bulge())* //position
    s*s;                                      //r^2 (cone)
}



Basic_prior::Basic_prior() :
  R0(8330.), z0(0.), iRd_thin(1./2600.), izd_thin(1./300.),
 norm_thick(0.0596),
  iRd_thick(1./3600.), izd_thick(1./900.),
  norm_halo(3.774e9),
  rcut_halo(3000.),
  power_halo(-3.39),
  qb(0.5), gammab(0.),betab(1.8),
  r0b(0.075*1000.), rcutb( 2.1*1000.)
{
  double Md = 8*3.141592*((1./izd_thin)*powf(1./iRd_thin,2)*1. +
                          (1./izd_thin)*powf(1./iRd_thin,2)*norm_thick);
  double MBonMS = 9.23/54.3;
  double Mb = MBonMS*Md/(1.-MBonMS);
  rho0b = 1e-9*97.4*(0.51/0.5)/8.96*Mb;
  //std::cerr << Md << ' ' << Mb << '\n';
}







double Basic_prior::prior_thin() {
  //  if(age>10.) return 0.; // age out of range
  const double isrttpi = 1./ 2.506628274631000502415765284811;
  double exponent = - R*iRd_thin-fabs(z)*izd_thin ;
  double out = isrttpi*exp(exponent);
  return out;
}

double Basic_prior::prior_thick() {
  //  if(age<8. || age>12.) return 0.; // age out of range
  const double isrttpi = 1./ 2.506628274631000502415765284811;
  double exponent = - R*iRd_thick-fabs(z)*izd_thick;       // position
  double out = norm_thick*isrttpi*exp(exponent);
  return out;
}

double Basic_prior::prior_halo() {
  if(r<rcut_halo) return norm_halo*pow(rcut_halo,power_halo);
  // power law diverges at low r - Carollo model applied where it shouldn't
  double out = norm_halo*pow(r,power_halo);
  return out;
}

double Basic_prior::prior_bulge() {
  double rho = rho0b, m = sqrt(R*R+(z/qb*z/qb)), m0 = m/r0b;
  if(gammab != 0.) rho *= powf(m0,-gammab);
  m0+=1.;
  rho *= powf(m0,-betab);
  //cerr << R << ' ' << rho << '\n';
  if(rcutb) rho *= exp(-powf(m/rcutb,2));
  //cerr << R << ' ' << rho << '\n';
  return rho;
}








void GRVS_prior::setup_flat_DR2() {

  Mmin = -12;
  Mmax = 7;

  div0_photprior  =  -2.34936054702 ;
  div1_photprior  =  1.04283066069 ;
  div2_photprior  =  3.6743145648 ;
  m0_photprior  =  0.636749185704 ;
  m1_photprior  =  0.1865447312 ;
  m2_photprior  =  0.599726656521 ;
  m3_photprior  =  0.144988560077 ;
  RC_cen  =  -0.435918445572 ;
  RC_sig  =  0.23199035803 ;
  RC_norm  =  8185.80862358 ;
  c1_photprior   =  3.40503108665 ;

  c2_photprior = c1_photprior+(m1_photprior- m2_photprior)*div1_photprior;
  c3_photprior = c2_photprior+(m2_photprior - m3_photprior)*div2_photprior;
  c0_photprior = c1_photprior+(m1_photprior- m0_photprior)*div0_photprior;
}
void GRVS_prior::setup_thin_DR2() {

  Mmin = -12;
  Mmax = 7;



div0_photprior  =  -2.02022468025 ;
div1_photprior  =  1.27873640026 ;
div2_photprior  =  3.5343533639 ;
m0_photprior  =  0.678238934106 ;
m1_photprior  =  0.213637402684 ;
m2_photprior  =  0.670052484961 ;
m3_photprior  =  0.141224575915 ;
RC_cen  =  -0.349737628723 ;
RC_sig  =  0.150875363233 ;
RC_norm  =  9919.71867715 ;
c1_photprior   =  3.3209869255 ;
  c2_photprior = c1_photprior+(m1_photprior- m2_photprior)*div1_photprior;
  c3_photprior = c2_photprior+(m2_photprior - m3_photprior)*div2_photprior;
  c0_photprior = c1_photprior+(m1_photprior- m0_photprior)*div0_photprior;

}




GRVS_prior::GRVS_prior() {
// Anything I need to do to set up the GRVS prior
  setup_thin_DR2();
}

GRVS_prior::GRVS_prior(string type) {
// Anything I need to do to set up the GRVS prior
  //if(type=="thin") setup_thin();
  //else if(type=="flat") setup_flat();
  //else if(type=="thin_DR2_old") setup_thin_DR2_old();
  //else
  if(type=="thin_DR2")      setup_thin_DR2();
  else if(type=="flat_DR2") setup_flat_DR2();
}
void GRVS_prior::prepare_sightline(double l, double b) {
    Non_mag_prior.prepare_sightline(l,b);
}

double GRVS_prior::prior(double s_in, double GRVS_in) {
  M_GRVS = GRVS_in - 5*log10(s_in/10.);
  return prior_mag() * Non_mag_prior.prior(s_in);
}

double GRVS_prior::prior_mag() {
  double log_pm;
  if(M_GRVS<Mmin) return 0.;
  else if(M_GRVS<div0_photprior) {
    log_pm = c0_photprior+m0_photprior*M_GRVS;
  }
  else if(M_GRVS<div1_photprior) {
    log_pm = c1_photprior+m1_photprior*M_GRVS;
  } else if(M_GRVS<div2_photprior) {
    log_pm = c2_photprior+m2_photprior*M_GRVS;
  } else if(M_GRVS<Mmax) {
    log_pm = c3_photprior+m3_photprior*M_GRVS;
  } else if(M_GRVS>=Mmax) {
    log_pm = c3_photprior+m3_photprior*Mmax;
  }
  double pm = pow(10.,log_pm);
  pm += RC_norm * exp(-pow(M_GRVS-RC_cen,2)/(2.*RC_sig*RC_sig));

  return pm;

}
