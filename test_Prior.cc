/*-----------------------------------------------------------------------------
*
* test_Prior.cc
* Author: PJM
* email: paul@astro.lu.se
* Github: https://github.com/PaulMcMillan-Astro/
*
* Description: Simple tests for the priors used to find distance estimates
*   N.B. these just give outputs that can be plotted and considered, i.e. they
* are not tests that the compiled code can pass/fail. Just checks that the user
* can consider
*
*-----------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
#include <sstream>

#include "Pi.h"
#include "Prior.h"

using namespace std;



const double deg2rad  = Pi/180.;

int main(int argc, char *argv[] ) {

  if(argc<3){
    cerr << "I need two output filenames\n";
    exit(0);
  }
  ofstream to(argv[1]), to2(argv[2]);

  Basic_prior Prior;
  GRVS_prior Prior_mag("thin_DR2");

  // test Basic Prior

  double l=10.,b=30.,s, G_RVS=11.;
  Prior.prepare_sightline(l*deg2rad,b*deg2rad);
  for(s=10;s<10000;s+=100) {
    to << s << ' ' << Prior.prior(s) << '\n';
    //Prior.prior(s);
    //to << s << ' ' << Prior.prior_thin() << ' ' << Prior.prior_bulge() << '\n';
  }

    // test GRVS_Prior
  Prior_mag.prepare_sightline(l*deg2rad,b*deg2rad);
  for(s=10;s<100000;s+=10) {
    double tmp = Prior_mag.prior(s,G_RVS);
    // Just checking the magnitude part, not anything else
    to2 << s << ' ' << Prior_mag.M_GRVS << ' ' << Prior_mag.prior_mag() << '\n';
  }

  // Just a check that the bulge component (which is new) has been implimented
  // correctly
  double d=20., mb = 0., md=0.;
  // Integrate density dx, dy, dz (over 1/8th of the distribution)
  for(double x=d;x<15000.;x+=d)
    for(double y=d;y<15000.;y+=d)
      for(double z=d;z<2000.;z+=d) {
        // 8330 pc copied from the Basic_prior
        double yt = y-8330.;
        l = atan2(x,yt);
        b = atan2(z,sqrt(x*x+yt*yt));
        // Dumb - but the only way to put it into the code
        Prior.prepare_sightline(l,b);
        s = sqrt(x*x+yt*yt+z*z);
        Prior.prior(s);
        // Bulge mass
        mb += Prior.prior_bulge()*d*d*d*8;
        // Disc mass
        md += (Prior.prior_thin()+Prior.prior_thick())*d*d*d*8;
      }
  // Output bulge 'mass', disc 'mass' and ratio
  // Note that the actual masses are meaningless - only close to
  // reality because density near the sun is not too far from 1 M_sun/pc^3
  cerr << mb << ' ' << md << ' ' << mb/md << '\n';

  return 0;

}
