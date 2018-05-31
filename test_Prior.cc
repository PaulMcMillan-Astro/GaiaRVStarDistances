
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
  }
  ofstream to(argv[1]), to2(argv[2]);

  Basic_prior Prior;
  GRVS_prior Prior_mag("thin_DR2");

  // test Basic Prior

  double l=0.,b=0.,s, G_RVS=11.;
  Prior.prepare_sightline(l*deg2rad,b*deg2rad);
  for(s=10;s<10000;s+=100) {
    //to << s << ' ' << Prior.prior(s) << '\n';
    Prior.prior(s);
    to << s << ' ' << Prior.prior_thin() << ' ' << Prior.prior_bulge() << '\n';
  }

  Prior_mag.prepare_sightline(l*deg2rad,b*deg2rad);
  for(s=10;s<100000;s+=10) {
    double tmp = Prior_mag.prior(s,G_RVS);
    to2 << s << ' ' << Prior_mag.M_GRVS << ' ' << Prior_mag.prior_mag() << '\n';
  }

  //completely dumb mass integration test
  double d=20., mb = 0., md=0.;
  for(double x=d;x<15000.;x+=d)
    for(double y=d;y<15000.;y+=d)
      for(double z=d;z<2000.;z+=d) {
        double yt = y-8330.;
        l = atan2(x,yt);
        b = atan2(z,sqrt(x*x+yt*yt));
        Prior.prepare_sightline(l,b);
        s = sqrt(x*x+yt*yt+z*z);
        Prior.prior(s);
        mb += Prior.prior_bulge()*d*d*d*8;
        md += (Prior.prior_thin()+Prior.prior_thick())*d*d*d*8;
      }
  cerr << mb << ' ' << md << ' ' << mb/md << '\n';

  return 0;

}
