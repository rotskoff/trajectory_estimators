/* gaussian_mixture.h
 
Copyright 2018 Grant M. Rotskoff

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */

#include "particles.h"

class gaussian_mixture : public particles {
public:

  int n_well;
  mat mus;
  vec sigmas;
  vec depths;

  double gaussian(vec, double, double);
  void force();

  void initialize(int, int);
  void compute_gradient();
  double compute_conf_energy();
  double compute_hamiltonian_energy();

};

const double pi = 4*atan(1.0);

double gaussian_mixture::gaussian(vec mu, double sigma, double depth)
{
  double s2 = sigma*sigma;
  return exp(-dot(qs-mu,qs-mu)/(2*s2) + depth);
}


// note that the physical "force" differs by a sign
void gaussian_mixture::force()
{
  double Ni = 0, Ntot = 0;
  for (int i=0; i<n_well; i++) {
    Ni = gaussian(mus.col(i), sigmas(i), depths(i));
    fs += (qs-mus.col(i))/(sigmas(i)*sigmas(i))*Ni;
    Ntot += Ni;
  }
  fs = fs/Ntot;
}

void gaussian_mixture::initialize(int dim0, int n_well0)
{
  dim = dim0;
  n_well = n_well0;
  mus = randu(dim, n_well)-0.5;
  sigmas = 0.025*randu(n_well);
  depths = 50*randu(n_well);
  qs = zeros(dim);
  ps = zeros(dim);
  fs = zeros(dim);
}

void gaussian_mixture::compute_gradient()
{
  fs*=0.;
  force();
}

// compute configurational
double gaussian_mixture::compute_conf_energy()
{
  double etot = 0;
  for (int i=0; i<n_well; i++) {
    etot += gaussian(mus.col(i), sigmas(i), depths(i));
  }
  return -log(etot);
}

double gaussian_mixture::compute_hamiltonian_energy()
{
  return dot(ps,ps) + compute_conf_energy();
}
