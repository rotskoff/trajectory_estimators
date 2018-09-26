/* dos_estimate_quartic.cc
 
Copyright 2018 Grant M. Rotskoff

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */

#include <armadillo>
#include "particles.h"
#include "quartic.h"


int main(int argc, char **argv)
{

  arma_rng::set_seed(45632426);
  quartic g;
  int dim = atoi(argv[1]);
  double gamma = atof(argv[2]);
  g.initialize(dim, gamma);
  double hmin = 0.;
  double hmax = 2. * dim;
  double hres = 0.5;

  g.qmin = -pow(4*hmax,0.25);
  g.qmax = pow(4*hmax,0.25);

  // set the integration parameters
  double dt = 1e-3;
  double tol = 1e-6;
  int max_iter = 0;
  bool use_mc = false;

  int n_traj = atoi(argv[3]);
  g.initialize_integration_variables(dt, gamma);
  g.initialize_estimator(hmax, hres, hmin, n_traj);
  // run the estimation trajectories
  for (int i=0; i<n_traj; i++) {
    g.run_estimation_trajectory(tol, hmax, max_iter, i, langevin, use_mc);
  }
  char dos_filename[80];
  sprintf(dos_filename, "quartic_dim=%03d_gamma=%05.3f_ntraj=%04d.dat", dim, gamma, n_traj);
  FILE *dosf = fopen(dos_filename, "w");
  g.dump_dos_estimate(n_traj, dosf);
}
