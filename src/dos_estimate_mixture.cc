/* dos_estimate_mixture.cc

Copyright 2018 Grant M. Rotskoff

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */
#include <armadillo>
#include "particles.h"
#include "gaussian_mixture.h"

int main(int argc, char **argv)
{

  // usage
  if (argc != 5)
  {
    fprintf(stderr,"usage ./run_mixture dim gamma n_traj n_well\n");
    exit(-1);
  }
  int dim = atoi(argv[1]);
  double gamma = atof(argv[2]);
  int n_traj = atoi(argv[3]);
  int n_well = atoi(argv[4]);

  arma_rng::set_seed(983948524);

  gaussian_mixture g;

  g.initialize(dim, n_well);

  // load in the parameters for the mog
  char buffer[50];
  sprintf(buffer, "../data/sigmas_%03d_%03d.dat", dim, n_well);
  g.sigmas.load(buffer);
  sprintf(buffer, "../data/depths_%03d_%03d.dat", dim, n_well);
  g.depths.load(buffer);
  sprintf(buffer, "../data/mus_%03d_%03d.dat", dim, n_well);
  g.mus.load(buffer);
  g.mus = g.mus.t();
  // set the integration parameters
  double dt = 1e-5;
  double tol = 1e-5;
  double hmin = -g.depths.max();
  double hmax = 450.;
  double hres = 0.2;
  int max_iter = 0;
  fprintf(stderr,"Setting hmin=%4.3f\n",hmin);

  // reinjection parameters
  g.qmin = -1.;
  g.qmax = 1.;
  g.pmin = 0; // (NB: using mc reinjection)
  g.pmax = 0;
  g.mc_stepsize = 1e-2;
  g.beta = 5;
  bool use_mc = true;

  g.initialize_integration_variables(dt, gamma);
  g.initialize_estimator(hmax, hres, hmin, n_traj);
  fprintf(stderr,"Initialized estimator\n");
  fprintf(stderr, "Running Mixture of Gaussians with");
  fprintf(stderr, "\n\t n_well=%d",g.n_well);
  fprintf(stderr, "\n\t dim=%d\n",g.dim);

  // save everything
  for (int i=0; i<n_traj; i++) {
    g.run_estimation_trajectory(tol, hmax, max_iter, i, langevin, use_mc);
    fprintf(stderr, "Ran traj %d\n", i);
  }

  char dos_filename[80];
  sprintf(dos_filename, "mog_dim=%04d_gamma=%05.3f_ntraj=%04d.dat", dim, gamma, n_traj);
  FILE *dosf = fopen(dos_filename, "w");
  g.dump_dos_estimate(n_traj, dosf);
}
