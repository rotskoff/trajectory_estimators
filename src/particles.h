/* particles.h
 
Copyright 2018 Grant M. Rotskoff

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */

#ifndef PARTICLES_H
#define PARTICLES_H

#include <armadillo>
#include <math.h>
using namespace arma;

enum Method {langevin, gd};

class particles {
public:

  int dim; /**< Dimensionality of configuration space. */
  mat qs;  /**< Particle positions. */
  mat ps;  /**< Particle momenta. */
  mat fs;  /**< Particle forces. */

  /* reinjection */
  double qmin;
  double qmax;
  double pmin;
  double pmax;
  double beta;
  double mc_stepsize;
  virtual void reinject();
  void reinject_mc();


  /* integration */
  double dt;
  double gamma;
  double edt;
  double edtu;
  double edtr;
  double edtru;

  /* DOS estimator */
  double hmax;
  double hres;
  double hmin;
  double base_volume;
  double prior_volume;
  vec bins;
  vec tbs;
  vec tfs;
  vec dos;
  vec doserr;
  vec final_hs;
  void initialize_estimator(double hmax0, double hres0, double hmin0, int n_traj);
  void update_dos_estimator_forward(double h, double t, int init_bin);
  void update_dos_estimator_backward(double h, double t, int init_bin);

  particles();
  ~particles();

  virtual void compute_gradient();
  virtual double compute_conf_energy();
  virtual double compute_hamiltonian_energy();

  void initialize_integration_variables(double dt0, double gamma0);


  virtual void langevin_step_backward();
  virtual void langevin_step_forward();
  virtual void gd_step_backward();
  virtual void gd_step_forward();
  void run_traj_forward(double tol, int max_iter, int traj, Method);
  void run_traj_backward(double hmax, int max_iter, Method);
  void run_estimation_trajectory(double tol, double hmax, int max_iter, int traj, Method, bool mc);
  void compute_dos_estimate();
  void dump_dos_estimate(int n_traj, FILE *dosf);
};


#endif
