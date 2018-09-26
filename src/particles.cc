/* particles.cc
 
Copyright 2018 Grant M. Rotskoff

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */
#include "particles.h"
#include <assert.h>

particles::particles(){}

particles::~particles(){}

void particles::compute_gradient(){}

double particles::compute_conf_energy(){ return 0.; }

double particles::compute_hamiltonian_energy(){ return 0.; }

void particles::initialize_estimator(double hmax0, double hres0, double hmin0, int n_traj)
{
  hmax = hmax0;
  hres = hres0;
  hmin = hmin0;
  base_volume = 0;
  prior_volume = pow(qmax-qmin,dim)*pow(pmax-pmin,dim);
  bins = regspace<vec>(hmin+0.5*hres, hres, hmax+0.5*hres);
  tbs = bins*0.;
  tfs = bins*0.;
  dos = bins*0.;
  doserr = bins*0.;
  final_hs = zeros<vec>(n_traj);
}

void particles::update_dos_estimator_forward(double h, double t, int init_bin)
{
  int bin_index = ceil((h-bins(0))/hres);
  // record the last time at which energy is above bins(bin_index)
  if (bin_index<(int)bins.size() && bin_index>=0) {
    if (tfs(bin_index)==0 && bin_index!=init_bin) {
      tfs(bin_index) = t-dt;
    }
  }
}

void particles::update_dos_estimator_backward(double h, double t, int init_bin)
{
  int bin_index = ceil((h-bins(0))/hres);
  // record the last time at which energy is below bins(bin_index)
  if (bin_index<(int)bins.size() && bin_index>=0)
  {
    if (tbs(bin_index)==0 && bin_index!=init_bin) {
      tbs(bin_index) = t+dt;
    }
  }
}

void particles::reinject()
{
  qs = qmin + randu(size(qs))*(qmax-qmin);
  ps = pmin + randu(size(ps))*(pmax-pmin);
  double h = compute_hamiltonian_energy();
  double n_attempt = 0.;
  while (h>hmax)
  {
    qs = qmin + randu(size(qs))*(qmax-qmin);
    ps = pmin + randu(size(ps))*(pmax-pmin);
    h = compute_hamiltonian_energy();
    n_attempt += 1.0;
  }
  base_volume += 1.0/n_attempt;
}

void particles::reinject_mc()
{
  qs = qmin + randu(size(qs))*(qmax-qmin);
  ps = pmin + randu(size(ps))*(pmax-pmin);

  double h0 = compute_hamiltonian_energy();
  while (h0>5*hmax)
  {
    qs = qmin + randu(size(qs))*(qmax-qmin);
    ps = pmin + randu(size(ps))*(pmax-pmin);
    h0 = compute_hamiltonian_energy();
  }
  while (h0>1.01*hmax)
  {
    gd_step_forward();
    h0 = compute_hamiltonian_energy();
  }
  vec dqs = 0.*qs;
  while (h0>hmax)
  {
    dqs = dt * (randu(size(qs))-0.5);
    qs = qs + dqs;
    double h1 = compute_hamiltonian_energy();
    if (randu() < exp(-beta * (h1-h0)))
    {
      h0 = h1;
    }
    else {
      qs = qs - dqs;
    }
  }
}

void particles::initialize_integration_variables(double dt0, double gamma0) {
  dt = dt0;
  gamma = gamma0;
  edt = exp(-gamma*dt);
  edtu = (1-exp(-gamma*dt))/gamma;
  edtr = exp(gamma*dt);
  edtru = (exp(gamma*dt)-1)/gamma;
}

void particles::langevin_step_backward()
{
  qs = qs - dt*ps;
  compute_gradient();
  ps = edtr*ps + edtru*fs;
  qs = qs - dt*ps;
}

void particles::langevin_step_forward()
{
  qs = qs + dt*ps;
  compute_gradient();
  ps = edt*ps - edtu*fs;
  qs = qs + dt*ps;
}

void particles::gd_step_backward()
{
  compute_gradient();
  qs = qs + dt*fs;
}

void particles::gd_step_forward()
{
  compute_gradient();
  qs = qs - dt*fs;
}

void particles::run_traj_forward(double tol, int max_iter, int traj, Method method)
{
  int iter = 0;
  double t = 0;
  double h = compute_hamiltonian_energy();
  int init_bin = ceil((h-bins(0))/hres);
  compute_gradient();
  while (h>(hmin+hres) && norm(fs,2)>tol)
  {
    if (method==langevin) {
      langevin_step_forward();
    } else if (method==gd) {
      gd_step_forward();
    }
    t += dt;
    iter += 1;
    h = compute_hamiltonian_energy();
    update_dos_estimator_forward(h, t, init_bin);
  }
  final_hs(traj) = h;
}

void particles::run_traj_backward(double hmax, int max_iter, Method method)
{
  int iter=0;
  double t=0;
  double h=compute_hamiltonian_energy();
  int init_bin = floor((h-bins(0))/hres);
  compute_gradient();
  while (h<hmax) // && iter < max_iter)
  {
    if (method==langevin) {
      langevin_step_backward();
    } else if (method==gd) {
      gd_step_backward();
    }
    t-=dt;
    iter+=1;
    h=compute_hamiltonian_energy();
    update_dos_estimator_backward(h,t,init_bin);
  }
}

void particles::run_estimation_trajectory(double tol, double hmax,
                                          int max_iter, int traj, Method method,
                                          bool mc)
{
  tfs = tfs*0.;
  tbs = tbs*0.;
  if (mc) {
    reinject_mc();
  }
  else {
    reinject();
  }
  vec qs_init = qs;
  vec ps_init = ps;
  compute_gradient();
  run_traj_backward(hmax, max_iter, method);
  qs = qs_init;
  ps = ps_init;
  run_traj_forward(tol, max_iter, traj, method);
  compute_dos_estimate();
}

void particles::compute_dos_estimate()
{
  double dgam = dim*gamma;
  int nbins = bins.size();
  vec tdiffs = (tfs+tbs-tbs(nbins-1)) % (tfs>0 || tbs<0);
  dos += exp(-dgam*tdiffs) % (tfs>0 || tbs<0);
  doserr += dos % dos;
  //tdiffs.t().print();
}

void particles::dump_dos_estimate(int n_traj, FILE *dosf)
{
  dos = dos / n_traj;
  doserr = doserr / n_traj;
  doserr = sqrt((doserr - dos % dos)/n_traj);
  int nbins = bins.size();
  for (int i=1; i<nbins; i++){
    fprintf(dosf, "%e %e %e\n", bins(i), dos(i), doserr(i));
  }
  fclose(dosf);
  final_hs.save("final_hs.dat", raw_ascii);
}
