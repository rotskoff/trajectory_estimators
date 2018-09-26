/* quartic.h
 
Copyright 2018 Grant M. Rotskoff

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */
#include "particles.h"

class quartic : public particles {
public:
  void initialize(int dim0, double gamma0);
  void compute_gradient();
  double compute_conf_energy();
  double compute_hamiltonian_energy();
  void reinject();

};

void quartic::initialize(int dim0, double gamma0)
{
  dim = dim0;
  gamma = gamma0;
  qs = zeros(dim);
  ps = zeros(dim);
  fs = zeros(dim);
}

void quartic::compute_gradient()
{
  fs = pow(qs,3);
}

double quartic::compute_conf_energy()
{
  return 0.25*accu(pow(qs,4));
}

double quartic::compute_hamiltonian_energy()
{
  return dot(ps,ps) + compute_conf_energy();
}

void quartic::reinject()
{
  qs = qmin + randu(size(qs))*(qmax-qmin);
  double h = compute_conf_energy();
  qs = qs / pow(h, 0.25) * pow(hmax, 0.25);
}
