#ifndef EVOLUTION_H
#define EVOLUTION_H
#include "constant.h"

int evolution(class Hamiltonian & Hami_sim, class timer & time_sim, class state & stat_sim);
arma::cx_mat dsigma_dt(const class Hamiltonian & Hami_sim, double t, const class state & stat_sim);
std::complex<double> k_calculation(const class Hamiltonian & Hami_sim, int j, int k, double t, int N);

#endif