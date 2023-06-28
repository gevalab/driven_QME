#ifndef EVOLUTION_H
#define EVOLUTION_H
#include "constant.h"

int evolution(class Hamiltonian & Hami_sim, class timer & time_sim, class state & stat_sim);
arma::cx_mat dsigma_dt(const class Hamiltonian & Hami_sim, double t, const class state & stat_sim);
arma::cx_mat k_calculation_full(const class Hamiltonian & Hami_sim, double t, int N);
arma::cx_mat I_calculation_full(const class Hamiltonian & Hami_sim, double t);
std::complex<double> I_calculation(const class Hamiltonian & Hami_sim, int j, int k, double t);
std::complex<double> k_calculation(const class Hamiltonian & Hami_sim, int j, int k, double t, int N);
int steady_state(class Hamiltonian & Hami_sim, class timer & time_sim);

#endif