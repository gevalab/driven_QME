#ifndef CONSTANT_H
#define CONSTANT_H

#include <iostream>
#include <armadillo>
#include <complex>

using namespace std;

const double h_ = 1.0;         // reduced plank constant, in atomic unit
const double pi = 3.14159265358979323846; //pi

class Hamiltonian
{
    public:
    static double  w;          // frequency of the field
    static double  w01;        // frequency of the 2-level system
    static double  epsilon;    // system-field coupling
    static double  Lambda;     // system-bath coupling
    static int     N;          // number of bath modes (for GOA model: primary = 1, secondary = N-1)
    static double  Omega;      // frequency of the primary mode
    static double  y0;         // shift of the primary mode
    static double  wc;         // cutoff frequency of Ohmic spectral density
    static double  eta;        // coupling strength of Ohmic spectral density
    static double  beta;       // 1/kBT
    static int     N_point;    // number of integral points 
    arma::vec      w_;         // frequency of the normalized bath modes
    arma::vec      Req;        // equilibrium position of the normalized bath modes
    void print() const;
    int Ohmic_mode_generator();
    std::complex<double> V10_V01(double t, double tau) const;
    std::complex<double> V01_V10(double t, double tau) const;
};

class timer
{
    public:
    double T;            // total simulation time
    double dt;           // timestep
    double t;            // current time, default=0
    timer(double=1,double=0.1,double=0);
    void print() const;
};

class state
{
    public:
    arma::cx_mat sigma = arma::cx_mat(2,2,arma::fill::zeros);  // density matrix
    void print() const;
    int check_normal();
    class state evolution(double dt, arma::cx_mat dsigma_dt) const;
    class state & operator= (const class state & A);
    double sigma_x();
    double sigma_y();
    double sigma_z();
    double sigma_tr();
};

#endif