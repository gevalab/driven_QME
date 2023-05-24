#include "initialize.h"

int initialize(class Hamiltonian & Hami_sim, class timer & time_sim, class state & stat_sim){
    cout << "initialize start" << endl;
    // generate the Ohmic bath modes
    Hami_sim.Ohmic_mode_generator();
    cout << "Ohmic bath modes generated!" << endl;
    // check the densitry matrix
    if(stat_sim.check_normal() == 1){
        cout << "density matrix check completed!" << endl;
    }else{
        cout << "error: unnormalized density matrix sigma!" << endl;
        stat_sim.print();
        return 1;
    }
    cout << "initialize over" << endl;
    return 0;
}

int check_correlation_function(class Hamiltonian & Hami_sim){
    double t_check = 0;
    double tau_check = 1;
    cout << "V1001(" << t_check << "," << tau_check << ") = " << Hami_sim.V10_V01(t_check,tau_check) << endl;
    cout << "V0110(" << t_check << "," << tau_check << ") = " << Hami_sim.V01_V10(t_check,tau_check) << endl;
    return 0;
}