#ifndef INITIALIZE_H
#define INITIALIZE_H
#include "constant.h"

int initialize(class Hamiltonian & Hami_sim, class timer & time_sim, class state & stat_sim){
    cout << "initialize start" << endl;
    //generate the Ohmic bath modes
    Hami_sim.Ohmic_mode_generator();
    cout << "initialize over" << endl;
    return 0;
}


#endif