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
