/////////////////////////////
//    driven_GQME 1.0      //
// edit by Zongwei Huang   //
//       2022.12.6         //
/////////////////////////////
// This is a population-only non-condon version for the driven system QME

#include <iostream>
#include "constant.h"
#include "input.h"
#include "initialize.h"

using namespace std;

int main(){
    class Hamiltonian Hami_sim;
    class timer time_sim;
    class state stat_sim;
    Hami_sim.print();
    time_sim.print();
    input(Hami_sim,time_sim,stat_sim);
    initialize(Hami_sim,time_sim,stat_sim);//hhhh
    //evolution();
    //analysis();
    //output();
    
    return 0;
}
