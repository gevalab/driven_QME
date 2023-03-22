/////////////////////////////
//    driven_GQME 1.0      //
// edit by Zongwei Huang   //
//       2022.12.6         //
/////////////////////////////
// This is a population-only non-condon version for the driven system QME

#include <iostream>
#include "constant.h"
#include "ioput.h"
#include "initialize.h"
#include "evolution.h"

int main(){
    class Hamiltonian Hami_sim;
    class timer time_sim;
    class state stat_sim;
    int error_flag = 0;
    input(Hami_sim,time_sim,stat_sim);
    error_flag = initialize(Hami_sim,time_sim,stat_sim);
    if(error_flag == 1){return 0;}
    evolution(Hami_sim,time_sim,stat_sim);
    //analysis();
    //output();
    
    return 0;
}
