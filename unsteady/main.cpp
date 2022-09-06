#include <iostream>
#include "Variables.h"
#include "schemes.h"
#include "progressbar-master/include/progressbar.hpp"

//This code models a low speed 2D flow around a square cylinder using a finite volume approach with
//explicit Pressure-Correction and explicit Forward Euler time stepping. The analysis focuses on the effects
//of the Reynolds number, as determined by the input velocity parameters, on the flow around the square.
//cylinder. The code is set to run both a 2D channel/pipe flow and a flow with a square cylinder, both
//of which could be validated to analytical solutions and literature results to ensure code accuracy.
int main() {

    // Parameter and variable initialize
    initialize();

    progressbar bar(100);
    for (int titer=1;titer<numb_timst;titer++){
        bar.update();
        // time step calculation
        t = t + dt; //time step for analysis
        tt = tt +1; //time step for file writing
        std::cout<<"time = "<<t<<std::endl;


        for (int iter=1;iter<numb_iter;iter++) {
            std::cout << "iteration = " << iter << std::endl;
            std::cout << std::endl;
            u_momentum();
            v_momentum();
            pressure();
            std::cout << std::endl;

        }
        u0 = u;
        v0 = v;

        // Writing output file every 10th timestep
        if ((tt%10) == 0){

            write();
        }

    }



}

