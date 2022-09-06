//
// Created by tkaravasilhs on 11/8/2022.
//


#ifndef UNSTEADY_VARIABLES_H
#define UNSTEADY_VARIABLES_H
#include "eigen-3.4.0/Eigen/Dense"
constexpr int nx=60 , ny = 30;


extern double xmin, xmax, ymin, ymax,mu, rho, numb_iter,\
              numb_timst, dt, omegau, omegav, omegap, beta,\
              uInlet, dx, dy, t;

extern int  tt;



extern Eigen::MatrixXd u,v,p,uold,vold,pp, \
        uu,vv,pcorner,apu,apv,app,source,ae,aw,an,as, \
        u0, v0;

void initialize();


#endif //UNSTEADY_VARIABLES_H
