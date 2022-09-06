//
// Created by tkaravasilhs on 11/8/2022.
//

#include "Variables.h"

double xmin =0, xmax = 5, ymin = 0, ymax = 1, mu = 0.05, rho = 1.0, numb_iter = 20, numb_timst = 1200, dt = 0.01,\
       omegau = 0.3, omegav = 0.7, omegap = 0.7, beta = 0.95, uInlet = 1.0, dx, dy,t=0;

int  tt = 0;

Eigen::MatrixXd  u,v,p,uold,vold,pp, \
        uu,vv,pcorner,apu,apv,app,source,ae,aw,an,as, \
        u0, v0 ;

void initialize(){
    u.resize(nx+2,ny+2);
    v.resize(nx+2,ny+2);
    p.resize(nx+2,ny+2);
    uold.resize(nx+2,ny+2);
    vold.resize(nx+2,ny+2);
    uu.resize(nx+2,ny+2);
    vv.resize(nx+2,ny+2);
    pp.resize(nx+2,ny+2);
    app.resize(nx+2,ny+2);
    apu.resize(nx+2,ny+2);
    apv.resize(nx+2,ny+2);
    ae.resize(nx+2,ny+2);
    aw.resize(nx+2,ny+2);
    an.resize(nx+2,ny+2);
    as.resize(nx+2,ny+2);
    u0.resize(nx+2,ny+2);
    v0.resize(nx+2,ny+2);
    source.resize(nx+2,ny+2);
    pcorner.resize(nx+2,ny+2);


    //initialize variables for channel flow
    u.setConstant(0), v.setConstant(0), p.setConstant(0);

    u(Eigen::all,0).setConstant(0);
    u(Eigen::all, Eigen::last).setConstant(0);

    u(1,Eigen::all).setConstant(uInlet);
    u(Eigen::last, Eigen::all) = u(1, Eigen::all);

    //compute mesh spacing
    dx=xmax/double(nx);
    dy=ymax/double(ny);

    //initialize uold, vold
    uold = u, vold =v;
    //initialize old time step results
    u0=u, v0=v;

    //     initialize solver coefficients
    app.setConstant(1);
    apu.setConstant(1);
    apv.setConstant(1);
    ae.setConstant(0);
    aw.setConstant(0);
    an.setConstant(0);
    as.setConstant(0);

}