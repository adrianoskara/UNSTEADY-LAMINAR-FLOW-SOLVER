//
// Created by tkaravasilhs on 11/8/2022.
//
#include "iostream"
#include "fstream"
#include <iomanip>
#include "eigen-3.4.0/Eigen/Dense"
#include "Variables.h"

// Function that deals with the u-momentum equations
void u_momentum() {
    double mdote, mdotw, mdotn, mdots, res, min, mout;

    Eigen::Matrix<double, nx + 2, ny + 2> dcu, dcc;
    dcc.setConstant(0);
    dcu.setConstant(0);

    //set u_momentum equation
    for (int i = 2; i < nx + 1; i++) {
        for (int j = 1; j < ny + 1; j++) {
            mdote = 0.5 * rho * dy * (uold(i + 1, j) + uold(i, j));
            mdotw = 0.5 * rho * dy * (uold(i - 1, j) + uold(i, j));
            mdotn = 0.5 * rho * dx * (vold(i, j + 1) + vold(i - 1, j + 1));
            mdots = 0.5 * rho * dx * (vold(i, j) + vold(i - 1, j));


            ae(i, j) = std::max(-mdote, 0.0);
            aw(i, j) = std::max(mdotw, 0.0);
            an(i, j) = std::max(-mdotn, 0.0);
            as(i, j) = std::max(mdots, 0.0);

            apu(i, j) = ae(i, j) + aw(i, j) + an(i, j) + as(i, j);


            dcu(i, j) = -(ae(i, j) * u(i + 1, j) + aw(i, j) * u(i - 1, j) + \
                                an(i, j) * u(i, j + 1) + as(i, j) * u(i, j - 1)) + apu(i, j) * u(i, j);

            dcc(i, j) = 0.5 * (mdote * (u(i + 1, j) + u(i, j)) - mdotw * (u(i, j) + u(i - 1, j)) + \
                              mdotn * (u(i, j + 1) + u(i, j)) - mdots * (u(i, j) + u(i, j - 1)));


            ae(i, j) = ae(i, j) + mu * dx / double((dy));
            aw(i, j) = aw(i, j) + mu * dx / double((dy));
            an(i, j) = an(i, j) + mu * dx / double((dy));
            as(i, j) = as(i, j) + mu * dx / double((dy));

        }
    }


    //overwrite boundary coefficients along north/south walls with half cell (dy) size
    // first order only


    for (int i = 2; i < nx + 1; i++) {
        mdotn = 0.5 * rho * dx * (vold(i, ny + 1) + vold(i - 1, ny + 1));
        mdots = 0.5 * rho * dx * (vold(i, 1) + vold(i - 1, 1));
        an(i, ny) = std::max(-mdotn, 0.0) + mu * dx / double((dy / 2.));
        as(i, 1) = std::max(mdots, 0.0) + mu * dx / double((dy / 2.));
    }

    apu = ae + aw + an + as;
    apu.array() += +rho * dx * dy / double(dt);
    apu = apu / omegau;

    //square cylinder addition
    apu(Eigen::seq(20,31),Eigen::seq(10,15).setConstant(1.e30);

    //iterate x-momentum equations
    for (int l = 1; l < 11; l++) {
        for (int i = 2; i < nx + 1; i++) {
            for (int j = 1; j < ny + 1; j++) {
                u(i, j) = (1. - omegau) * uold(i, j) + 1. / double(apu(i, j)) * (ae(i, j) * u(i + 1, j) + \
                                  aw(i, j) * u(i - 1, j) + an(i, j) * u(i, j + 1) + as(i, j) * u(i, j - 1) + \
                                  dy * (p(i - 1, j) - p(i, j)) - beta * (dcc(i, j) - dcu(i, j)) +
                                                                                 rho * dx * dy * u0(i, j) / double(dt));
            }
        }
        u(nx + 1, Eigen::all) = u(nx, Eigen::all); //du/dx=0 at outlet
    }
    //In order to ensure mass conservation
    min = 1.0;
    mout = 0.0;

    for (int i = 1; i < ny+1; i++) {
        mout = mout + dy * rho * u(nx + 1, i);
    }
    u(nx + 1, Eigen::all) = u(nx + 1, Eigen::all) * min / double(mout);

    // Computing x-momentum residual
    res = 0.0;
    for (int i = 2; i < nx + 1; i++) {
        for (int j = 1; j < ny + 1; j++) {
            res = res + pow((u(i, j) - uold(i, j)), 2.0);
        }
    }
    res = sqrt(res);
    std::cout<<"x-momentum residual = "<<res<<std::endl;
}

void v_momentum(){

    double mdote, mdotw, mdotn, mdots, res, min, mout;

    Eigen::Matrix<double, nx + 2, ny + 2> dcu, dcc;
    dcc.setConstant(0);
    dcu.setConstant(0);

    //set v_momentum equation
    for (int i = 1; i < nx + 1; i++) {
        for (int j = 2; j < ny + 1; j++) {
            mdote=0.5*rho*dy*(uold(i+1,j) + uold(i+1,j-1));
            mdotw=0.5*rho*dy*(uold(i,j-1) + uold(i,j));
            mdotn=0.5*rho*dx*(vold(i,j+1) + vold(i,j));
            mdots=0.5*rho*dx*(vold(i,j-1) + vold(i,j));
            ae(i,j)=std::max(-mdote,0.0);
            aw(i,j)=std::max( mdotw,0.0);
            an(i,j)=std::max(-mdotn,0.0);
            as(i,j)=std::max( mdots,0.0);
            apv(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j);

            dcu(i,j)=-(ae(i,j)*v(i+1,j)+aw(i,j)*v(i-1,j) + an(i,j)*v(i,j+1)+ \
                                as(i,j)*v(i,j-1)) +apv(i,j)*v(i,j);
            dcc(i,j)=0.5*( mdote*(v(i+1,j)+v(i,j)) - mdotw*(v(i,j)+v(i-1,j)) + \
                              mdotn*(v(i,j+1)+v(i,j)) - mdots*(v(i,j)+v(i,j-1)) );

            ae(i,j)=ae(i,j)+mu*dy/double(dx);
            aw(i,j)=aw(i,j)+mu*dy/double(dx);
            an(i,j)=an(i,j)+mu*dx/double(dy);
            as(i,j)=as(i,j)+mu*dx/double(dy);
        }
    }


    //overwrite boundary coefficients along east/west boundaries due to half cell (dx) size
    //first order only

    for (int j=2;j<ny+1;j++){
        mdote=0.5*rho*dy*(uold(nx+1,j) + uold(nx+1,j-1));
        mdotw=0.5*rho*dy*(uold(1,j) +  uold(1,j-1));
        ae(nx,j)=std::max(-mdote,0.0)+mu*dy/double((dx/2.0));
        aw(1,j)=std::max( mdotw,0.0)+mu*dy/double((dx/2.0));
        //      ae(nx,j)=-mdote/2. + mu*dy/(dx/2.0);
        //       aw(1,j)=  mdotw/2. + mu*dy/(dx/2.0);
    }


    apv=ae+aw+an+as;
    apv.array() += +rho * dx * dy /double(dt);
    apv=apv/omegav;


    //square cylinder addition
    apv(Eigen::seq(20,30),Eigen::seq(10,16).setConstant(1.e30);

    //iterate y-momentum equations
    for (int l = 1; l < 11; l++) {
        for (int i = 1; i < nx + 1; i++) {
            for (int j = 2; j < ny + 1; j++) {
                v(i,j)=(1.-omegav)*vold(i,j)+1./double(apv(i,j))*(ae(i,j)*v(i+1,j)+\
                                aw(i,j)*v(i-1,j)+an(i,j)*v(i,j+1)+as(i,j)*v(i,j-1)+ \
                                dx*(p(i,j-1)-p(i,j)) - beta*(dcc(i,j)-dcu(i,j)) + rho*dx*dy*v0(i,j)/double(dt));
            }
        }
        v(nx + 1, Eigen::all) = v(nx, Eigen::all); //dv/dx=0 at outlet

    }

    //calculate residual
    res = 0.0;
    for (int i=1;i<nx+1;i++){
        for (int j=2;j<ny+1;j++){
            res = res + pow((v(i,j)-vold(i,j)),2);
        }
    }
    res = sqrt(res);
    std::cout<<"y-momentum residual = "<<res<<std::endl;


}

//PRESSURE CORRECTION ROUTINE
void pressure(){

    double mdote, mdotw, mdotn, mdots, res, min, mout;

    Eigen::Matrix<double, nx + 2, ny + 2> dcu, dcc;



    //set coefficients
    for (int i = 1;i<nx+1;i++){
        for (int j = 1;j<ny+1;j++){
            ae(i,j)=rho*pow(dy,2)/double(apu(i+1,j));
            aw(i,j)=rho*pow(dy,2)/double(apu(i,j));
            an(i,j)=rho*pow(dx,2)/double(apv(i,j+1));
            as(i,j)=rho*pow(dx,2)/double(apv(i,j));
        }
    }

    //set boundary values for coefficients
    ae(nx,Eigen::all).setConstant(0.0);
    aw(1,Eigen::all).setConstant(0.0);
    an(Eigen::all,ny).setConstant(0.0);
    as(Eigen::all,1).setConstant(0.0);

    app=ae+aw+an+as;

    pp.setConstant(0);

    //     compute the mass source term before corrections
    source.setConstant(0);
    for (int i=1;i<nx+1;i++){
        for (int j=1;j<ny+1;j++){
            source(i,j)=rho*dy*(u(i+1,j)-u(i,j)) +rho*dx*(v(i,j+1)-v(i,j));
        }
    }


    // SOR iterations to solve for pressure correction ppl
    for (int l=1;l<100;l++){
        for (int j=1;j<ny+1;j++){
            for (int i=1;i<nx+1;i++){
                pp(i,j)=pp(i,j)+1.7/double(app(i,j))*(ae(i,j)*pp(i+1,j)+ \
                        aw(i,j)*pp(i-1,j)+an(i,j)*pp(i,j+1)+as(i,j)*pp(i,j-1)- \
                                source(i,j)-pp(i,j)*app(i,j));
            }
        }
    }

    //Apply corrections to pressure

    for (int i=1;i<nx+1;i++){
        for (int j=1;j<ny+1;j++){
            p(i,j)=p(i,j)+omegap*pp(i,j);
        }
    }


    //Apply corrections to u-component
    for (int i=2;i<nx+1;i++){
        for (int j=1;j<ny+1;j++) {
            u(i,j)=u(i,j)+dy/double(apu(i,j))*(pp(i-1,j)-pp(i,j));
        }
    }


    //Apply corrections to v-component
    for (int i=1;i<nx+1;i++){
        for (int j=2;j<ny+1;j++) {
            v(i,j)=v(i,j)+dx/apv(i,j)*(pp(i,j-1)-pp(i,j));
        }
    }


    //update velocity variables
    uold=u;
    vold=v;

}

void write(){
    double xx=0.0,yy=0.0, velmag = 0.0;
    std::ofstream fileOut;

    fileOut.open("results.csv."+std::to_string(tt));



    //interpolate velocities to scalar cell corners
    for (int i=2;i<nx+1;i++){
        for (int j=2;j<ny+1;j++) {
            uu(i,j)=0.5*(u(i,j-1)+u(i,j));
            vv(i,j)=0.5*(v(i-1,j)+v(i,j));
            pcorner(i,j)=0.25*(p(i-1,j-1)+p(i-1,j)+p(i,j-1)+p(i,j));
        }
    }

    // fill in east and west boundaries
    uu(1,Eigen::seq(2,ny))=0.5*(u(1,Eigen::seq(1,ny-1))+u(1,Eigen::seq(2,ny)));
    uu(nx+1,Eigen::seq(2,ny))=0.5*(u(nx+1,Eigen::seq(1,ny-1))+u(nx+1,Eigen::seq(2,ny)));
    vv(1,Eigen::seq(2,ny))=v(0,Eigen::seq(2,ny));
    vv(nx+1,Eigen::seq(2,ny))=v(nx+1,Eigen::seq(2,ny));
    pcorner(1,Eigen::seq(2,ny))=0.5*(p(1,Eigen::seq(1,ny-1))+p(1,Eigen::seq(2,ny)));
    pcorner(nx+1,Eigen::seq(2,ny))=0.5*(p(nx,Eigen::seq(1,ny-1))+p(nx,Eigen::seq(2,ny)));

   //    fill in north and south boundaries
    uu(Eigen::seq(2,nx),1)=u(Eigen::seq(2,nx),0);
    uu(Eigen::seq(2,nx),ny+1)=u(Eigen::seq(2,nx),ny+1);
    vv(Eigen::seq(2,nx),1)=0.5*(v(Eigen::seq(1,nx-1),1)+v(Eigen::seq(2,nx),1));
    vv(Eigen::seq(2,nx),ny+1)=0.5*(v(Eigen::seq(1,nx-1),ny+1)+v(Eigen::seq(2,nx),ny+1));
    pcorner(Eigen::seq(2,nx),1)=0.5*(p(Eigen::seq(1,nx-1),1)+p(Eigen::seq(2,nx),1));
    pcorner(Eigen::seq(2,nx),ny+1)=0.5*(p(Eigen::seq(1,nx-1),ny)+p(Eigen::seq(2,nx),ny));

    //fill in southwest corner
    uu(1,1)=0;
    vv(1,1)=0;
    pcorner(1,1)=p(1,1);

    //fill in southeast corner
    uu(nx+1,1)=0;
    vv(nx+1,1)=0;
    pcorner(nx+1,1)=p(nx,1);

    //fill in northeast corner
    uu(nx+1,ny+1)=0;
    vv(nx+1,ny+1)=0;
    pcorner(nx+1,ny+1)=p(nx,ny);

    // fill in northwest corner
    uu(1,ny+1)=0;
    vv(1,ny+1)=0;
    pcorner(1,ny+1)=p(1,ny);

    fileOut<<"x, y, z, pressure, velmag, vx, vy, vz"<<std::endl;

    //Square cylinder addition
    uu(Eigen::seq(20,31),Eigen::seq(10,15)).setConstant(0);
    vv(Eigen::seq(20,30),Eigen::seq(10,16)).setConstant(0);

    for (int i=1;i<nx+2;i++){
        xx = double(i-1)*dx;
        for (int j=1;j<ny+2;j++){
            yy = double(j-1)*dy;
            velmag=sqrt(pow(uu(i,j),2)+pow(vv(i,j),2));
            fileOut<<std::setprecision(4)<<xx<<","<<yy<<","<<0.0<<","<<pcorner(i,j)<<","<<velmag<<","<<uu(i,j)<<","<<vv(i,j)<<","<<0.0<<std::endl;
        }
    }
    fileOut.close();


}