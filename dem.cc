#include <iostream>
#include <cmath>
#include <fstream>

// grain structure with the following members:radius, mass density, x of center, y of center, translational velocity and acceleration in 2 directions, rotation angle, rotational velocity, rotational acceleration, and mass
struct grain{
  double r,rho,xc,yc,vx,vy,ax,ay,theta,omega,alpha;
  double m=rho*M_PI*r*r;
};

const int ngrains=2;
struct grain grains[ngrains];
const double kn=1e6;

// function to calculate the forces between particles jm and js and update the acceleration of jm (master particle) in 2 directions.
void update_acceleration(int jm,int js){

  double distance;
  double delta;
  double normal[2];
  double force_normal;
  double force[2];

  distance = sqrt(pow(grains[jm].xc-grains[js].xc,2)+pow(grains[jm].yc-grains[js].yc,2)); 
  delta = distance - grains[jm].r -grains[js].r;
  if (delta<0.){
    force_normal=-kn*delta;
    normal[0]=(grains[jm].xc-grains[js].xc)/distance;
    normal[1]=(grains[jm].yc-grains[js].yc)/distance;
    force[0]=force_normal*normal[0];
    force[1]=force_normal*normal[1];
    grains[jm].ax=force[0]/grains[jm].m;
    grains[jm].ay=force[1]/grains[jm].m;
  }
}

// function to update the location and velocity of particle jm using dt as timestep
void update_kinematics(int jm,double dt){
  grains[jm].xc+=grains[jm].vx*dt+0.5*grains[jm].ax*dt*dt;
  grains[jm].yc+=grains[jm].vy*dt+0.5*grains[jm].ay*dt*dt;
  grains[jm].vx+=grains[jm].ax*dt;
  grains[jm].vy+=grains[jm].ay*dt;
}

// function to caluclate the coordnates of npoints on the perimeter of all the grains and write to file for plotting
void write_for_plotting(int i){
  std::ofstream file;
  int npoints=32;
  double x_perimeter[npoints];
  double y_perimeter[npoints];
  double angle=2*M_PI/npoints;
  
  file.open("stats" + std::to_string(i) + ".dat");
  for (int j=0; j<ngrains;++j){
    file <<"grain="+ std::to_string(j)<<std::endl;
    for (int np=0; np<npoints;++np){
      x_perimeter[np]=grains[j].xc+grains[j].r*cos(angle*np);
      y_perimeter[np]=grains[j].yc+grains[j].r*sin(angle*np);
      file <<x_perimeter[np]<<"\t"<<y_perimeter[np]<<std::endl;
    }
    file<<"\n\n";
  }
  file.close();
}


int main(){
  
  grains[0]= {1,2650.,2.,5.,-5};
  grains[1]= {1,2650.,-2,5.,5};
  double dt=1e-4;
  double t=1.;
  int iwrite=100;
  double nsteps=int(t/dt);

  for (int i=0; i<nsteps;++i){
    //std::cout<<i<<std::endl;
    for (int jm=0; jm<ngrains; ++jm){
      for (int js=0; js<ngrains; ++js){
        if (jm!=js){
        update_acceleration(jm,js);
        }
      }
    }
    for (int jm=0; jm<ngrains; ++jm){
      update_kinematics(jm,dt);
    }
    if (i%iwrite==0){
      write_for_plotting(i);
    }
  }
}

