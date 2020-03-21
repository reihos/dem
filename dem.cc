#include <cmath>
#include <iostream>
#include <fstream>

// grain structure with the following members:radius, mass density, x of center, y of center, translational velocity and acceleration in 2 directions, rotation angle, rotational velocity, rotational acceleration, mass, moment of inertia, and gravity
struct grain{
  double r,rho,xc,yc,vx,vy,ax,ay,theta,omega,alpha;
  double m=rho*M_PI*r*r;
  double inertia=0.5*m*r*r;
  double gx,gy;
  bool fix;
};

const int ngrains=2;
struct grain grains[ngrains];
double kn=1e8;             //normal stiffness
double kt=1e6;             //tangent stiffness
double damp_ratio=0.05;    //damping ratio
double fric=0.5;           // friction coefficient (mu)
double rollfric=0.5;       // rolling friiction coefficient (mu_r)
double c_critical=2*sqrt(kn*grains[0].m); // critical viscosity
double cn=damp_ratio*c_critical; // normal viscosity
double dt;                 // timestep

// function to calculate the forces between particles jm and js and update the acceleration of jm (master particle) in 2 directions.
void update_acceleration(int jm,int js){

  double distance;
  double delta;
  double normal[2];
  double tangent[2];
  double rel_vel_normal;
  double rel_vel_tangent;
  double rel_theta;
  double r_avg;
  double force_normal;
  double force_tangent;
  double rolling_moment;
  double force[2];

  distance = sqrt(pow(grains[jm].xc-grains[js].xc,2)+pow(grains[jm].yc-grains[js].yc,2)); 
  delta = distance - grains[jm].r -grains[js].r;
  if (delta<0.){
    // caluclating the normal and tangent unit vectors
    normal[0]=(grains[jm].xc-grains[js].xc)/distance;
    normal[1]=(grains[jm].yc-grains[js].yc)/distance;
    tangent[0]=-normal[1];
    tangent[1]=normal[0];
    // calculating the normal and tangent relative velocities
    rel_vel_normal=(grains[jm].vx-grains[js].vx)*normal[0]+(grains[jm].vy-grains[js].vy)*normal[1];
    rel_vel_tangent=(grains[jm].vx-grains[js].vx)*tangent[0]+(grains[jm].vy-grains[js].vy)*tangent[1];
    // caluclating the relative rotation angle and average radius
    rel_theta=(grains[jm].omega+grains[js].omega)*dt;
    r_avg=(grains[jm].r+grains[js].r)/2;
    // caluclating normal and tangent forces and rolling moment
    force_normal=-kn*delta-cn*rel_vel_normal;
    if (force_normal<0){force_normal=0;}
    force_tangent=-kt*rel_vel_tangent*dt;
    if (abs(force_tangent)>fric*force_normal){
      if (force_tangent>0){force_tangent=fric*force_normal;}
      if (force_tangent<0){force_tangent=-fric*force_normal;}
    }
    rolling_moment=-kt*r_avg*r_avg*rel_theta;
    if (abs(rolling_moment)>rollfric*force_normal*r_avg){
      if (rolling_moment>0){rolling_moment=rollfric*force_normal*r_avg;}
      if (rolling_moment<0){rolling_moment=-rollfric*force_normal*r_avg;}
    }
    // caluclating forces in x and y direction
    force[0]=force_normal*normal[0]+force_tangent*tangent[0];
    force[1]=force_normal*normal[1]+force_tangent*tangent[1];
    // updatig accerelations
    grains[jm].ax=force[0]/grains[jm].m;
    grains[jm].ay=force[1]/grains[jm].m;
    grains[jm].alpha=(-force_tangent*grains[jm].r+rolling_moment)/grains[jm].inertia;
  }
  // adding the assigned gravity to the accelerations
  grains[jm].ax+=grains[jm].gx;
  grains[jm].ay+=grains[jm].gy;
}

// function to update the location and velocity of particle jm
void update_kinematics(int jm){
  if (grains[jm].fix==false){
    grains[jm].xc+=grains[jm].vx*dt+0.5*grains[jm].ax*dt*dt;
    grains[jm].yc+=grains[jm].vy*dt+0.5*grains[jm].ay*dt*dt;
    grains[jm].vx+=grains[jm].ax*dt;
    grains[jm].vy+=grains[jm].ay*dt;
  }
}

// function to caluclate the coordinates of npoints on the perimeter of all the grains and write to file for plotting
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

// function to set the gravity of the grains in 2 directions to the assigned values
void set_gravity(double gx, double gy){
  for (int j=0; j<ngrains;++j){
    grains[j].gx=gx;
    grains[j].gy=gy;
  }
}

// function to fix the velocity or position of grain j (for now I have only coded "all" which fixes both velocity and acceleration
void fix_grain(int j, std::string fix_type){
  if (fix_type=="all"){
   grains[j].fix = true;
  }
}

// function to calculate the time step for iterations
void calculate_timestep(){
  double dt_critical=2*(sqrt(1+damp_ratio*damp_ratio)-damp_ratio)/sqrt(kn/grains[0].m);
  dt=0.1*dt_critical;
  std::cout<<"time step = "<<dt<<std::endl;
}

int main(){
  
  grains[0]= {1,2650.,0.,-1.};
  grains[1]= {1,2650.,0.,5.};
  set_gravity(0.,-10);
  fix_grain(0,"all");
  calculate_timestep();
  double t=1.;
  int iwrite=100;
  int nsteps=int(t/dt);

  for (int i=0; i<nsteps;++i){
    for (int jm=0; jm<ngrains; ++jm){
      for (int js=0; js<ngrains; ++js){
        if (jm!=js){
        update_acceleration(jm,js);
        }
      }
    }
    for (int jm=0; jm<ngrains; ++jm){
      update_kinematics(jm);
    }
    if (i%iwrite==0){
      write_for_plotting(i);
    }
  }
}
