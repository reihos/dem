#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

// grain structure with the following members:radius, mass density, x of center,
// y of center, translational velocity and acceleration in 2 directions,
// rotation angle, rotational velocity, rotational acceleration, mass, moment of
// inertia, and gravity
struct grain {
  double r, rho, xc, yc, vx, vy, ax, ay, theta, omega, alpha;
  double m = rho * M_PI * r * r;
  double inertia = 0.5 * m * r * r;
  double gx, gy;
  bool fix;
};

const int ngrains = 2;
struct grain grains[ngrains];

struct contact {
  double kn;  // normal stiffness
  double kt;  // tangent stiffness
  // double damp_ratio;    //damping ratio
  double cor;       // coefficient of restitution
  double fric;      // friction coefficient (mu)
  double rollfric;  // rolling friiction coefficient (mu_r)
  double gamma =
      -log(cor) /
      sqrt(M_PI * M_PI + log(cor) * log(cor));  // needed for calculating cn
                                                // based on cor (Cleary,2000)
};

struct contact grain_grain;
double dt = 1e-6;  // timestep

// function to calculate the forces between particles jm and js and update the
// acceleration of jm (master particle) in 2 directions.
void update_acceleration(int jm, int js) {

  double distance;
  double delta;
  double normal[2];
  double tangent[2];
  double rel_vel_normal;
  double rel_vel_tangent;
  double rel_theta;
  double r_avg;
  double cn;
  double force_normal;
  double force_tangent;
  double rolling_moment;
  double force[2];

  grains[jm].ax = grains[jm].gx;
  grains[jm].ay = grains[jm].gy;

  distance = sqrt(pow(grains[jm].xc - grains[js].xc, 2) +
                  pow(grains[jm].yc - grains[js].yc, 2));
  delta = distance - grains[jm].r - grains[js].r;
  if (delta < 0.) {
    // calculating the normal and tangent unit vectors
    normal[0] = (grains[jm].xc - grains[js].xc) / distance;
    normal[1] = (grains[jm].yc - grains[js].yc) / distance;
    tangent[0] = -normal[1];
    tangent[1] = normal[0];
    // calculating the normal and tangent relative velocities
    rel_vel_normal = (grains[jm].vx - grains[js].vx) * normal[0] +
                     (grains[jm].vy - grains[js].vy) * normal[1];
    rel_vel_tangent = (grains[jm].vx - grains[js].vx) * tangent[0] +
                      (grains[jm].vy - grains[js].vy) * tangent[1];
    // calculating the relative rotation angle and average radius
    rel_theta = (grains[jm].omega + grains[js].omega) * dt;
    r_avg = (grains[jm].r + grains[js].r) / 2;
    // calculating normal and tangent forces and rolling moment
    cn = 2 * grain_grain.gamma *
         sqrt(grain_grain.kn * (grains[jm].m * grains[js].m) /
              (grains[jm].m + grains[js].m));  //(Cleary,2000)
    force_normal = -grain_grain.kn * delta - cn * rel_vel_normal;
    // if (force_normal<0){force_normal=0;}
    force_tangent = -grain_grain.kt * rel_vel_tangent * dt;
    if (abs(force_tangent) > grain_grain.fric * force_normal) {
      if (force_tangent > 0) {
        force_tangent = grain_grain.fric * force_normal;
      }
      if (force_tangent < 0) {
        force_tangent = -grain_grain.fric * force_normal;
      }
    }
    rolling_moment = -grain_grain.kt * r_avg * r_avg * rel_theta;
    if (abs(rolling_moment) > grain_grain.rollfric * force_normal * r_avg) {
      if (rolling_moment > 0) {
        rolling_moment = grain_grain.rollfric * force_normal * r_avg;
      }
      if (rolling_moment < 0) {
        rolling_moment = -grain_grain.rollfric * force_normal * r_avg;
      }
    }
    // calculating forces in x and y direction
    force[0] = force_normal * normal[0] + force_tangent * tangent[0];
    force[1] = force_normal * normal[1] + force_tangent * tangent[1];
    // updatig accerelations
    grains[jm].ax = force[0] / grains[jm].m;
    grains[jm].ay = force[1] / grains[jm].m;
    grains[jm].alpha =
        (-force_tangent * grains[jm].r + rolling_moment) / grains[jm].inertia;
    // std::cout<<delta<<"\t"<<force_normal<<"\t"<<force[1]<<"\t"<<grains[jm].ay<<std::endl;
    // adding the assigned gravity to the accelerations
    grains[jm].ax += grains[jm].gx;
    grains[jm].ay += grains[jm].gy;
  }
}

// function to update the location and velocity of particle jm
void update_kinematics(int jm) {
  if (grains[jm].fix == false) {
    grains[jm].xc += grains[jm].vx * dt + 0.5 * grains[jm].ax * dt * dt;
    grains[jm].yc += grains[jm].vy * dt + 0.5 * grains[jm].ay * dt * dt;
    grains[jm].vx += grains[jm].ax * dt;
    grains[jm].vy += grains[jm].ay * dt;
  }
}

// function to calculate the coordinates of npoints on the perimeter of all the
// grains and write to file for gnuplotting
void write_for_plotting(int i) {
  std::ofstream file;
  int npoints = 32;
  double x_perimeter[npoints];
  double y_perimeter[npoints];
  double angle = 2 * M_PI / npoints;

  file.open("stats" + std::to_string(i) + ".dat");
  for (int j = 0; j < ngrains; ++j) {
    file << "grain=" + std::to_string(j) << std::endl;
    for (int np = 0; np < npoints; ++np) {
      x_perimeter[np] = grains[j].xc + grains[j].r * cos(angle * np);
      y_perimeter[np] = grains[j].yc + grains[j].r * sin(angle * np);
      file << x_perimeter[np] << "\t" << y_perimeter[np] << std::endl;
    }
    file << "\n\n";
  }
  file.close();
}

// function to calculate the coordinates of npoints on the perimeter of all the grains and write a PostScript 
void write_ps(int step, int max_steps){
  std::ofstream file;
  std::stringstream file_name;
  int npoints=32;
  double x_perimeter;
  double y_perimeter;
  double angle=2*M_PI/npoints;
  
  //file_name.str(std::string());
  file_name<<"stats";
  file_name.fill('0');
  int digits=log10(max_steps)+1;
  file_name.width(digits);
  file_name<<step;
  file_name<<".ps";
  file.open(file_name.str());
  file<<"%!PS \n";
  file<<"%%BoundingBox: -200 -100 200 300\n"; 
  for (int j=0; j<ngrains;++j){
    for (int np=0; np<npoints;++np){
      x_perimeter=grains[j].xc+grains[j].r*cos(angle*np);
      y_perimeter=grains[j].yc+grains[j].r*sin(angle*np);
      if (np==0){
        file<<x_perimeter*1000<<" "<<y_perimeter*1000<<" "<<"newpath moveto ";
      }
      else{
        file<<x_perimeter*1000<<" "<<y_perimeter*1000<<" "<<"lineto ";
      }
    }
    file<<"closepath gsave ";
    file<<"0.9 0.7 0.5 setrgbcolor ";
    file<<"fill grestore ";
    file<<"0.8 0.4 0.1 setrgbcolor ";
    file<<"stroke\n";
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
void fix_grain(int j, const std::string& fix_type){
  if (fix_type=="all"){
   grains[j].fix = true;
  }
}

// function to calculate the time step for iterations
// void calculate_timestep(){
//  double
//  dt_critical=2*(sqrt(1+damp_ratio*damp_ratio)-damp_ratio)/sqrt(kn/grains[0].m);
//  dt=0.1*dt_critical;
//  std::cout<<"time step = "<<dt<<std::endl;
//}

int main() {

  grains[0] = {0.025, 1e8, 0., -0.025};
  grains[1] = {0.025, 33.33, 0., 0.2};
  set_gravity(0., -9.81);
  fix_grain(0, "all");
  grain_grain = {1e4, 0., 0.7, 0., 0.};
  // calculate_timestep();
  double t = 1.;
  int iwrite = 1000;
  int nsteps = int(t / dt);

  for (int i = 0; i < nsteps; ++i) {
    for (int jm = 0; jm < ngrains; ++jm) {
      for (int js = 0; js < ngrains; ++js) {
        if (jm != js) {
          update_acceleration(jm, js);
        }
      }
    }
    for (int jm = 0; jm < ngrains; ++jm) {
      update_kinematics(jm);
    }
    if (i%iwrite==0){
      write_ps(i, nsteps);
    }
    
  //  std::ofstream file;
  //  file.open("cpp_dem_results.txt", std::fstream::app);
  //  file <<i*dt<<"\t"<<grains[1].yc<<std::endl;
  //  file.close();

  }
}
