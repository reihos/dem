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

const int ngrains = 1;
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

struct wall {
  double normal[2];  // coordinates of the normal unit vector of the wall
  double point[2];   // coordinates of a point on the wall
  double tangent[2] = {-normal[1],
                       normal[0]};  // coordinates of the tangent unit vector
  double tan_slope =
      tangent[1] / tangent[0];  // tangent of the slope of the wall
};

const int nwalls = 1;
struct wall walls[nwalls];

// function to initialize the acceleration of grain jm
void initialize_acceleration(int jm) {
  grains[jm].ax = grains[jm].gx;
  grains[jm].ay = grains[jm].gy;
  grains[jm].alpha = 0.;
}

// function to calculate the forces between particles jm and js and update the
// acceleration of jm (master particle) in 2 directions.
void update_acceleration_grain_grain(int jm, int js) {

  const double distance = sqrt(pow(grains[jm].xc - grains[js].xc, 2) +
                               pow(grains[jm].yc - grains[js].yc, 2));
  const double delta = distance - grains[jm].r - grains[js].r;
  if (delta < 0.) {
    // calculating the normal and tangent unit vectors
    double normal[2];
    double tangent[2];
    normal[0] = (grains[jm].xc - grains[js].xc) / distance;
    normal[1] = (grains[jm].yc - grains[js].yc) / distance;
    tangent[0] = -normal[1];
    tangent[1] = normal[0];
    // calculating the normal and tangent relative velocities
    const double rel_vel_normal = (grains[jm].vx - grains[js].vx) * normal[0] +
                                  (grains[jm].vy - grains[js].vy) * normal[1];
    const double rel_vel_tangent =
        (grains[jm].vx - grains[js].vx) * tangent[0] +
        (grains[jm].vy - grains[js].vy) * tangent[1] -
        grains[jm].r * grains[jm].omega - grains[js].r * grains[js].omega;
    // calculating the relative rotation angle and average radius
    const double rel_theta = (grains[jm].omega + grains[js].omega) * dt;
    const double r_avg = (grains[jm].r + grains[js].r) / 2;
    // calculating normal and tangent forces and rolling moment
    const double cn = 2 * grain_grain.gamma *
                      sqrt(grain_grain.kn * (grains[jm].m * grains[js].m) /
                           (grains[jm].m + grains[js].m));  //(Cleary,2000)
    double force_normal = -grain_grain.kn * delta - cn * rel_vel_normal;
    // if (force_normal<0){force_normal=0;}
    double force_tangent = -grain_grain.kt * rel_vel_tangent * dt;
    if (std::abs(force_tangent) > grain_grain.fric * force_normal) {
      if (force_tangent > 0) {
        force_tangent = grain_grain.fric * force_normal;
      }
      if (force_tangent < 0) {
        force_tangent = -grain_grain.fric * force_normal;
      }
    }
    double rolling_moment = -grain_grain.kt * r_avg * r_avg * rel_theta;
    if (std::abs(rolling_moment) >
        grain_grain.rollfric * force_normal * r_avg) {
      if (rolling_moment > 0) {
        rolling_moment = grain_grain.rollfric * force_normal * r_avg;
      }
      if (rolling_moment < 0) {
        rolling_moment = -grain_grain.rollfric * force_normal * r_avg;
      }
    }
    // calculating forces in x and y direction
    double force[2];
    force[0] = force_normal * normal[0] + force_tangent * tangent[0];
    force[1] = force_normal * normal[1] + force_tangent * tangent[1];
    // updatig accerelations
    grains[jm].ax += force[0] / grains[jm].m;
    grains[jm].ay += force[1] / grains[jm].m;
    grains[jm].alpha +=
        (-force_tangent * grains[jm].r + rolling_moment) / grains[jm].inertia;
  }
}

// function to calculate the forces between particle jm and wall w and update
// the acceleration of jm in 2 directions.
void update_acceleration_grain_wall(int jm, int w) {

  const double distance =
      (walls[w].tan_slope * grains[jm].xc - grains[jm].yc -
       walls[w].tan_slope * walls[w].point[0] + walls[w].point[1]) /
      sqrt(walls[w].tan_slope * walls[w].tan_slope + 1.);
  double delta;
  if (walls[w].tangent[0] > 0.) {
    delta = distance - grains[jm].r;
  } else if (walls[w].tangent[0] < 0.) {
    delta = -distance - grains[jm].r;
  }
  if (delta < 0.) {
    // calculating the normal and tangent relative velocities
    const double rel_vel_normal =
        grains[jm].vx * walls[w].normal[0] + grains[jm].vy * walls[w].normal[1];
    const double rel_vel_tangent = grains[jm].vx * walls[w].tangent[0] +
                                   grains[jm].vy * walls[w].tangent[1] -
                                   grains[jm].r * grains[jm].omega;
    // calculating the relative rotation angle and average radius
    const double rel_theta = grains[jm].omega * dt;
    const double r_avg = grains[jm].r;
    // calculating normal and tangent forces and rolling moment
    const double cn =
        2 * grain_grain.gamma * sqrt(grain_grain.kn * grains[jm].m);
    double force_normal = -grain_grain.kn * delta - cn * rel_vel_normal;
    // if (force_normal<0){force_normal=0;}
    double force_tangent = -grain_grain.kt * rel_vel_tangent * dt;
    if (std::abs(force_tangent) > grain_grain.fric * force_normal) {
      if (force_tangent > 0.) {
        force_tangent = grain_grain.fric * force_normal;
      }
      if (force_tangent < 0.) {
        force_tangent = -grain_grain.fric * force_normal;
      }
    }
    double rolling_moment = -grain_grain.kt * r_avg * r_avg * rel_theta;
    if (std::abs(rolling_moment) >
        grain_grain.rollfric * force_normal * r_avg) {
      if (rolling_moment > 0.) {
        rolling_moment = grain_grain.rollfric * force_normal * r_avg;
      }
      if (rolling_moment < 0.) {
        rolling_moment = -grain_grain.rollfric * force_normal * r_avg;
      }
    }
    // calculating forces in x and y direction
    double force[2];
    force[0] =
        force_normal * walls[w].normal[0] + force_tangent * walls[w].tangent[0];
    force[1] =
        force_normal * walls[w].normal[1] + force_tangent * walls[w].tangent[1];
    // updatig accerelations
    grains[jm].ax += force[0] / grains[jm].m;
    grains[jm].ay += force[1] / grains[jm].m;
    grains[jm].alpha +=
        (-force_tangent * grains[jm].r + rolling_moment) / grains[jm].inertia;
    std::cout << "normalized tangent force = "
              << force_tangent / (grains[jm].m * 9.81 * sin(45. / 180. * M_PI))
              << std::endl;
    std::cout << "normalized normal force = "
              << force_normal / (grains[jm].m * 9.81 * cos(45. / 180. * M_PI))
              << std::endl;
  m
}
// function to update the location and velocity of particle jm
void update_kinematics(int jm) {
  if (grains[jm].fix == false) {
    grains[jm].xc += grains[jm].vx * dt + 0.5 * grains[jm].ax * dt * dt;
    grains[jm].yc += grains[jm].vy * dt + 0.5 * grains[jm].ay * dt * dt;
    grains[jm].vx += grains[jm].ax * dt;
    grains[jm].vy += grains[jm].ay * dt;
    grains[jm].theta +=
        grains[jm].omega * dt + 0.5 * grains[jm].alpha * dt * dt;
    grains[jm].omega += grains[jm].alpha * dt;
  }
}

// function to calculate the coordinates of npoints on the perimeter of all the
// grains and write to file for gnuplotting
void write_for_plotting(int i) {
  std::ofstream file;
  int npoints = 32;
  double x_perimeter;
  double y_perimeter;
  double angle = 2 * M_PI / npoints;

  file.open("stats" + std::to_string(i) + ".dat");
  for (int j = 0; j < ngrains; ++j) {
    file << "grain=" + std::to_string(j) << std::endl;
    for (int np = 0; np < npoints; ++np) {
      x_perimeter = grains[j].xc + grains[j].r * cos(angle * np);
      y_perimeter = grains[j].yc + grains[j].r * sin(angle * np);
      file << x_perimeter << "\t" << y_perimeter << std::endl;
    }
    file << "\n\n";
  }
  file.close();
}

// function to calculate the coordinates of npoints on the perimeter of all the
// grains and write a PostScript
void write_ps(int step, int max_steps) {
  std::ofstream file;
  std::stringstream file_name;
  int npoints = 32;
  double x_perimeter;
  double y_perimeter;
  double angle = 2 * M_PI / npoints;

  // file_name.str(std::string());
  file_name << "stats";
  file_name.fill('0');
  int digits = log10(max_steps) + 1;
  file_name.width(digits);
  file_name << step;
  file_name << ".ps";
  file.open(file_name.str());
  file << "%!PS \n";
  file << "%%BoundingBox: -200 -100 200 300\n";
  for (int j = 0; j < ngrains; ++j) {
    for (int np = 0; np < npoints; ++np) {
      x_perimeter = grains[j].xc + grains[j].r * cos(angle * np);
      y_perimeter = grains[j].yc + grains[j].r * sin(angle * np);
      if (np == 0) {
        file << x_perimeter * 10 << " " << y_perimeter * 10 << " "
             << "newpath moveto ";
      } else {
        file << x_perimeter * 10 << " " << y_perimeter * 10 << " "
             << "lineto ";
      }
    }
    file << "closepath gsave ";
    file << "0.9 0.7 0.5 setrgbcolor ";
    file << "fill grestore ";
    file << "0.8 0.4 0.1 setrgbcolor ";
    file << "stroke\n";
  }
  file.close();
}

// function to set the gravity of the grains in 2 directions to the assigned
// values
void set_gravity(double gx, double gy) {
  for (int j = 0; j < ngrains; ++j) {
    grains[j].gx = gx;
    grains[j].gy = gy;
  }
}

// function to fix the velocity or position of grain j (for now I have only
// coded "all" which fixes both velocity and acceleration
void fix_grain(int j, const std::string& fix_type) {
  if (fix_type == "all") {
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

  double slope_angle = 45. / 180. * M_PI;
  walls[0] = {{-sin(slope_angle), cos(slope_angle)}, {0., 0.}};
  double x_start = 0.5;
  grains[0] = {0.01, 130.};
  grains[0].xc = x_start - grains[0].r * sin(slope_angle);
  grains[0].yc = x_start * tan(slope_angle) + grains[0].r * cos(slope_angle);
  set_gravity(0., -9.81);
  grain_grain = {1e5, 1e7, 0.5, tan(60. / 180. * M_PI), 1.};
  // calculate_timestep();
  double t = 0.5;
  int iwrite = 1000;
  int nsteps = int(t / dt);
  for (int i = 0; i < nsteps; ++i) {
    for (int jm = 0; jm < ngrains; ++jm) {
      initialize_acceleration(jm);
      for (int js = 0; js < ngrains; ++js) {
        if (jm != js) {
          update_acceleration_grain_grain(jm, js);
        }
      }
      for (int w = 0; w < nwalls; ++w) {
        update_acceleration_grain_wall(jm, w);
      }
    }
    for (int jm = 0; jm < ngrains; ++jm) {
      update_kinematics(jm);
    }
    if (i % iwrite == 0) {
      write_for_plotting(i);
    }

    //  std::ofstream file;
    //  file.open("cpp_dem_results.txt", std::fstream::app);
    //  file << i * dt << "\t" << grains[1].yc << std::endl;
    //  file.close();
  }
}
