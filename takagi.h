#ifndef _TAKAGI_H_
#define _TAKAGI_H_

#define _USE_MATH_DEFINES
#include "vector3d.hpp"
#include <cmath>
#include <complex>
using namespace std;

typedef float angle_num_t;
typedef double real_num_t;
typedef complex<real_num_t> num_t;
typedef int int_t;

struct real_coordinate_function_args {
  num_t chi;
  real_num_t bh;
  real_num_t z_d;
  real_num_t y_d;
};

typedef num_t (*coordinate_function) (int_t x, int_t z, void* cf_args);
typedef num_t (*real_coordinate_function) (real_num_t x, real_num_t y, real_num_t z, void* cf_args);

int parallel_plate_dislocation(real_coordinate_function_args chi0r_args,
		   real_coordinate_function_args chihr_args,
		   real_coordinate_function_args chihcr_args,
		   real_num_t exposition_x_min,
		   real_num_t exposition_x_max,
		   int_t Nx,
		   real_num_t exposition_y_min,
		   real_num_t exposition_y_max,
		   real_num_t dy,
		   num_t alpha,
		   real_num_t plate_turn,
		   real_num_t plate_depth,
		   real_num_t k0,
		   real_num_t h,
		   real_num_t chi_order,
		   real_num_t precision,
		   num_t** D0_result, 
		   num_t** Dh_result );

int parallel_plate_single_y(num_t chi0,
		   num_t chih,
		   num_t chihc,
		   real_num_t exposition_x_min,
		   real_num_t exposition_x_max,
		   int_t* Nx,
		   int_t* Nz,
		   real_num_t exposition_y,
		   vector3d<real_num_t> dislocation_position,
		   real_num_t alpha,
		   real_num_t plate_turn,
		   real_num_t plate_depth,
		   real_num_t k0,
		   real_num_t h,
		   real_num_t precision,
		   num_t** D0_result, 
		   num_t** Dh_result );

static num_t ii(0.0, 1.0);

inline angle_num_t _bragg_angle(real_num_t h, real_num_t k0) { return asin(h/(2*k0)); };
inline num_t ds_coef (real_num_t k, num_t chi) { return num_t(M_PI*k, 0)*chi/ii; };

inline void mesh_size(const real_num_t& plate_depth, const real_num_t& plate_turn, const real_num_t& precision, const real_num_t& k0, const real_num_t& h, const real_num_t& x_length, const num_t chih, int& Nx, int& Nz) {
  angle_num_t bragg_angle = asin(h/(2*k0));
  real_num_t effective_plate_depth = abs((plate_depth) / cos(plate_turn));
  real_num_t ds = 1/(abs(ds_coef(k0, chih))*precision);
  real_num_t dz =     ds * cos(bragg_angle);
  real_num_t dx = 2 * ds * sin(bragg_angle);
  Nx = x_length/dx + 1;
  Nz = abs(plate_depth)/cos(plate_turn)/dz;
}

void init_tabulations();

#endif
