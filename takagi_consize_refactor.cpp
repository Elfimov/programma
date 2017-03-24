#include "takagi.h"
#include "tabulate.hpp"
#include "vector3d.hpp"
#include <algorithm>
#include <iostream>
#include <functional>
#include "green_tensor.hpp"
#include "straight_dislocation.hpp"

tabulation2d<matrix3d<real_num_t>, real_num_t>* green_function_tabulation;
tabulation2d<rang3t3d<real_num_t>, real_num_t>* green_function_derivative_tabulation;
rang4t3d<real_num_t> D;

template <typename R> R phi (const vector3d<R>& d) {
  return atan2(d.y, d.x); }

template <typename real_num_t> matrix3d<real_num_t> dislocation_strain(const rang4t3d<real_num_t>& elastic_tensor, const vector3d<real_num_t>& point, 
		std::function<vector3d<real_num_t>(const real_num_t&) > dislocation_line, 
		std::function<rang3t3d<real_num_t>(const vector3d<real_num_t>& ) > green_function_derivative, 
		const vector3d<real_num_t>& burgers) {
    const int integration_length = 1000;
    const real_num_t line_derivative_step = 1.0f/integration_length;
    vector3d<real_num_t> v0( 0, 0, 0);
    matrix3d<real_num_t> m0(v0,v0,v0);
    rang3t3d<real_num_t> dislocation_line_integral(m0,m0,m0);
    rang3t3d<real_num_t> bC = reduce(burgers, elastic_tensor); // indeces: jkl
    
    for (int i=0; i< integration_length; i++) {
	real_num_t t = i*line_derivative_step;
	vector3d<real_num_t> line_derivative = (dislocation_line(t+line_derivative_step) - dislocation_line(t-line_derivative_step))/2.0f;
	rang3t3d<real_num_t> gfd = green_function_derivative(dislocation_line(t)-point);
	matrix3d<real_num_t> bc_gfd = matrix3d<real_num_t>(reduce2(bC.x, gfd), reduce2(bC.y, gfd), reduce2(bC.z, gfd)); // complexity N^4
	dislocation_line_integral += tensor(line_derivative, bc_gfd);
    };
    rang3t3d<real_num_t> epsilon = levichita3d((real_num_t*) NULL);
    return reduce2(epsilon, dislocation_line_integral);
}

class deformation_t {
public:
//	static constexpr real_num_t nu=0.25;
	real_num_t    nu;
	real_num_t    bh;
	real_num_t    bth;
	real_num_t    bph;
	real_num_t    bxh;
	real_num_t    bzh;
	real_num_t    y;
	real_num_t    dx;
	real_num_t    dz;
	real_num_t    dsh;
	real_num_t    K;
	real_num_t    x_min;
	real_num_t    z_min;
	real_num_t    Nz;
	num_t    A;
	num_t    B;
	num_t    C;
	angle_num_t   rotation;
	vector3d<real_num_t> h_dislocation_basis;
	vector3d<real_num_t> dsh_dislocation_basis;
	vector3d<real_num_t> dislocation_position;
	matrix3d<real_num_t> dislocation_basis;
	matrix3d<real_num_t> lattice_basis;
	vector3d<real_num_t> h_crystal_basis;
	vector3d<real_num_t> b_crystal_basis;
	vector3d<real_num_t> t_lattice_basis;
	vector3d<real_num_t> dsh_crystal_basis;
	
	num_t chi_transform(num_t chi, num_t step, num_t K) { return chi * step * K * num_t(M_PI, 0) / ii; };
	deformation_t (vector3d<real_num_t> h, vector3d<real_num_t> n, vector3d<real_num_t> b, vector3d<real_num_t> t, vector3d<real_num_t> dislocation_position, real_num_t y, real_num_t dx, real_num_t dz, real_num_t dsh, real_num_t K, real_num_t x_min, real_num_t z_min, real_num_t Nz, angle_num_t rotation, num_t chih, num_t chihc, real_num_t alpha) {
		this->nu=0.25;
		this->nu=0.0;
		this->y = y;
		this->dx = dx;
		this->dz = dz;
		this->dsh = dsh;
		this->K = K;
		this->x_min = x_min;
		this->z_min = z_min;
		this->Nz = Nz;
		this->rotation = rotation;
		this->dislocation_position = dislocation_position;
		this->bh=dot(b,h);
		this->bth=dot(b-normalize(t)*dot(b,normalize(t)),h);
		this->bph=dot(cross(normalize(t),b),h);
		this->h_crystal_basis = h;
		this->b_crystal_basis = b;
		this->t_lattice_basis = t;
		this->lattice_basis = matrix3d<real_num_t>(normalize(h), cross(normalize(n), normalize(h)), normalize(n));
		this->dsh_crystal_basis = reduce(vector3d<real_num_t>(dx/2, dz*sin(rotation), dz*cos(rotation)), lattice_basis);
		vector3d<real_num_t> ex,ey,ez;

		ez=reduce(t,conj(lattice_basis));
		ez=normalize(ez);
		ex=normalize(reduce(b,conj(lattice_basis)));
		ex=ex-ez*dot(ez,ex);
//		ex=normalize(reduce(b,conj(lattice_basis))-ez*dot(ez, reduce(b,conj(lattice_basis))));
		if (isnan(ex) || norm(ex)<1e-8) ex=normalize(vector3d<real_num_t>(M_PI,M_PI*M_PI,M_PI*M_PI*M_PI)-ez*dot(ez, vector3d<real_num_t>(M_PI,M_PI*M_PI,M_PI*M_PI*M_PI)));
		ex=normalize(ex);
		ey=cross(ez,ex);
		this->bxh=norm(h)*dot(reduce(b,conj(lattice_basis)),ex);
		this->bzh=norm(h)*dot(reduce(b,conj(lattice_basis)),ez);	

		dislocation_basis = matrix3d<real_num_t>(ex,ey,ez);
		this->h_dislocation_basis = conj(dislocation_basis)[0]; /* due to h={h, 0, 0} */
		this->dsh_dislocation_basis = reduce(vector3d<real_num_t>(dx/2, dz*sin(rotation), dz*cos(rotation)),conj(dislocation_basis)); /* dsh = {dx/2, 0, dz} */
		A = chi_transform(chih,  dsh, K);
		B = chi_transform(chihc, dsh, K);
		C = chi_transform(alpha, dsh, K);
	};
	deformation_t() {};
	void indeces_to_beam (const int_t & x_indeces, const int_t & z_indeces, real_num_t & x_beam, real_num_t & z_beam) const {
	  x_beam = x_min - (Nz - (z_indeces + 0.5)) * dx / 2  + x_indeces * dx;
	  z_beam = z_min + (z_indeces + 0.5) * dz;
	}
	void beam_to_crystal (const real_num_t & y_beam, const real_num_t & z_beam, real_num_t & y_crystal, real_num_t & z_crystal) const {
	  y_crystal = y_beam/cos(rotation) + z_beam*sin(rotation);
	  z_crystal = z_beam*cos(rotation);  
	}
	num_t var_alpha_indeces (const int_t & x_indeces, const int_t & z_indeces) const {
/* z-> z+1/2*/
	  real_num_t x_beam, z_beam;
	  real_num_t y_crystal, z_crystal;
	  indeces_to_beam (x_indeces, z_indeces, x_beam, z_beam);
	  beam_to_crystal(y, z_beam, y_crystal, z_crystal);
	  return var_alpha_dislocation(x_beam, y_crystal, z_crystal);
	}
	// calcutaes for deformation projection for non-straight dislocation. Uncomment	lines to calculate non-straight dislocation strain
	num_t var_alpha_dislocation (real_num_t x, real_num_t y, real_num_t z) const {
	    vector3d<real_num_t> position(x,y,z);
	    vector3d<real_num_t> crystal_coordinates = reduce(position,lattice_basis);
	    vector3d<real_num_t> t = this->t_lattice_basis;
	    std::function <vector3d<real_num_t>(const real_num_t&) > straight_dislocation_line = [t] (const real_num_t& p) { return t*(p-0.5)*100; }; // dislocation line function !!!
	    std::function<const matrix3d<real_num_t>(const vector3d<real_num_t>& )> green_tensor_tabulated = [green_function_tabulation] (const vector3d<real_num_t>& d) { return green_function_tabulation->bilinear(d.z, phi(d)); };
	    std::function<const rang3t3d<real_num_t>(const vector3d<real_num_t>& )> green_tensor_derivative_tabulated = [green_function_derivative_tabulation] (const vector3d<real_num_t>& d) { return green_function_derivative_tabulation->bilinear(d.z, phi(d)); };
	    std::function<rang3t3d<real_num_t>(const vector3d<real_num_t>& )> green_tensor_derivative_tabulated_callback = [D, green_tensor_tabulated, green_tensor_derivative_tabulated] (const vector3d<real_num_t>& d) 
		{ return full_green_tensor_derivative(d, green_tensor_tabulated, green_tensor_derivative_tabulated); };
//	    matrix3d<real_num_t> strain = dislocation_strain(D, crystal_coordinates, straight_dislocation_line, green_tensor_derivative_tabulated_callback, b_crystal_basis);
//	    real_num_t alpha = dot(h_crystal_basis, reduce(dsh_crystal_basis, strain))*2.0f*M_PI;
	    num_t alpha_old = var_alpha_dislocation_old(x, y, z); // comment this not to calculate straight dislocation strain
//	    return num_t(0, -alpha);
	    return alpha_old;
	}

	num_t var_alpha_dislocation_old (real_num_t x, real_num_t y, real_num_t z) const {
	  vector3d<real_num_t> dislocation_coordinates = reduce(vector3d<real_num_t>(x,y,z)-dislocation_position,conj(dislocation_basis));
	  num_t v(dislocation_coordinates[0], dislocation_coordinates[1]);
	  real_num_t alpha = alpha_local(real(v), imag(v));
	  return num_t(0, -alpha);
	}
	matrix3d<real_num_t> deformation_derivatives(const real_num_t& x, const real_num_t& y) const {
	  matrix3d<real_num_t> J;
	  real_num_t iR2 = 1/(x*x+y*y);
	  if (iR2==0) return matrix3d<real_num_t>(vector3d<real_num_t>(0,0,0),vector3d<real_num_t>(0,0,0),vector3d<real_num_t>(0,0,0));
	  real_num_t C = bxh/(2*(1-nu))*iR2;
	  real_num_t inu=1/(2*(1-nu));
//	  (du_i/dx_k)ki

	  J[0][0]=-C*y*(1-2*nu+2*x*x*iR2);
	  J[0][1]=-C*x*(2-2*nu-(x*x-y*y)*iR2);
	  J[0][2]=-bzh*y*iR2;
	  J[1][0]= C*x*(3-2*nu-2*y*y*iR2);
	  J[1][1]=-C*y*( -2*nu-(x*x-y*y)*iR2);
	  J[1][2]= bzh*x*iR2;
	  J[2]=vector3d<real_num_t>(0,0,0);

	  matrix3d<real_num_t> I = straight_dislocation_strain(vector3d<real_num_t>(1,0,0), vector3d<real_num_t>(bzh,bxh,0), vector3d<real_num_t>(0,x,y), nu);
	  J[0][0] = I[1][1]; J[0][1] = I[1][2]; J[0][2] = I[1][0];
	  J[1][0] = I[2][1]; J[1][1] = I[2][2]; J[1][2] = I[2][0];
	  J[2][0] = I[0][1]; J[2][1] = I[0][2]; J[2][2] = I[0][0];
	  return J;
	}
	vector3d<real_num_t> hu_derivatives(const real_num_t& x, const real_num_t& y) const {
	  real_num_t inu = 1/(2*(1-nu));
	  real_num_t iR2 = 1/(x*x+y*y);
	  if (iR2==0) return vector3d<real_num_t>(0,0,0);
	  real_num_t dx = bh*(-y)*iR2+bth*y*(y*y-x*x)*inu*iR2*iR2-bph*(x*iR2-x*(x*x-y*y)*inu*iR2*iR2);
	  real_num_t dy = bh*x*iR2+bth*x*(x*x-y*y)*inu*iR2*iR2-bph*((-2*nu)*y*iR2*inu-y*(x*x-y*y)*inu*iR2*iR2);
	  return vector3d<real_num_t>(dx, dy, 0);
	}
	real_num_t alpha_local(const real_num_t& x, const real_num_t& y) const {
//	  return dot(h_dislocation_basis*deformation_derivatives(x,y), dsh_dislocation_basis);
	  matrix3d<real_num_t> J = deformation_derivatives(x,y);
	  vector3d<real_num_t> hu_derivatives1(dot(J[0],h_dislocation_basis),dot(J[1],h_dislocation_basis),0);
	  vector3d<real_num_t> hu_derivatives2 = hu_derivatives(x,y);
	  real_num_t alpha = dot(hu_derivatives2, dsh_dislocation_basis);
	  return alpha;
	}

};

// this function solves finite difference scheme
// with known deformation field (deformation)
int takagi_solve_midpoint_fullout_dimless (int_t Nx, int_t Nz, num_t C, num_t* D0_0, num_t* Dh_0, num_t* D0, num_t* Dh, num_t A, num_t B, const deformation_t & deformation) {
  int_t Nx_z_current;
  num_t *D0_current_line, *Dh_current_line;
  num_t *D0_next_line, *Dh_next_line;
  num_t var_alpha_value;
  num_t rhs1, rhs2, m11, m12, m21, m22, det;
  
  D0_current_line = D0;
  Dh_current_line = Dh;
  D0_next_line = &(D0[Nx+Nz-1]);
  Dh_next_line = &(Dh[Nx+Nz-1]);
  
  // copy the initial values into the first row
  Nx_z_current = Nz + Nx - 1;
  for (int_t x = 0; x < Nx_z_current; x++) {
    D0_current_line[x]=D0_0[x];
    Dh_current_line[x]=Dh_0[x];
  }

  // get solution for the other rows using simpson's rule
  for (int_t z = 1; z < Nz; z++) {
    Nx_z_current = Nx + (Nz - z - 1); 
    for (int_t x = 0; x < Nx_z_current; x++) {
      /*
       * m11*D0(x,z)+m12*Dh(x,z)=rhs1,
       * m21*D0(x,z)+m22*Dh(x,z)=rhs2;
       */

      var_alpha_value = deformation.var_alpha_indeces(x, z); // this is where the deformations enter!!!

      rhs1=D0_current_line[x+1]+						    Dh_current_line[x+1]*B*num_t(0.5, 0.0);
      rhs2=Dh_current_line[x]*(num_t(1.0, 0.0)-num_t(0.5, 0.0)*(C-var_alpha_value))+D0_current_line[x  ]*A*num_t(0.5, 0.0);
      m11 = num_t(1.0,0.0);
      m12 =		   - num_t(0.5,0.0)*B;
      m21 =		   - num_t(0.5,0.0)*A;
      m22 = num_t(1.0,0.0) + num_t(0.5,0.0)*(C-var_alpha_value);      

      det = m11*m22-m12*m21;

      D0_next_line[x] = m22/det*rhs1 - m12/det*rhs2;
      Dh_next_line[x] = m11/det*rhs2 - m21/det*rhs1;
    }
    D0_current_line=&(D0_current_line[Nx+Nz-1]);
    Dh_current_line=&(Dh_current_line[Nx+Nz-1]);
    D0_next_line=&(D0_next_line[Nx+Nz-1]);
    Dh_next_line=&(Dh_next_line[Nx+Nz-1]);
  }
  
  return 0;
}

int parallel_plate_single_y(num_t chi0,
		   num_t chih,
		   num_t chihc,
		   real_num_t exposition_x_min,
		   real_num_t exposition_x_max,
		   int_t* Nx,
		   int_t* Nz,
		   real_num_t y,
		   vector3d<real_num_t> dislocation_position,
		   real_num_t alpha,
		   real_num_t plate_turn,
		   real_num_t plate_depth,
		   real_num_t k0,
		   real_num_t h,
		   real_num_t precision,
		   num_t** D0_result, 
		   num_t** Dh_result )
{
  real_num_t effective_plate_depth;
  real_num_t bragg_angle;
  real_num_t dz, ds, dx;
  real_num_t real_x_max, real_x_min;
  num_t *D0, *Dh, *D0_0, *Dh_0, *D0_hires, *Dh_hires;
  deformation_t args;
  
  bragg_angle = asin(h/(2*k0));
  
  effective_plate_depth = abs((plate_depth) / cos(plate_turn));
  ds = 1/(abs(ds_coef(k0, chih))*precision);
  dz =     ds * cos(bragg_angle);
  dx = 2 * ds * sin(bragg_angle);
  *Nx = (exposition_x_max-exposition_x_min)/dx + 1;
  *Nz = abs(plate_depth)/cos(plate_turn)/dz;
  
  D0_0 = new num_t[(*Nx)+(*Nz)-1];
  Dh_0 = new num_t[(*Nx)+(*Nz)-1];
  
  D0 = new num_t[((*Nx)+(*Nz)-1)*(*Nz)];
  Dh = new num_t[((*Nx)+(*Nz)-1)*(*Nz)];
  
  real_num_t y_shift = y*sin(plate_turn)*tan(bragg_angle);
  for (int_t x = 0; x < (*Nx)+(*Nz)-1; x++) {
    D0_0[x] = exp(num_t(0.0, ((x-(*Nx+*Nz-1)/2)*dx-y_shift*2)*((x-(*Nx+*Nz-1)/2)*dx-y_shift*2)*k0*1e-5)); // spherical waves
//    D0_0[x] = num_t(1.0, 0.0); // plane waves
    Dh_0[x] = num_t(0.0, 0.0);
  }

  vector3d<real_num_t> hv(2,2,0); 			// Change here for diffraction vector orientation
  vector3d<real_num_t> b(0.5,0.5,0);			// Burgers vector
  vector3d<real_num_t> n(1,-1,-1);			// Plane normal (must be perpendicular to hv)
  vector3d<real_num_t> t(1,0,1); 			// Straight dislocation line direction
/*  vector3d<real_num_t> t(1,0,0);
  vector3d<real_num_t> hv(0,0,1);
  vector3d<real_num_t> b(0, 0, -1);
  vector3d<real_num_t> n(-1, 0, 0); */

  real_num_t x_min = exposition_x_min - y_shift;
  args = deformation_t(hv,n,b,t,dislocation_position, y, dx, dz, ds, k0, x_min, -effective_plate_depth/2, *Nz, plate_turn, chih, chihc, alpha);
    
  takagi_solve_midpoint_fullout_dimless(*Nx, *Nz, args.C, D0_0, Dh_0, D0, Dh, args.B, args.A, args);
  for (int_t z=0; z<(*Nz); z++) {
    num_t factor = exp(-M_PI*k0*z*dz/cos(bragg_angle)*chi0*ii);
    for (int_t x=0; x<(*Nx+*Nz-1); x++) {
      D0[z*(*Nx+*Nz-1)+x] = D0[z*(*Nx+*Nz-1)+x]*factor;
      Dh[z*(*Nx+*Nz-1)+x] = Dh[z*(*Nx+*Nz-1)+x]*factor;
    }
  }
  //std::cerr << *Nx << endl << *Nz << endl;
  *D0_result = D0;
  *Dh_result = Dh;
  
  delete[] Dh_0;
  delete[] D0_0;

  return 0;
}

void init_tabulations() {
  real_num_t nu = 0.0;
  real_num_t c12=nu; real_num_t c44 = (1-nu*2.0)/2.0; real_num_t c11 = 1-nu;
  c11 = 165.6; // silicon
  c12 = 63.9;
  c44 = 79.5;
  D = cubic_elastic_tensor(c11,c12,c44);
    
  std::function<const matrix3d<real_num_t>(const real_num_t&, const real_num_t& )> green_tensor_by_angles = [D] (const real_num_t& zcos, const real_num_t& phi) {vector3d<real_num_t> d(sqrt(1-zcos*zcos)*cos(phi),sqrt(1-zcos*zcos)*sin(phi),zcos); return green_tensor(d, D);};
  std::function<const rang3t3d<real_num_t>(const real_num_t&, const real_num_t& )> green_tensor_derivative_by_angles = [D] (const real_num_t& zcos, const real_num_t& phi) { vector3d<real_num_t> d(sqrt(1-zcos*zcos)*cos(phi),sqrt(1-zcos*zcos)*sin(phi),zcos); return green_tensor_derivative(d, D);};

  green_function_tabulation = new tabulation2d<matrix3d<real_num_t>, real_num_t>(green_tensor_by_angles, -1, 1, 10, -M_PI, M_PI, 10);
  green_function_derivative_tabulation = new tabulation2d<rang3t3d<real_num_t>, real_num_t>(green_tensor_derivative_by_angles, -1, 1, 10, -M_PI, M_PI, 10);
}
