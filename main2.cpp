#include "takagi.h"
#include "vector3d.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstring>

template <class I, class O> void output_complex_float(const std::string & filename_r, const std::string & filename_i, complex<I>* data, int data_size, O* output_type) {
  O *r = new O[data_size];
  O *i = new O[data_size];
  std::ofstream r_handle, i_handle;

  for (int it = 0; it < data_size; it++) {
    r[it] = real(data[it]);
    i[it] = imag(data[it]);
  }

  r_handle.open(filename_r);
  r_handle.write((char*)r, sizeof(O)*data_size);
  r_handle.close();

  i_handle.open(filename_i);
  i_handle.write((char*)i, sizeof(O)*data_size);
  i_handle.close();



  delete[] r;
  delete[] i;
}

template <class I, class O> void output_real_float(const std::string & filename_r, I* data, int data_size, O* output_type) {
  O *r = new O[data_size];
  std::ofstream r_handle;

  for (int it = 0; it < data_size; it++)
    r[it] = real(data[it]);

  r_handle.open(filename_r);
  r_handle.write((char*)r, sizeof(O)*data_size);
  r_handle.close();

  delete[] r;
}

int main(int argc, char* argv[]) {
  num_t chi0(-0.317e-5, -0.16060E-07);
  num_t chih(0.192e-5, -0.15497E-07);

  real_num_t K = 1e5/0.7093; // wave number in inverse 10 micron units. Simply 1e5/ wavelength in Angstrom.
  real_num_t h = 1e5/1.9201; // 1e5/ interplanar distance of reflection 
  real_num_t alpha_max = 2e-4*h/K; // alpha interval maximum 

//  real_num_t bh = 1;
  
  init_tabulations();
  const int_t max_k_id = 1000; // alpha number (1000 x 2 +1 = 2001 point in alpha interval)
  for (int_t phi = -30; phi <= 30; phi += 30) { // phi points
    int_t Nx, Nz;
    mesh_size(100, phi/(180/M_PI), 20, K, h, 80, chih, Nx, Nz);
// first argument: crystal depth in 10 micron units. 100 units -> 100 units x 10 mkm/units = 1000 microns = 1 mm
// second argument better not touch it (crystal plate turn in radians, if you need to change that better change the for loop)
// third argument: precision. The more this argument, the less the step in the finite difference scheme.
// K: wave number
// h
// sixth argument: topograph x width in 10 micron units
// seventh argument: vospriimchivost ot Stepanova.


    for (int_t y_id = -20; y_id <= 20; y_id += 1) {
      std::ostringstream filename_transmitted_int;
      std::ostringstream filename_diffracted_int;

      vector3d<real_num_t> dislocation_position(0,0,0);
      real_num_t *I0_alpha_sum, *Ih_alpha_sum;


      I0_alpha_sum = new real_num_t[(Nx+Nz-1)*Nz];
      Ih_alpha_sum = new real_num_t[(Nx+Nz-1)*Nz];
      memset(I0_alpha_sum, 0, sizeof(real_num_t)*(Nx+Nz-1)*Nz);
      memset(Ih_alpha_sum, 0, sizeof(real_num_t)*(Nx+Nz-1)*Nz);

      #pragma omp parallel for
      for (int_t alpha_id = -abs(max_k_id); alpha_id < abs(max_k_id)+1; alpha_id++) {
  	num_t *D0, *Dh;

	real_num_t alpha = (alpha_max*alpha_id)/max_k_id;

        parallel_plate_single_y(chi0, chih, chih, -40, 40, &Nx, &Nz, y_id, dislocation_position, alpha,
		  phi/(180/M_PI), 100, K, h, 20, &D0, &Dh);
	//  1-3: vospriimchivost ot stepanova
	// 4: x minimum
	// 5: x maximum. (5)-(4) must give sixth argument of mesh_size call.
	// 6, 7: Nx, Nz
	// 8: y value.
	// 9: don't touch it, probably it's useless.
	// 10: alpha in TT equations
	// 11: plate turn in radians (the same as in mesh_size call)
	// 12: plate depth in 10 micron units (the same as in mesh_size)
	// 13, 14: wave number, diffraction vector length
	// 16: precision (the same as in mesh_size)
	// 17, 18: return values

	if (alpha_id == 0) {
	  std::ostringstream filename_transmitted_real;
	  std::ostringstream filename_transmitted_imag;
	  std::ostringstream filename_diffracted_real;
	  std::ostringstream filename_diffracted_imag;

	  filename_transmitted_real<<"D0r_y"<<std::setfill('0')<<std::setw(4)<<y_id<<"_a"<<std::setfill('0')<<std::setw(4)<<alpha_id<<"_p"<<std::setfill('0')<<std::setw(2)<<phi<<".data";
	  filename_transmitted_imag<<"D0i_y"<<std::setfill('0')<<std::setw(4)<<y_id<<"_a"<<std::setfill('0')<<std::setw(4)<<alpha_id<<"_p"<<std::setfill('0')<<std::setw(2)<<phi<<".data";
	  filename_diffracted_real<<"Dhr_y"<<std::setfill('0')<<std::setw(4)<<y_id<<"_a"<<std::setfill('0')<<std::setw(4)<<alpha_id<<"_p"<<std::setfill('0')<<std::setw(2)<<phi<<".data";
	  filename_diffracted_imag<<"Dhi_y"<<std::setfill('0')<<std::setw(4)<<y_id<<"_a"<<std::setfill('0')<<std::setw(4)<<alpha_id<<"_p"<<std::setfill('0')<<std::setw(2)<<phi<<".data";

	  output_complex_float(filename_transmitted_real.str(), filename_transmitted_imag.str(), D0, Nz*(Nx+Nz-1), (float*) NULL);
	  output_complex_float(filename_diffracted_real.str(),  filename_diffracted_imag.str(),  Dh, Nz*(Nx+Nz-1), (float*) NULL);

	}

	#pragma omp critical
	for (int p=0; p<(Nx+Nz-1)*Nz; p++) {
            real_num_t I0, Ih;
	    I0=real(D0[p])*real(D0[p])+imag(D0[p])*imag(D0[p]);
	    Ih=real(Dh[p])*real(Dh[p])+imag(Dh[p])*imag(Dh[p]);
//	    #pragma omp atomic
	    I0_alpha_sum[p]+=I0;
//	    #pragma omp atomic
	    Ih_alpha_sum[p]+=Ih;
	}

        delete[] D0;
        delete[] Dh;
	std::cerr<<"-";
      }

      filename_transmitted_int<<"I0_y"<<std::setfill('0')<<std::setw(4)<<y_id<<"_p"<<std::setfill('0')<<std::setw(2)<<phi<<".data";
      filename_diffracted_int<<"Ih_y"<<std::setfill('0')<<std::setw(4)<<y_id<<"_p"<<std::setfill('0')<<std::setw(2)<<phi<<".data";

      output_real_float(filename_transmitted_int.str(), I0_alpha_sum, Nz*(Nx+Nz-1), (float*) NULL);
      output_real_float(filename_diffracted_int.str(),  Ih_alpha_sum, Nz*(Nx+Nz-1), (float*) NULL);
      delete[] I0_alpha_sum;
      delete[] Ih_alpha_sum;
      std::cerr<<endl;
      std::cerr << "Plate turn: " << phi << "; Dimensions: Nx=" << Nx << " Nz="<< Nz << std::endl;
    }
  }
  return 0;
}

