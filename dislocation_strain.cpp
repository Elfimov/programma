#include <functional>
#include "vector3d.hpp"
#include "green_tensor.hpp";
#include "straight_dislocation.hpp";
#include "tabulate.hpp";

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

template <typename real_num_t> real_num_t alpha_var_dislocation(real_num_t x, real_num_t y, real_num_t z) {
    auto straight_dislocation_line = [] (const real_num_t& p) { return vector3d<real_num_t>((p-0.5)*100,0,0); };
    real_num_t nu = 0.25;
    real_num_t c12=nu; real_num_t c44 = (1-nu*2.0)/2.0; real_num_t c11 = 1-nu;
    rang4t3d<real_num_t> C = cubic_elastic_tensor(c11,c12,c44); // isotropic elastic tensor
    
    std::function<const matrix3d<real_num_t>(const real_num_t&, const real_num_t& )> green_tensor_by_angles = [C] (const real_num_t& zcos, const real_num_t& phi) {vector3d<real_num_t> d(sqrt(1-zcos*zcos)*cos(phi),sqrt(1-zcos*zcos)*sin(phi),zcos); return green_tensor(d, C);};
    std::function<const rang3t3d<real_num_t>(const real_num_t&, const real_num_t& )> green_tensor_derivative_by_angles = [C] (const real_num_t& zcos, const real_num_t& phi) { vector3d<real_num_t> d(sqrt(1-zcos*zcos)*cos(phi),sqrt(1-zcos*zcos)*sin(phi),zcos); return green_tensor_derivative(d, C);};
    
    tabulation2d<matrix3d<real_num_t>, real_num_t> green_function_tabulation(green_tensor_by_angles, -1, 1, 400, -M_PI, M_PI, 400);
    tabulation2d<rang3t3d<real_num_t>, real_num_t> green_function_derivative_tabulation(green_tensor_derivative_by_angles, -1, 1, 400, -M_PI, M_PI, 400);
    
    std::function<const matrix3d<real_num_t>(const vector3d<real_num_t>& )> green_tensor_tabulated = [green_function_tabulation] (const vector3d<real_num_t>& d) { return green_function_tabulation.bilinear(d.z, phi(d)); };
    std::function<const rang3t3d<real_num_t>(const vector3d<real_num_t>& )> green_tensor_derivative_tabulated = [green_function_derivative_tabulation] (const vector3d<real_num_t>& d) { return green_function_derivative_tabulation.bilinear(d.z, phi(d)); };
    auto green_tensor_derivative_tabulated_callback = [C, green_tensor_tabulated, green_tensor_derivative_tabulated] (const vector3d<real_num_t>& d) 
	{ return full_green_tensor_derivative(d, green_tensor_tabulated, green_tensor_derivative_tabulated); };
    matrix3d<real_num_t> strain = dislocation_strain(C, vector3d<real_num_t>(x,y,z), straight_dislocation_line, green_tensor_derivative_tabulated_callback, b_crystal_basis);
    return reduce(h_crystal_basis, reduce(dsh, strain));
}