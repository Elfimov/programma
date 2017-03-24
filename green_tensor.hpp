#include "vector3d.hpp"
#include <cmath>
#include <functional>

template <typename R> matrix3d<R> green_tensor(const vector3d<R>& input_direction, const rang4t3d <R>& elastic_tensor) {
    const int N = 200;
    vector3d<R> r1a = vector3d<R> (1,0,0);
    vector3d<R> r2a = vector3d<R> (0,1,0);
    vector3d<R> rotation_axis;
    matrix3d<R> sum = matrix3d<R>(vector3d<R>(0,0,0), vector3d<R>(0,0,0), vector3d<R>(0,0,0));
    vector3d<R> direction = normalize(input_direction);
    
    if (abs(dot(r1a,direction))<abs(dot(r2a,direction)))
	r1a = normalize(r1a-direction*dot(r1a,direction));
    else
	r1a = normalize(r2a-direction*dot(r2a,direction));
    
    r2a = cross(r1a, direction);
    
    for (int i = 0;  i < N; i++) {
	R phi = (2.0f*M_PI*i)/N;
	vector3d<R> idirection = r1a*cos(phi) + r2a*sin(phi);
	matrix3d<R> reduced_elastic_tensor = reduce(idirection, transpose(reduce(idirection, elastic_tensor)));
	sum += invert(reduced_elastic_tensor);
    };
    return sum*(1.0f/N);
};

template <typename R> rang3t3d<R> green_tensor_derivative(const vector3d<R>& direction, const rang4t3d <R>& elastic_tensor) {
    R h = 0.00001;
    vector3d<R> dx = vector3d<R> (h,0,0);
    vector3d<R> dy = vector3d<R> (0,h,0);
    vector3d<R> dz = vector3d<R> (0,0,h);
    
    rang3t3d<R> out = rang3t3d<R> ((green_tensor(direction+dx, elastic_tensor)-green_tensor(direction-dx, elastic_tensor))/(2*h),
				   (green_tensor(direction+dy, elastic_tensor)-green_tensor(direction-dy, elastic_tensor))/(2*h),
				   (green_tensor(direction+dz, elastic_tensor)-green_tensor(direction-dz, elastic_tensor))/(2*h));
    return out;
};

template <typename R> rang3t3d<R> full_green_tensor_derivative(const vector3d<R>& position, 
		std::function<const matrix3d<R>(const vector3d<R>&)> green_tensor_by_direction, 
		std::function<const rang3t3d<R>(const vector3d<R>&)> green_tensor_derivative_by_direction) {
    R distance = norm(position);
    rang3t3d<R> result1 = -tensor(position,green_tensor_by_direction(normalize(position)))/(distance*distance*distance);
    rang3t3d<R> result2 = green_tensor_derivative_by_direction(normalize(position))/(distance*distance);
    return (result1+result2)/(4.0f*M_PI);
};

template <typename R> matrix3d<R> isotropic_green_tensor(const vector3d<R>& input_direction, const R& c12, const R& c44) {
    return (c12+3*c44)/(2*c44*(c12+2*c44))*kronecker((R*) NULL)+(c12+c44)/(2*c44*(c12+2*c44))*tensor(input_direction, input_direction); }
    
template <typename R> rang4t3d<R> cubic_elastic_tensor(const R& c11, const R& c12, const R& c44) {
    vector3d<R> v0( 0, 0, 0);
    matrix3d<R> m0(v0,v0,v0);
    rang3t3d<R> r0(m0,m0,m0);
    rang4t3d<R> C( r0,r0,r0);
    C.x.x.x.x=C.y.y.y.y=C.z.z.z.z=c11;
    C.x.x.y.y=C.x.x.z.z=C.y.y.x.x=C.y.y.z.z=C.z.z.x.x=C.z.z.y.y=c12;
    C.x.y.x.y=C.x.z.x.z=C.y.x.y.x=C.y.z.y.z=C.z.x.z.x=C.z.y.z.y=
    C.x.y.y.x=C.x.z.z.x=C.y.x.x.y=C.y.z.z.y=C.z.x.x.z=C.z.y.y.z=c44;
    return C;
};