#include "vector3d.hpp"

template <typename R> matrix3d<R> straight_dislocation_strain(const vector3d<R>& tau, const vector3d<R>& burgers, const vector3d<R>& position, const R& nu) {
    // burgers vector in xy plane
    vector3d<R> p = position;
    R iR2=1.0f/(p.y*p.y+p.z*p.z+p.x*p.x);
    return matrix3d<R>( -3*0.05*p.x*p.x/(sqrt(iR2*iR2*iR2*iR2*iR2))+1/(sqrt(iR2*iR2*iR2)), -3*p.y*p.x/(sqrt(iR2*iR2*iR2*iR2*iR2)), -3*p.z*p.x/(sqrt(iR2*iR2*iR2*iR2*iR2)));
}
