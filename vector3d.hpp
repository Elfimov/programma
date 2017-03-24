#ifndef _VECTOR3D_HPP_
#define _VECTOR3D_HPP_

#include "matrixnd.hpp"
#include "cmath"

template <typename R> class vector3d{
public:
	vector3d() {};
	vector3d(const R& _x, const R& _y, const R& _z) {x = _x; y = _y; z = _z;};
	R x;
	R y;
	R z;
        R& operator[] (const int& i) { switch (i) { case 0: return x; case 1: return y; case 2: return z; }; };
	const R& operator[] (const int& i) const { switch (i) { case 0: return x; case 1: return y; case 2: return z; }; }; 
};

template <typename R> using matrix3d = vector3d<vector3d<R> >;
template <typename R> using rang3t3d = vector3d<matrix3d<R> >;
template <typename R> using rang4t3d = vector3d<rang3t3d<R> >;

template <typename R> matrix3d<R> transpose(const matrix3d<R>& in) {
	matrix3d<R> out;
	out.x.x = in.x.x;	out.x.y = in.y.x;	out.x.z = in.z.x;
	out.y.x = in.x.y;	out.y.y = in.y.y;	out.y.z = in.z.y;
	out.z.x = in.x.z;	out.z.y = in.y.z;	out.z.z = in.z.z;
	return out;
};

template <typename R> R trace (const matrix3d<R>& in) {
	return in.x.x + in.y.y + in.z.z;
};

template <typename R> R dot (const vector3d<R>& v1, const vector3d<R>& v2) {
	return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
};

template <typename R> vector3d<R> cross (const vector3d<R>& v1, const vector3d<R>& v2) {
	return vector3d<R> (v1.y*v2.z-v2.y*v1.z, v1.z*v2.x-v2.z*v1.x, v1.x*v2.y-v2.x*v1.y); };

template <typename R, typename R2> vector3d<R2> tensor(const vector3d<R>& v1, const R2& v2) {
	return vector3d<R2>(v1.x*v2, v1.y*v2, v1.z*v2); };
template <typename R, typename R2> R2 tensor(const R& v1, const R2& v2) { return v1*v2; };

template <template<typename> class R1, template<typename> class R2, typename R> R1<R2<R> > reduce(const vector3d<R1<R> >& v1, const vector3d<R2<R> >& v2) {
	return tensor(v1.x, v2.x) + tensor(v1.y, v2.y) + tensor(v1.z, v2.z); };
	
template <template<typename> class R1, template<typename> class R2, typename R> R1<R2<R> > reduce2(const matrix3d<R1<R> >& v1, const matrix3d<R2<R> >& v2) {
	return  tensor(v1.x.x, v2.x.x) + tensor(v1.x.y, v2.x.y) + tensor(v1.x.z, v2.x.z) + 
		tensor(v1.y.x, v2.y.x) + tensor(v1.y.y, v2.y.y) + tensor(v1.y.z, v2.y.z) + 
		tensor(v1.z.x, v2.z.x) + tensor(v1.z.y, v2.z.y) + tensor(v1.z.z, v2.z.z); };
		
template <class R, class R2> R2 reduce(const vector3d<R>& v1, const vector3d<R2>& v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z; };
	
template <class R, class R2> R2 reduce2(const matrix3d<R>& v1, const matrix3d<R2>& v2) {
	R2 r = v1.x.x*v2.x.x + v1.x.y*v2.x.y + v1.x.z*v2.x.z + 
	       v1.y.x*v2.y.x + v1.y.y*v2.y.y + v1.y.z*v2.y.z + 
	       v1.z.x*v2.z.x + v1.z.y*v2.z.y + v1.z.z*v2.z.z; 
	return r; };

template <typename R> const rang3t3d<R> levichita3d(R*) {
    return rang3t3d<R>(matrix3d<R>(vector3d<R>(0,0,0),vector3d<R>(0,0,1),vector3d<R>(0,-1,0)),
		       matrix3d<R>(vector3d<R>(0,0,-1),vector3d<R>(0,0,0),vector3d<R>(1,0,0)),
		       matrix3d<R>(vector3d<R>(0,1,0),vector3d<R>(-1,0,0),vector3d<R>(0,0,0))); };

template <typename R> const matrix3d<R> kronecker(R*) {
    return matrix3d<R>(vector3d<R>(1,0,0),vector3d<R>(0,1,0),vector3d<R>(0,0,1)); };
    
template <typename R, typename R2> vector3d<R2> operator*(const R& s, const vector3d<R2>& T) {
    return vector3d<R2>(s*T.x, s*T.y, s*T.z); };
template <typename R, typename R2> vector3d<R2> operator*(const vector3d<R2>& T, const R& s) {
    return vector3d<R2>(s*T.x, s*T.y, s*T.z); };
template <typename R, typename R2> vector3d<R2> operator/(const vector3d<R2>& T, const R& s) {
    return vector3d<R2>(T.x/s, T.y/s, T.z/s); };


template <typename R> vector3d<R> operator* (const vector3d<R>& v, const R& s) { return vector3d<R>(v.x*s, v.y*s, v.z*s); };
template <typename R> vector3d<R> operator* (const R& s, const vector3d<R>& v) { return v*s; };
template <typename R> vector3d<R> operator/ (const vector3d<R>& v, const R& s) { return vector3d<R>(v.x/s, v.y/s, v.z/s); };
template <typename R> vector3d<R> operator+ (const vector3d<R>& v1, const vector3d<R>& v2) { return vector3d<R>(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z); };
template <typename R> vector3d<R> operator- (const vector3d<R>& v1, const vector3d<R>& v2) { return vector3d<R>(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z); };
template <typename R> vector3d<R> operator- (const vector3d<R>& v) {return vector3d<R>(-v.x, -v.y, -v.z);};
template <typename R> vector3d<R> operator+= (vector3d<R>& v1, const vector3d<R>& v2) { return vector3d<R>(v1.x+=v2.x, v1.y+=v2.y, v1.z+=v2.z); };

template <typename R> matrix3d<R> operator* (const matrix3d<R>& m, const R& s) {return matrix3d<R>(m.x*s, m.y*s, m.z*s); };
template <typename R> matrix3d<R> operator/ (const matrix3d<R>& m, const R& s) {return matrix3d<R>(m.x/s, m.y/s, m.z/s); };
template <typename R> vector3d<R> operator* (const matrix3d<R>& m, const vector3d<R>& v) { return vector3d<R>(m.x.x*v.x+m.x.y*v.y+m.x.z*v.z, m.y.x*v.x+m.y.y*v.y+m.y.z*v.z, m.z.x*v.x+m.z.y*v.y+m.z.z*v.z); };

template <typename R> R norm (const vector3d<R>& v) { return sqrt(v.x*v.x+v.y*v.y+v.z*v.z); };
template <typename R> vector3d<R> normalize (const vector3d<R>& v) {return v/norm(v); };
template <typename R> bool isnan (const vector3d<R>& v) {return isnan(v.x) || isnan(v.y) || isnan(v.z); };
template <typename R> matrix3d<R> conj(const matrix3d<R>& i) { return matrix3d<R>(vector3d<R>(i.x.x, i.y.x, i.z.x), vector3d<R>(i.x.y,i.y.y,i.z.y), vector3d<R>(i.x.z,i.y.z,i.z.z)); };

template <typename R> matrixnd<R, 3> tond(const matrix3d<R>& m) {
	matrixnd<R,3> o;
	o.e[0] = m.x.x;	o.e[1] = m.x.y; o.e[2] = m.x.z;
	o.e[3] = m.y.x;	o.e[4] = m.y.y; o.e[5] = m.y.z;
	o.e[6] = m.z.x;	o.e[7] = m.z.y; o.e[8] = m.z.z;
	return o;
}

template <typename R> matrix3d<R> to3d(const matrixnd<R, 3>& m) {
	return matrix3d<R>(vector3d<R>(m.e[0], m.e[1], m.e[2]), vector3d<R>(m.e[3], m.e[4], m.e[5]), vector3d<R>(m.e[6], m.e[7], m.e[8])); };
template <typename R> matrix3d<R> invert(const matrix3d<R>& m) { return to3d(invert((tond(m)))); };


#endif
