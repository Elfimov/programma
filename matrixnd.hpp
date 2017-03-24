#ifndef _MATRIXND_HPP_
#define _MATRIXND_HPP_

#include <cmath>

template <typename R, int n> class matrixnd {
public:
  R e[n*n];
};

template <typename R, int n> matrixnd<R, n> operator*(const matrixnd<R,n>& m1, const matrixnd<R,n>& m2) { 
  matrixnd<R, n> r;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      r.e[i*n+j] = 0;
      for (int k=0; k<n; k++) {
        r.e[i*n+j] += m1.e[i*n+k]*m2.e[k*n+j];
      }
    }
  };
  return r;
};

template <typename R> R determinant(const matrixnd<R, 1>& m) { return m.e[0]; };
template <typename R> R determinant(const matrixnd<R, 2>& m) { return m.e[0]*m.e[3]-m.e[1]*m.e[2]; };
 

template <typename R, int n> R determinant(const matrixnd<R, n>& m) {
  R d = 0;
  for (int i=0; i<n; i++) d += m.e[i]*algebraic_extension(m, 0, i);
  return d;
}

template <typename R, int n> R algebraic_extension(const matrixnd<R, n>& m, int e_i, int e_j) {
  matrixnd<R, n-1> r;
  for (int i=0, r_i = 0; i<n; i++) {
    if (i==e_i) continue;
    for (int j=0, r_j = 0; j<n; j++) {
      if (j==e_j) continue;
      r.e[r_i*(n-1)+r_j]=m.e[i*n+j];
    r_j++;
    }
  r_i++;
  }
  
  R c = -((e_i+e_j) % 2)*2+1;
  return c* determinant(r);
};


template <typename R, int n> matrixnd<R, n> invert(const matrixnd<R, n>& m) {
  matrixnd<R, n> I;
  R d = determinant(m);
  for (int i=0; i<n; i++) {
  for (int j=0; j<n; j++) {
    I.e[i+n*j] = algebraic_extension(m,i,j)/d;
  }
  }
  return I;
};

/*
template <typename R, int n> matrixnd<R, n-1>& invert(const matrixnd<R, n>& m) {
  matrixnd<R, n> I;
  // initialize inverse matrix with unity
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      I.e[i,j]=0;
    }
    I.e[i,i]=1;
  }
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      for (int k=0; k<n; k++) {
        R multiplier = 1/m.e[k,k];
      }
    }
  }
}*/


#endif