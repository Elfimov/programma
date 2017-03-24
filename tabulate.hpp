#include <functional>

template <typename d, typename R> class tabulation2d {
private:
    d* data;
    R x_min;
    R x_max;
    int Nx;
    R y_min;
    R y_max;
    int Ny;
public:
    tabulation2d(const std::function<d(const R&, const R&)>& function, const R& _x_min, const R& _x_max, const int& _Nx, const R& _y_min, const R& _y_max, const int& _Ny) {
	x_min = _x_min;
	x_max = _x_max;
	Nx = _Nx;
	y_min = _y_min;
	y_max = _y_max;
	Ny = _Ny;
	
	data = new d[Nx*Ny];
	for (int i=0; i<Nx; i++) {
	    R x = x_min + i*(x_max-x_min)/(Nx-1);
	    for (int j=0; j<Ny; j++) {
		R y = y_min + j*(y_max-y_min)/(Ny-1);
		data[j+i*Ny] = function(x, y);
	    }
	}
    }
    d bilinear(const R& x, const R& y) const {
	R ip = (x-x_min)/(x_max-x_min)*(Nx-1);
	R jp = (y-y_min)/(y_max-y_min)*(Ny-1);
	int i = floor(ip); R dip = ip-i;
	int j = floor(jp); R djp = jp-j;
	return data[j+i*Ny]*(1-dip)*(1-djp) + data[(j+1)+i*Ny]*dip*(1-djp) +
	       data[j+(i+1)*Ny]*(1-dip)*djp + data[(j+1)+(i+1)*Ny]*dip*djp;
    }
    ~tabulation2d() {};
    void destroy() { delete[] data; };
};