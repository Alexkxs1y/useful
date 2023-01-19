#include "mesh.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

int main(){
	const int Dim = 2;
	int A = 100, B = 100;
	double hx = 1, hy = 1;
	double x0 = 0, y0 = 0;
	double A0 = 1;
	double tau = 1e-1, t = 0;
	double Omega = 2 * M_PI / 200 / tau; 
	double c = 1;
	double k = c * c * tau * tau;
	double quadr_step_x = hx*hx, quadr_step_y = hy*hy;

	Ind<Dim> box = Ind<Dim>(A, B); 
	Vec<Dim> coord = Vec<Dim>(x0, y0);
	Vec<Dim> step = Vec<Dim>(hx, hy);
	Ind<Dim> source = Ind<Dim>(A/4, B/4);

	Mesh<double, Dim> m1 = Mesh<double, Dim>(box, coord, step, 2);
	Mesh<double, Dim> m2 = Mesh<double, Dim>(box, coord, step, 2);
	for(size_t i = 0; i < m1.prod(); i++){
		Ind<Dim> pos = i%box;
		m1[pos] = 0;
		m2[pos] = 0;
	}

	Ind<Dim> x_pos = Ind<Dim>(1,0);
	Ind<Dim> y_pos = Ind<Dim>(0,1);
	
	for (int i = 0; i < 1000; i ++){
		for(size_t j = 0; j < m1.prod(); j++){
			Ind<2> pos = j % box;
			m2[pos] = - m2[pos] + 2*m1[pos] + 
			k * ((m1[pos - x_pos] + m1[pos + x_pos] - 2 * m1[pos])/quadr_step_x + (m1[pos - y_pos] + m1[pos + y_pos] - 2 * m1[pos])/quadr_step_y);
			if(pos == source){
				m2[pos] = A0*cos(Omega * t);
			}
		}
		t += tau;
		m1.swap_data(m2);
		if(i % 10 == 0){
			std::stringstream ss;
			ss << std::setw(4) << std::setfill('0') << i;
			std::string sss = ss.str();
			std::string s = "res/" + sss + ".dat";
			const char* path = s.c_str();
			m1.out2dat(path);
		}
	}
	return 1;
}
