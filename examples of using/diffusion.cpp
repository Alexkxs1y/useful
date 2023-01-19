#include <omp.h>
#include "mesh.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

int main(){
	const int Dim = 2;
	double hx = 0.1, hy = 0.1;
	double x0 = 1, y0 = 1;
	int A = 100, B = 100;
	double quadr_stepx = hx*hx , quadr_stepy = hy*hy; 
	Ind<Dim> box = Ind<Dim>(A, B); 
	Vec<Dim> coord = Vec<Dim>(x0, y0);
	Vec<Dim> step = Vec<Dim>(hx, hy);  
	Mesh<double, Dim> m1 = Mesh<double, Dim>(box, coord, step);
	Mesh<double, Dim> m2 = Mesh<double, Dim>(box, coord, step);
	for(size_t i = 0; i < m1.prod(); i++){
		Ind<Dim> pos = i%box;
		m1[pos] = 0;
	}
	Ind<Dim> pos = Ind<Dim>(A/2, B/2);
	m1[pos] = 10;

	double D = 1;
	double T = 1e-3;
	double coef = D*T;
	Ind<Dim> x_pos = Ind<Dim>(1,0);
	Ind<Dim> y_pos = Ind<Dim>(0,1);
	//m1.out2dat("res/0.dat");
	size_t sz = m1.prod() - A;
	
	double t0 = omp_get_wtime();

	for (int i = 0; i < 1000; i ++){
//# pragma omp parallel for 
		for(size_t j = 0; j < sz; j++){
			Ind<2> pos = j % box;
			m2[pos] = m1[pos] + coef*
			((m1[pos - x_pos] + m1[pos + x_pos] - 2 * m1[pos])/quadr_stepx +
			(m1[pos - y_pos] + m1[pos + y_pos] - 2 * m1[pos])/quadr_stepy);
			/*m2[j] = m1[j] + coef*
			((m1[j-1] + m1[j+1] - 2*m1[j])/quadr_stepx +
			 (m1[j-A] + m1[j+A] - 2*m1[j])/quadr_stepy); */
		}
		/*if(i % 10 == 0){
			std::stringstream ss;
			ss << std::setw(4) << std::setfill('0') << i;
			std::string sss = ss.str();
			//std::string tmp = std::to_string(i);
			std::string s = "res/" + sss + ".dat";
			const char* path = s.c_str();
			m1.out2dat(path);
		}
		m1.swap_data(m2);c*/
	}
	m1.out2dat("res/0.dat");

	std::cout << omp_get_wtime() - t0 << '\n';
	return 1;
}
