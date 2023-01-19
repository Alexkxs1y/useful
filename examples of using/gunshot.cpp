#include "Range_Kutta.hpp"
#include "4401-81.hpp"
#include <math.h>

GHOST4401 gost;
double h = 130, m = 0.000679, V0 = 69.27, d = 0.0045, Cxa = 0.45, u = 0.0107;
double ro = gost.get_density(h);
double g = gost.get_g(h);

// a = (x,y,V,tetta)
double Xa(Vec<4> a){ return Cxa * ro*a[2]*a[2]/2 * M_PI*d*d/4;}

double dx_dt(Vec<4> a){ return a[2]*cos(a[3]);}

double dy_dt(Vec<4> a){ return a[2]*sin(a[3]);}

double dV_dt(Vec<4> a){ return -Xa(a)/m - g*sin(a[3]);}

double dtetta_dt(Vec<4> a){ return -g*cos(a[3])/a[2];}

double interpol(double x, double x1, double y, double y1){return x + (x1 - x) * y/(y - y1);}

//type == 1 - вывод всех параметров, type == 0 - отсутсвие вывода
Vec<2> gunshot (double dt, double x_0, double y_0, double V_0, double tetta_0, int type, int print_iter = 1){     
	double t = 0;
	double y_max = 0;
	Vec<4> parametrs = Vec<4>(x_0, y_0, V_0, tetta_0);
	Vec<4, double (*)(Vec<4>)> diff = Vec<4, double (*)(Vec<4>)>(&dx_dt, &dy_dt, &dV_dt, &dtetta_dt);
	Vec<4> tmp;
	int i = 0;
	while (true){
		if((type == 1) && (i%print_iter == 0))
			std::cout << t << ' ' << parametrs << '\n';
		tmp = parametrs;
		range_kutta(parametrs, diff, dt);
		t+=dt;
		i++;
		if (parametrs[1] > y_max)
			y_max = parametrs[1];
		if (parametrs[1]< 0){
			double y1 = parametrs[1];
			for (int i = 0; i < 4; i++){
				parametrs[i] = interpol(tmp[i], parametrs[i], tmp[1], y1);
			}
			break;
		}
	}
	if(type == 1)
		std::cout << t << ' ' << parametrs << '\n';
	Vec<2> max_xy = Vec<2>(parametrs[0], y_max);
	return max_xy;
}

int main(int argc, char** argv){
	std::cout.precision(4);
	double dt = atof(argv[1]);
	double d_tetta = (M_PI/2) / 90 / 20 ; // 5 сотых градуса в радианах
	double y_pr = 1, delta = 0.01, tetta_pr = 0, x_pr = 0, tetta_max = 0, x_max = 0, y_max = 0;;
	//поиск угла прямого встрела
	for(int i = 0; i < 90*20; i++){
		Vec<2> xy = gunshot(dt, 0, 0, V0, i*d_tetta, 0);
		if((-delta < (xy[1] - y_pr)) &&  ((xy[1] - y_pr) < delta)) {
			tetta_pr = i * d_tetta;
			y_pr = xy[1];
			x_pr = xy[0];
			break; 
		}
	}
	//поиск оптимального угла выстрела
	for(int i = 30*20; i < 50*20; i++){
		Vec<2> xy = gunshot(dt, 0, 0, V0, i*d_tetta, 0);
		if(xy[0] > x_max){
			x_max = xy[0];
			y_max = xy[1];
			tetta_max = i * d_tetta;
		}
	}

	gunshot(dt, 0, 0, V0, u, 1, 10);
	std::cout << '\n' << '\n' << '\n';
	gunshot(dt, 0, 0, V0, tetta_pr, 1, 30);
	std::cout << '\n' << '\n' << '\n';
	gunshot(dt, 0, 0, V0, tetta_max, 1, 100);
	std::cout << '\n' << '\n' << '\n';
	std::cout << "tetta_pr =  " << tetta_pr * 180 / M_PI << "; " << "x_pr = " << x_pr << "; y_pr = "<< y_pr << '\n';
	std::cout << "tetta_max =  " << tetta_max * 180 / M_PI << "; x_max = " << x_max << "; y_max = "<< y_max << '\n';
	std::cout << ro << ' ' << g << '\n';
	return 1;
}