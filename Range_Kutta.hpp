// Для использования:
// Нахождение количества дифферинциальных уравнений в системе
// В исполняемом файле завести Вектор переменных, для которых нужно считать производную Vec<D>, где D - их количество и проинициализировать их, соблюдая определенный самим тобой порядок
// Реализовать функции, которые получают на вход вектор параметров системы и считают производные для каждого параметра. Пример: double dx_dt(Vec<D>){return...};
// Сделать вектор ссылок на эти функции, соблюдая порядок какой был установлен для переменных
// Теперь можно использовать функцию range_kutta(вектор пересенных, вектор указателей на функции, шаг);
// Он возвращает изменненый вектор парметров, сделав один шаг

#include "vec.hpp"

const double sixth = 1.0 / 6.0;

template <int D, typename T1>
Vec<D, T1> range_kutta(Vec<D,T1> &a, Vec<D, double (*)(Vec<D,T1>)> const &fx_i, double dt){
	Vec<D,T1> k1;
	Vec<D,T1> half_k1;
	Vec<D,T1> k2;
	Vec<D,T1> half_k2;
	Vec<D,T1> k3; 
	Vec<D,T1> k4;
	for(int i = 0; i < D; i++){
		k1[i] = fx_i[i](a) * dt;
		half_k1[i] = k1[i] * 0.5;

	}

	Vec <D,T1> tmp = a + half_k1; 

	for(int i = 0; i < D; i++){
		k2[i] = fx_i[i](tmp) * dt;
		half_k2[i] = k2[i] * 0.5;
	}

	tmp = a + half_k2;

	for(int i = 0; i < D; i++){
		k3[i] = fx_i[i](tmp) * dt;
	}

	tmp = a + k3;

	for(int i = 0; i < D; i++){
		k4[i] = fx_i[i](tmp) * dt;
	}

	a = a + sixth * (k1 + 2*(k2 + k3) + k4);

	return a;
}