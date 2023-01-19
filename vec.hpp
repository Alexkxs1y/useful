#ifndef VEC_H
#define VEC_H
#include <iostream>
#include <math.h>

template <int D, typename T = double> class Vec{
	private:
	T p[D];
        void set_x(){}
        template <typename T2, typename ... Args> void set_x(const T2& x, const Args&... xxx){ p[D-1-sizeof...(Args)] = x; set_x(xxx...); }
	public:
    	T operator [](int i) const { return p[i]; }
    	T& operator [](int i) { return p[i]; }
    	Vec(T x=0){ for(int i=0; i<D; i++) p[i] = x; }
        template <typename ... Args> explicit Vec(const Args&... xxx){
                static_assert(sizeof...(Args)==D, "illegal parametrs count!"); 
                set_x(xxx...);
        }
	T abs () const{
		T quadr = 0;
		for (int i = 0 ; i < D; i++ )
			quadr += p[i]*p[i];
		return std::sqrt(quadr);
	}
};

template <int D, typename T1, typename T2> 
Vec<D, decltype(T1()*T2())> operator + (const Vec<D, T1> &a, const Vec<D, T2> &b){
	Vec<D, decltype(T1()*T2())> res;
	for (int i = 0; i < D; i++){
		res[i] = a[i] + b[i];
	}
	return res;
} 

template <int D, typename T1, typename T2> 
Vec<D, decltype(T1()*T2())> operator * (const Vec<D, T1> &a , T2 b){
	Vec<D, decltype(T1()*T2())> res;
	for (int i = 0; i < D; i++){
		res[i] = a[i] * b;
	}
	return res;
}

template <int D, typename T1, typename T2> 
decltype(T1()*T2()) operator * (const Vec<D, T1> &a, const Vec<D, T2> &b){
	decltype(T1()*T2()) res = 0;
	for (int i = 0; i < D; i++){
		res += a[i] * b[i];
	}
	return res;
}

template <typename T1, typename T2>
Vec<3, decltype(T1()*T2())> operator % (const Vec<3, T1> &a, const Vec<3,T2> &b){
	return Vec<3>(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2] , a[0]*b[1] - a[1]*b[0]);
}

template <int D, typename T1, typename T2> 
Vec<D, decltype(T1()*T2())> operator * (T2 b, const Vec<D, T1> &a ){
	Vec<D, decltype(T1()*T2())> res;
	for (int i = 0; i < D; i++){
		res[i] = a[i] * b;
	}
	return res;
}

template <int D, typename T1, typename T2> 
Vec<D, decltype(T1()*T2())> operator - (const Vec<D, T1> &a, const Vec<D, T2> &b){
	Vec<D, decltype(T1()*T2())> res;
	for (int i = 0; i < D; i++){
		res[i] = a[i] - b[i];
	}
	return res;
}

template <int D, typename T1>
Vec<D, T1> operator - (const Vec<D, T1> &a){
	Vec<D, T1> res;
	for (int i = 0; i < D; i++)
		res[i] = -a[i];
	return res;
}

template <int D, typename T1>
std::ostream & operator << (std::ostream &S, const Vec<D, T1> &a){
	for (int i = 0; i < D; i++){
		S << a[i] << ' ';
	}
	return S;
}

template <int D, typename T1, typename T2> 
bool operator < (const Vec<D, T1> &a, const Vec<D, T2> &b){
	for (int i = 0; i < D; i++){
		if(a[i] >= b[i])
			return false; 
	}
	return true;
}

template <int D, typename T1, typename T2> 
bool operator > (const Vec<D, T1> &a, const Vec<D, T2> &b){
	for (int i = 0; i < D; i++){
		if(a[i] <= b[i])
			return false; 
	}
	return true;
}

template <int D, typename T1, typename T2> 
bool operator == (const Vec<D, T1> &a, const Vec<D, T2> &b){
	for (int i = 0; i < D; i++){
		if(a[i] != b[i])
			return false; 
	}
	return true;
}

template <int D, typename T1, typename T2> 
bool operator != (const Vec<D, T1> &a, const Vec<D, T2> &b){
	for (int i = 0; i < D; i++){
		if(a[i] != b[i])
			return true; 
	}
	return false;
}


template<int D, typename T1, typename T2>
Vec<D, decltype(T1()*T2())> operator &(const Vec<D,T1> &a, const Vec<D,T2> &b){
	Vec<D, decltype(T1()*T2())> res;
	for(int i = 0; i < D; i++)
		res[i] = a[i]*b[i];
	return res;
}

template<int D> using Ind = Vec<D, int>;  // Ind<3>

template<int D>
Ind<D> operator % (size_t x, const Ind<D> &p){
	Ind<D> r;
    for(int i=0; i<D; i++){ r[i] = x%p[i]; x /= p[i]; }
    return r;
}

template<int D> class DRange{
	Ind<D> a,b;
public:
	DRange(Ind<D> a_, Ind<D> b_){
		a = a_;
		b = b_;
	}
	struct iterator{
		Ind<D> pos;
		const DRange* range = nullptr;
		const Ind<D>& operator *() const { return pos; }
		bool operator != (const iterator& ){return pos!=range->b;}
		Ind<D> operator ++ (){
			if(pos == (range->b - Ind<D>(1))){
				for(int i = 0; i < D; i++){
					pos[i] = range->b[i];
				}
				return pos;
			}
			for(int i = 0; i < D; i++){
				if(pos[i] + 1 < range->b[i]){
					pos[i] += 1;
					return pos;
				} else {
					pos[i] = 0;
				}
			}
			std::cout<<"trouble in iterator++ \n";
			return range->b;
		}
	};

	iterator begin(){
		iterator i;
		i.pos = a;
		i.range = this;
		return i;
	}
	iterator end(){
		return iterator();
	}	
};

#endif