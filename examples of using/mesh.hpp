#ifndef MESH_H
#define MESH_H
#include "vec.hpp"
#include <vector>
#include <fstream>

template <typename T, int D> class Mesh{
	std::vector<T> data;
	Ind<D> box_, mul;
	int bc_type; // 1 - продолжение, 2 - переодика
	Vec<D> step, coord0;
public:

	const Ind<D> & box() const { return box; }

	T& operator [] (const Ind<D> &pos){
		if((1<<0)&bc_type){
			Ind<D> new_pos;
			for (int i = 0; i < D; i++){
				new_pos[i] = pos[i];	
				if ( pos[i] >= box_[i])
					new_pos[i] = box_[i] - 1;
				if (pos[i] < 0)
					new_pos[i] = 0;
			}
			return data[new_pos * mul];
		} else { 
			Ind<D> new_pos; 
			for(int i = 0; i < D; i++){
				if(pos[i] < 0)
					new_pos[i] = (pos[i] + (1 - pos[i]/box_[i]) * box_[i]) % box_[i];
				new_pos[i] = pos[i] % box_[i];
			}
			return data[new_pos * mul];
		} 
	} 

	T operator [](size_t f) const {return data[f];}
	T& operator [](size_t f){return data[f];}

	Mesh(const Ind<D> &new_box, Vec<D> new_step, Vec<D> new_coord0, int new_bc = 1){
		step = new_step;
		coord0 = new_coord0;
		bc_type = new_bc;
		box_ = new_box;
		size_t sz = box_[0]; mul[0] = 1;
		for(int i=1; i<D; i++){ sz *= box_[i]; mul[i] = mul[i-1]*box_[i-1]; }
		data.resize(sz);
	}

	T& operator () (const Vec<D> &coord){
		Ind<D> pos;
		for(int i = 0; i < D; i++){
			pos[i] = floor((coord[i] - coord0)/step);
		}
		return data[pos];
	}

	size_t prod (){
		size_t value = 1;
		for (int i = 0; i < D; i++)
			value *= box_[i];
		return value;
	}

	void out2dat(const char *path){
		std::ofstream out;
		out.open(path);
		for(size_t i = 0, sz = (*this).prod(); i < sz; i++ ){
			Ind<D> pos = i % box_;
			out << pos  <<' ' << (*this)[pos] << '\n';
		}
		out.close();
	}

	void swap_data(Mesh<T,D> &m){
		this->data.swap(m.data);
	}

};

template <int D> struct SmartInd{
	size_t f;
	Ind<D> pos, box; 
	Ind<2*D> d_out, d_in;
	bool bound;
	SmartInd (const Ind<D> &box_, int bc_type): f(0), pos(), box(box_), bound(true){
		d_in[0] = 1;
		d_in[1] = -1;
		for(int i = 1; i < D; i++){
			d_in[2*i] = d_[2*(i-1)] * box[i-1];
			d_in[2*i + 1] = - d_in[2*i]; 
		}
		if((1<<0)&bc_type){
			for(int i = 0; i < 2*D; i++)
				d_out[i] = 0;
		} else {
			int tmp1 = 1;
			int tmp2 = 1;
			for(int i = 0; i < D; i++){
				tmp1 *= box[i];
				d_out[2*i] = tmp1 - tmp2;
				d_out[2*i + 1] = -d_out[2*i];
				tmp2 *= box[i] 
			}
		}

	}


};

#endif