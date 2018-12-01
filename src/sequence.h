/*
 * sequence.h
 *
 *  Created on: Nov 24, 2018
 *      Author: ali
 */

#ifndef SRC_SEQUENCE_H_
#define SRC_SEQUENCE_H_

#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include "score.h"

typedef std::string sequence_t;
typedef std::vector<sequence_t> sequences_t;

typedef std::vector<std::size_t> positions_t;

typedef std::size_t index_t;
typedef unsigned long long permutation_t;

struct position_util_t{
	std::vector<size_t> lens;
	std::vector<size_t> cumulative_lens;
	enum direction_t {FORWARD , BACKWARD};

	position_util_t(const sequences_t& seqs){
		cumulative_lens.push_back(1);

		for (const auto& s : seqs){
			lens.push_back(s.size());
			cumulative_lens.push_back(cumulative_lens.back() * (s.size() + 1));
		}
	}

	index_t convert(const positions_t& pos){
		assert(pos.size() == lens.size());
		std::size_t ret = 0;
		for (auto i = 0 ; i < lens.size(); i++){
			ret += (pos[i]) * cumulative_lens[i];
		}
		return ret;
	}

	void convert_back(index_t index, positions_t& pos){
		pos.resize(lens.size());
		for (int d = lens.size() - 1; d >= 0; d--){
			pos[d] = index / cumulative_lens[d];
			index %= cumulative_lens[d];
		}
	}

	bool get_neighbor_by_pos(positions_t& pos, const permutation_t permutation, const direction_t direction, chars_t& chars){
		for (auto d = 0; d < lens.size(); d++){
			if (permutation & (1 << d)){
				if (pos[d] == 0){
					return false;
				} else {
					pos[d]--;
				}
			} else {
				chars[d] = GAP;
			}
		}
		return true;
	}

	std::pair<bool, index_t> get_backward_neighbor_by_index(index_t index, const permutation_t permutation, chars_t& chars, const sequences_t& seqs){
		index_t neighbor = index;
		for (int d = lens.size() - 1; d >= 0; d--){
			if (permutation & (1 << d)){
				if (cumulative_lens[d] > index){
					return {false, 0};
				} else {
					neighbor -= cumulative_lens[d];
					const auto pos_d = index/cumulative_lens[d];
					chars[d] = seqs[d][pos_d - 1];
				}
			} else {
				chars[d] = GAP;
			}
			index %= cumulative_lens[d];
		}
		return {true, neighbor};
	}

	std::pair<bool, index_t> get_forward_neighbor_by_index(index_t index, const permutation_t permutation, chars_t& chars, const sequences_t& seqs){
		index_t neighbor = index;
		for (int d = lens.size() - 1; d >= 0; d--){
			if (permutation & (1 << d)){
				const auto pos_d = index/cumulative_lens[d];
				if (pos_d >= lens[d]){
					return {false, 0};
				} else {
					neighbor += cumulative_lens[d];
					chars[d] = seqs[d][pos_d];
				}
			} else {
				chars[d] = GAP;
			}
			index %= cumulative_lens[d];
		}
		//std::cerr << "get_forward_neigbor_by_index for " << index << " with perm " << permutation  << " is " <<  neighbor  << std::endl;
		return {true, neighbor};
	}





	std::size_t space_size(){
		return cumulative_lens.back();
	}

	positions_t first_pos(){
		return positions_t(lens.size(), 0);
	}

	positions_t last_pos(){
		return lens;
	}

	index_t first_index(){
		return 0;
	}

	index_t last_index(){
		return space_size() - 1;
	}



	void get_next_pos(positions_t& pos){
		auto d = 0;
		while (d < lens.size()){
			if (pos[d] == lens[d]){
				pos[d] = 0;
				++d;
			} else {
				pos[d] ++;
				break;
			}
		}
		assert(d < lens.size());
	}
};


std::ostream& operator<<(std::ostream& out, const positions_t& pos){
	out << "(";
	for (auto i = 0; i < pos.size(); i++){
		out << pos[i] << ((i == pos.size() - 1) ? "" : ",");
	}
	out << ")";
	return out;
}

std::ostream& operator<<(std::ostream& out, const chars_t& chars){
	out << "{";
	for (auto i = 0; i < chars.size(); i++){
		out << chars[i] << ((i == chars.size() - 1) ? "" : ",");
	}
	out << "}";
	return out;
}

std::ostream& operator<<(std::ostream& out, const sequences_t& seqs){
	out << std::endl;
	for (auto i = 0; i < seqs.size(); i++){
		out << seqs[i] << std::endl;
	}
	return out;
}


#endif /* SRC_SEQUENCE_H_ */
