/*
 * naive_aligner.h
 *
 *  Created on: Nov 24, 2018
 *      Author: ali
 */

#ifndef SRC_NAIVE_ALIGNER_H_
#define SRC_NAIVE_ALIGNER_H_

#include "aligner.h"
#include <vector>
#include <cassert>
#include <algorithm>

/*
 * Implementation of the naive dynamic programming approach for multiple sequence alignment
 * O(k^2.2^k.n^k) in time and O(n^k) in memory for k sequences of length n
 */
class naive_aligner_t : public aligner_t {
private:

	struct position_util_t{
		std::vector<size_t> lens;
		std::vector<size_t> cumulative_lens;
		position_util_t(const sequences_t& seqs){
			cumulative_lens.push_back(1);

			for (const auto& s : seqs){
				lens.push_back(s.size());
				cumulative_lens.push_back(cumulative_lens.back() * (s.size() + 1));
			}
		}

		std::size_t convert(const positions_t pos){
			assert(pos.size() == lens.size());
			std::size_t ret = 0;
			for (auto i = 0 ; i < lens.size(); i++){
				ret += (pos[i]) * cumulative_lens[i];
			}
			return ret;
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



	typedef unsigned long long permutation_t;

public:
	naive_aligner_t(scoring_function_t& _scoring) : aligner_t(_scoring){}

	std::pair<sequences_t,score_t> get_alignment(const sequences_t& seqs) {
		//std::cerr << "starting to align" << std::endl;
		position_util_t pu(seqs);
		auto dims = seqs.size();
		score_t* best_scores = new score_t[pu.space_size()];
		permutation_t* parents = new permutation_t[pu.space_size()];

		assert(dims < 64);

		const auto&  first_pos = pu.first_pos();
		const auto& last_pos = pu.last_pos();

		positions_t pos{first_pos};
		best_scores[pu.convert(pos)] = 0;


		chars_t chars_base(dims);
		do {
			pu.get_next_pos(pos);
			//std::cerr << "at possition " << pos << std::endl;
			const auto pos_index = pu.convert(pos);
			best_scores[pos_index] = MINUS_INFINITY;
			parents[pos_index] = 0;

			for (auto d = 0; d < dims; d++){
				chars_base[d] = seqs[d][pos[d] - 1];
			}

			for (permutation_t permutation = 1; permutation < (1 << dims); permutation++){
				positions_t neighbor{pos};
				auto chars = chars_base;
				bool valid = true;

				for (auto d = 0; d < dims; d++){
					if (permutation & (1 << d)){
						if (neighbor[d] == 0){
							valid = false;
							break;
						} else {
							neighbor[d]--;
						}
					} else {
						chars[d] = GAP;
					}
				}

				if (valid){
					//std::cerr << "neighbor is " << neighbor << "chars are " << chars << std::endl;
					const auto nei_score = best_scores[pu.convert(neighbor)] + scoring.sp_score(chars);
					if (best_scores[pos_index] < nei_score){
						best_scores[pos_index] = nei_score;
						parents[pos_index] = permutation;
					}
				}

			}


		} while (pos != last_pos);

		score_t best_score = best_scores[pu.convert(last_pos)];
		sequences_t alignment{dims};

		//std::cerr << "best score is " << best_scores[pu.convert(last_pos)] << std::endl;

		//std::cerr << "traceback to construct an alignment" << std::endl;
		//traceback
		while (pos != first_pos){
			//std::cerr << "at possition " << pos << std::endl;
			auto permutation = parents[pu.convert(pos)];
			assert(permutation != 0);
			for (auto d = 0; d < dims; d++){
				if (permutation & (1 << d)){
					assert (pos[d] > 0);
					alignment[d] += seqs[d][pos[d] - 1];
					pos[d]--;
				} else {
					alignment[d] += GAP;
				}
			}
		}

		for (auto d = 0; d < dims; d++){
			std::reverse(alignment[d].begin(), alignment[d].end());
		}

		delete best_scores;
		delete parents;

		return {alignment, best_score};
	}

};



#endif /* SRC_NAIVE_ALIGNER_H_ */
