/*
 * efficient_aligner.h
 *
 *  Created on: Nov 24, 2018
 *      Author: ali
 */

#ifndef SRC_EFFICIENT_ALIGNER_H_
#define SRC_EFFICIENT_ALIGNER_H_

#include "aligner.h"
#include <set>
#include <map>
#include <queue>

/*
 * Implementation of multiple sequence alignment using A* with Carillo-Lipman heuristic
 */
class astar_aligner_t : public aligner_t {
private:

	struct node_t {
		index_t index;
		score_t f,g;
		permutation_t parent;
	};

	typedef std::priority_queue<node_t> openset_t;

	typedef std::unordered_map<index_t, index_t> closedset_t; //= closed set

	score_t get_optimal_pairwise_alignment_score(const sequence_t& s1, const std::size_t p1, const sequence_t& s2, const std::size_t p2){
		std::vector<score_t> best_scores(s2.size() + 1 - p2);
		best_scores[0] = 0;
		// i = 0;
		for (auto j = 1; j <= s2.size() - p2; j++)
			best_scores[j] = best_scores[j-1] + scoring.score(GAP, s2[j + p2 - 1]);
		for (auto i = 1 ; i <= s1.size() - p1; i++){
			//j=0
			score_t prev_row_best_score = best_scores[0];
			best_scores[0] = best_scores[0] + scoring.score(s1[i + p1 - 1], GAP);
			for (auto j = 1; j <= s2.size(); j++){
				const score_t _c = best_scores[i - 1] + scoring.score(GAP, s2[j + p2 - 1]);
				const score_t c_ = best_scores[i] + scoring.score(s1[i + p1 - 1], GAP);
				const score_t cc = prev_row_best_score + scoring.score(s1[i + p1 - 1], s2[j + p2 - 1]);
				prev_row_best_score = best_scores[i];
				best_scores[i] = std::max(cc, std::max(_c,c_));
			}
		}
		return best_scores[s2.size()];
	}

	score_t get_carrillo_lipman_score(const sequences_t& seqs, const positions_t& pos){
		score_t ret;
		for (auto i = 0; i < seqs.size(); i++){
			for (auto j = 0; j < seqs.size(); j++){
				ret += get_optimal_pairwise_alignment_score(seqs[i], pos[i] + 1, seqs[j], pos[j] + 1);
			}
		}
		return ret;
	}

public:
	astar_aligner_t(scoring_function_t& _scoring) : aligner_t(_scoring){}

	std::pair<sequences_t,score_t> get_alignment(const sequences_t& seqs) {
		//std::cerr << "starting to align" << std::endl;
		position_util_t pu(seqs);
		openset_t openset;
		closedset_t parents;
		auto dims = seqs.size();

		assert(dims < 64);


		const auto& first_pos = pu.first_pos();
		const auto& last_pos = pu.last_pos();

		const auto first_index = pu.first_index();
		const auto last_index = pu.last_index();


		openset.emplace(first_index, 0, get_carrillo_lipman_score(seqs, first_pos), 0);

		score_t best_score = MINUS_INFINITY;

		chars_t chars(dims);
		positions_t neighbor_pos(dims);
		while (true) {
			assert(not openset.empty());

			const auto node = openset.top();
			openset.pop();
			if (parents.find(node.index) != parents.end())
				continue;
			parents.emplace(node.index, node.parent);

			if (node.index == last_index){
				best_score = node.f;
				break;
			}

			for (permutation_t permutation = 1; permutation < (1 << dims); permutation++){
				index_t neighbor_index;
				bool valid;
				std::tie(valid, neighbor_index) = pu.get_forward_neighbor_by_index(node.index, permutation, chars, seqs);

				if (not valid or parents.find(neighbor_index) != parents.end())
					continue;

				pu.convert_back(neighbor_index, neighbor_pos);
				const auto cl_score = get_carrillo_lipman_score(seqs, neighbor_pos);
				const auto sp_score = scoring.sp_score(chars);
				openset.emplace(neighbor_index, node.f + sp_score, cl_score, permutation);
			}
		}

		assert (parents.find(last_index) != parents.end());

		sequences_t alignment{dims};

		//std::cerr << "best score is " << best_scores[pu.convert(last_pos)] << std::endl;

		//std::cerr << "traceback to construct an alignment" << std::endl;
		//traceback
		index_t index = last_index;
		while (index != first_index){
			//std::cerr << "at possition " << pos << std::endl;
			const auto permutation = parents.at(index);
			assert(permutation != 0);
			bool valid;
			std::tie(valid, index) = pu.get_backward_neighbor_by_index(index, permutation, chars, seqs);
			assert(valid);

			for (auto d = 0; d < dims; d++){
				alignment[d] += chars[d];
			}
		}

		for (auto d = 0; d < dims; d++){
			std::reverse(alignment[d].begin(), alignment[d].end());
		}

		return {alignment, best_score};
	}

};



#endif /* SRC_EFFICIENT_ALIGNER_H_ */
