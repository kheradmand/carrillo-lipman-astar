/*
 * efficient_aligner.h
 *
 *  Created on: Nov 24, 2018
 *      Author: ali
 */

#ifndef SRC_ASTAR_ALIGNER_H_
#define SRC_ASTAR_ALIGNER_H_

#include "aligner.h"
#include "carrillo_lipman.h"
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

	struct node_compare_t {
		bool operator()(const node_t& a, const node_t& b) const {
			if (a.f + a.g == b.f + b.g){
				if (a.index == b.index){
					return a.parent < b.parent;
				}
				return a.index < b.index;
			}
			return a.f + a.g < b.f + b.g;
		}
	};

	typedef std::priority_queue<node_t, std::vector<node_t>, node_compare_t> openset_t;

	typedef std::unordered_map<index_t, index_t> closedset_t; //= closed set


public:
	astar_aligner_t(scoring_function_t& _scoring) : aligner_t(_scoring){}

	std::pair<sequences_t,score_t> get_alignment(const sequences_t& seqs) {
		carrillo_lipman_score_t carillo_lipman_score(seqs, scoring);
		std::cerr << "starting to align" << std::endl;
		position_util_t pu(seqs);
		openset_t openset;
		closedset_t parents;
		auto dims = seqs.size();

		assert(dims < 64);


		const auto& first_pos = pu.first_pos();
		const auto& last_pos = pu.last_pos();

		const auto first_index = pu.first_index();
		const auto last_index = pu.last_index();


		openset.emplace(node_t{first_index, 0, 0, 0});

		score_t best_score = MINUS_INFINITY;

		chars_t chars(dims);
		positions_t neighbor_pos(dims);
		std::size_t visited = 0;
		while (true) {
			++visited;
			assert(not openset.empty());
			const auto node = openset.top();
			openset.pop();
			if (parents.find(node.index) != parents.end())
				continue;
			parents.emplace(node.index, node.parent);

			//std::cerr << "at index " << node.index << " with f " << node.f << " and g " << node.g << std::endl;

			if (node.index == last_index){
				best_score = node.f;
				break;
			}

			for (permutation_t permutation = 1; permutation < (1 << dims); permutation++){
				index_t neighbor_index;
				bool valid;
				std::tie(valid, neighbor_index) = pu.get_forward_neighbor_by_index(node.index, permutation, chars, seqs);
				//std::cerr << "valid: " << valid << " neighbor index " << neighbor_index << std::endl;

				if (not valid or parents.find(neighbor_index) != parents.end())
					continue;

				//std::cerr << "at neighbor index " << neighbor_index << std::endl;

				pu.convert_back(neighbor_index, neighbor_pos);
				const auto cl_score = carillo_lipman_score.get_score(seqs, neighbor_pos);
				const auto sp_score = scoring.sp_score(chars);
				openset.emplace(node_t{neighbor_index, node.f + sp_score, cl_score, permutation});
			}
		}

		assert (parents.find(last_index) != parents.end());

		std::cerr << "open set " << openset.size() << " closed set " << parents.size() << " sum " << openset.size() + parents.size() << std::endl;
		std::cerr << "space size " << pu.space_size() << " visited " << visited <<  std::endl;

		sequences_t alignment{dims};

		//std::cerr << "best score is " << best_score << std::endl;
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



#endif /* SRC_ASTAR_ALIGNER_H_ */
