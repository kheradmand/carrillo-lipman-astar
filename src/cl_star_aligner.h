/*
 * efficient_aligner.h
 *
 *  Created on: Nov 30, 2018
 *      Author: ali
 */

#ifndef SRC_CL_STAR_ALIGNER_H_
#define SRC_CL_STAR_ALIGNER_H_

#include "aligner.h"
#include "carrillo_lipman.h"
#include "dp_aligner.h"
#include <set>
#include <map>
#include <queue>

/*
 * Implementation of multiple sequence alignment using Carrillo-Lipman speed up method
 * Using tree alignment on center star (a 2-approximation of the optimal solution) to get the initial alignment score
 */

class cl_star_aligner_t : public aligner_t {
private:


	struct node_t {
		index_t index;
		score_t f;
		permutation_t parent;
	};

	struct node_compare_t {
		bool operator()(const node_t& a, const node_t& b) const {
			if (a.f == b.f){
				if (a.index == b.index){
					return a.parent > b.parent;
				}
				return a.index > b.index;
			}
			return a.f > b.f ;
		}
	};

	typedef std::priority_queue<node_t, std::vector<node_t>, node_compare_t> openset_t;
	typedef std::unordered_map<index_t, index_t> closedset_t; //= closed set



	score_t get_star_alignment_score(const sequences_t& seqs, const std::size_t c){
		dp_aligner_t dp_aligner(scoring);
		sequences_t final_alignment{seqs[c]};
		//std::cerr << "center is " << seqs[c] << std::endl;
		for (auto i = 0 ; i < seqs.size(); i++){
			if (i == c)
				continue;
			const auto center = final_alignment[0];
			const auto alignment = dp_aligner.get_alignment({center, seqs[i]});
			//std::cerr << "alignment of  " << center << " (center) and " << seqs[i] << " is " <<  alignment.first << std::endl;
			const auto aligned_center = alignment.first[0];
			auto old_p = 0;
			sequences_t new_alignment(final_alignment.size() + 1);
			for (auto p = 0; p < aligned_center.size(); p++){
				bool gap = false;
				if (aligned_center[p] != center[old_p]){
					assert(aligned_center[p] == GAP);
					gap = true;
				}
				for (auto j = 0; j < final_alignment.size(); j++){
					new_alignment[j] += (gap ? GAP : final_alignment[j][old_p]);
				}
				if (not gap){
					old_p++;
				}
			}
			new_alignment[final_alignment.size()] = alignment.first[1];
			final_alignment = new_alignment;
			//std::cerr << "alignment so far is " << final_alignment << std::endl;
		}
		score_t ret = 0;
		assert(final_alignment.size() == seqs.size());
		chars_t chars(final_alignment.size());
		for (auto p = 0; p < final_alignment[0].length(); p++){
			for (auto i = 0; i < final_alignment.size(); i++){
				chars[i] = final_alignment[i][p];
			}
			ret += scoring.sp_score(chars);
		}

		std::cerr << "star alignment score is  " << ret << std::endl;
		return ret;
	}


public:
	cl_star_aligner_t(scoring_function_t& _scoring) : aligner_t(_scoring){
	}

	std::pair<sequences_t,score_t> get_alignment(const sequences_t& seqs) {
		carrillo_lipman_score_t carillo_lipman_score(seqs, scoring);
		const auto c = carillo_lipman_score.find_center();
		const auto z = get_star_alignment_score(seqs, c);

		std::cerr << "starting to align" << std::endl;
		position_util_t pu(seqs);
		openset_t openset;
		closedset_t parents;
		auto dims = seqs.size();

		assert(dims < 64);

		const auto first_index = pu.first_index();
		const auto last_index = pu.last_index();


		openset.emplace(node_t{first_index, 0, 0});

		score_t best_score = INF;

		chars_t chars(dims);
		positions_t pos(dims);
		std::size_t visited = 0;
		while (true) {
			++visited;
			assert(not openset.empty());
			const auto node = openset.top();
			openset.pop();
			if (parents.find(node.index) != parents.end())
				continue;
			parents.emplace(node.index, node.parent);

			//std::cerr << "at index " << node.index << " with f " << node.f << std::endl;

			if (node.index == last_index){
				best_score = node.f;
				break;
			}

			pu.convert_back(node.index, pos);
			//std::cerr << "at pos " << pos << std::endl;

			const auto cl_score = carillo_lipman_score.get_score(seqs, pos);
			if (node.f + cl_score > z){
				//std::cerr << " skipping the node because " << node.f << " + " << cl_score << " > " << z << std::endl;
				continue;
			}

			for (permutation_t permutation = 1; permutation < (1 << dims); permutation++){
				index_t neighbor_index;
				bool valid;
				std::tie(valid, neighbor_index) = pu.get_forward_neighbor_by_index(node.index, permutation, chars, seqs);
				//std::cerr << "valid: " << valid << " neighbor index " << neighbor_index << std::endl;

				if (not valid or parents.find(neighbor_index) != parents.end())
					continue;

				//std::cerr << "at neighbor index " << neighbor_index << std::endl;
				//{positions_t neighbor_pos; pu.convert_back(neighbor_index, neighbor_pos); std::cerr << "neigbor pos " << neighbor_pos << std::endl;}

				const auto sp_score = scoring.sp_score(chars);
				openset.emplace(node_t{neighbor_index, node.f + sp_score, permutation});
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



#endif /* SRC_CL_STAR_ALIGNER_H_ */
