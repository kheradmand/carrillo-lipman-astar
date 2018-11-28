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

/*
 * Implementation of multiple sequence alignment using A* with Carillo-Lipman heuristic
 */
class efficient_aligner : public aligner_t {

private:
	struct
	std::set
	score_t get_optimal_pairwise_alignment_score(const sequence_t& s1, const std::size_t p1, const sequence_t& s2, const std::size_t p2){
		std::vector<score_t> best_scores(s2.size() - p2 + 1);
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

public:
	efficient_aligner(scoring_function_t& _scoring) : aligner_t(_scoring){}

	std::pair<sequences_t,score_t> get_alignment(const sequences_t& seqs) {

	}

};



#endif /* SRC_EFFICIENT_ALIGNER_H_ */
