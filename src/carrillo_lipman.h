/*
 * carrillo_lipman.h
 *
 *  Created on: Nov 30, 2018
 *      Author: ali
 */

#ifndef SRC_CARRILLO_LIPMAN_H_
#define SRC_CARRILLO_LIPMAN_H_

#include "score.h"

class carrillo_lipman_score_t {
	scoring_function_t& scoring;
	score_t** scores;
	std::size_t dims;

public:
	carrillo_lipman_score_t(const sequences_t& seqs, scoring_function_t& _scoring) : scoring(_scoring){
		dims = seqs.size();
		assert(dims > 1);
		scores = new score_t*[dims * dims];
		memset(scores, 0, dims * dims * sizeof(score_t*));

		for (auto i = 0; i < dims; i++){
			for (auto j = i + 1; j < dims; j++){
				scores[ i * dims + j ] = compute_pairwise_backward_alignment_scores(seqs[i], seqs[j]);
			}
		}
	}

	~carrillo_lipman_score_t(){
		for (auto i = 0; i < dims * dims; i++){
			delete scores[i];
		}
		delete scores;
	}

	inline std::size_t convert(const sequence_t& a, int i, int j){
		return j * (a.size() + 1) + i;
	}

	score_t* compute_pairwise_backward_alignment_scores(const sequence_t& a, const sequence_t& b){
		score_t* best_scores = new score_t[(a.size() + 1) * (b.size() + 1)];

		best_scores[convert(a, a.size(), b.size())] = 0;

		//i = a.size()
		for (int j = b.size() - 1; j >= 0 ; j--)
			best_scores[convert(a, a.size(), j)] = best_scores[convert(a, a.size(), j + 1)] + scoring.score(GAP, b[j]);
		for (int i = a.size() - 1; i >= 0; i--){
			//j = b.size()
			best_scores[convert(a, i, b.size())] = best_scores[convert(a, i + 1, b.size())] + scoring.score(a[i], GAP);
			for (int j = b.size() - 1; j >= 0; j--){
				const auto cc = best_scores[convert(a, i + 1, j + 1)] + scoring.score(a[i], b[j]);
				const auto c_ = best_scores[convert(a, i + 1, j)] + scoring.score(a[i], GAP);
				const auto _c = best_scores[convert(a, i, j + 1)] + scoring.score(GAP, b[j]);
				best_scores[convert(a, i, j)] = std::max(cc, std::max(_c,c_));
			}
		}

		//print_best_scores(a, b, best_scores);

		return best_scores;
	}

	void print_best_scores(const sequence_t& a, const sequence_t& b, score_t* best_scores){
		for (int i = 0 ; i <= a.size(); i++){
			for (int j = 0; j <= b.size(); j++){
				std::cerr << best_scores[convert(a, i, j)] << "\t";
			}
			std::cerr << std::endl;
		}
		std::cerr << std::endl;
	}



	score_t get_score(const sequences_t& seqs, const positions_t& pos){
		//std::cerr << "pos is " << pos << std::endl;
		score_t ret = 0;
		for (auto i = 0; i < dims; i++){
			for (auto j = i + 1; j < dims; j++){
				const auto score = scores[ i * dims + j ][convert(seqs[i], pos[i], pos[j])];
				//std::cerr << "optimal alignment for " << seqs[i].substr(pos[i], seqs[i].length() - pos[i]) << " and  " << seqs[j].substr(pos[j], seqs[j].length() - pos[j]) << " is " << score  << std::endl;
				ret += score;
			}
		}
		return ret;
	}

	std::size_t find_center(){
		std::size_t center = 0;
		score_t best_score = MINUS_INFINITY;
		for (auto i = 0; i < dims; i++){
			score_t score = 0;
			for (auto j = 0 ; j < i; j++){
				score += scores[j * dims + i][0];
			}
			for (auto j = i + 1; j < dims; j++){
				score += scores[i * dims + j][0];
			}
			if (score > best_score){
				center = i;
				best_score = score;
			}
		}
		std::cerr << "center index is " << center << " with score " << best_score << std::endl;
		return center;
	}

};



#endif /* SRC_CARRILLO_LIPMAN_H_ */
