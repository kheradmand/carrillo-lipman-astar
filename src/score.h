/*
 * score.h
 *
 *  Created on: Nov 24, 2018
 *      Author: ali
 */

#ifndef SRC_SCORE_H_
#define SRC_SCORE_H_

#include <limits>
#include <unordered_map>
#include "sequence.h"

const char GAP = '-';

typedef long long score_t;

typedef std::vector<char> chars_t;

const score_t MINUS_INFINITY = std::numeric_limits<score_t>::min();

struct scoring_function_t {
	score_t scoring_matrice[256][256];

	scoring_function_t(){
		memset(scoring_matrice, -1, sizeof(scoring_matrice));
		for (auto i = 0 ; i < 256; i++){
			scoring_matrice[i][i] = 1;
		}
		scoring_matrice[GAP][GAP] = 0;
	}

	score_t sp_score(const chars_t& chars){
		score_t ret = 0;
		for (auto i = 0; i < chars.size(); i++){
			for (auto j = i + 1; j < chars.size(); j++){
				if (scoring_matrice[chars[i]][chars[j]] == MINUS_INFINITY)
					return MINUS_INFINITY;
				ret += scoring_matrice[chars[i]][chars[j]];
			}
		}
		//std::cout << "returning " << ret << std::endl;
		return ret;
	}

	inline score_t score(const char c1, const char c2){
		//std::cout << "checking score of " << c1 << "," << c2 << std::endl;
		return scoring_matrice[c1][c2];
	}
};





#endif /* SRC_SCORE_H_ */
