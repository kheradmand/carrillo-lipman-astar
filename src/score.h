/*
 * score.h
 *
 *  Created on: Nov 24, 2018
 *      Author: ali
 */

#ifndef SRC_SCORE_H_
#define SRC_SCORE_H_

#include <limits>
#include <fstream>
#include <unordered_map>
#include "sequence.h"

const char GAP = '-';

typedef long long score_t;

typedef std::vector<char> chars_t;

//const score_t MINUS_INFINITY = std::numeric_limits<score_t>::min();
const score_t INF = std::numeric_limits<score_t>::max();

struct scoring_function_t {
	score_t scoring_matrice[256][256];

	void reset(score_t v){
		//memset to 1 does do not work since it works per byte
		for (auto i = 0 ; i < 256; i++){
			for (auto j = 0 ; j < 256; j++){
				scoring_matrice[i][j] = v;
			}
		}
	}

	void set_score(const char c1, const char c2, score_t s){
		//std::cout << "setting " << c1 << " " << c2 << " to " << s << std::endl;
		scoring_matrice[c1][c2] = s;
		scoring_matrice[c2][c1] = s;
	}

	scoring_function_t(){
		reset(1);
		for (auto i = 0 ; i < 256; i++){
			scoring_matrice[i][i] = 0;
		}
		scoring_matrice[GAP][GAP] = 0;
	}

	scoring_function_t(const std::string& filename){
		//reading the scores from file
		reset(INF);
		std::ifstream fin(filename);
		assert(fin.is_open());
		chars_t chars;
		char c;
		//std::cout << "reading "<< std::endl;
		while (fin >> c){
			if (c == '#') { //comment;
				std::string temp;
				std::getline(fin, temp);
				//std::cout << "skipping " << temp << std::endl;
			} else {
				chars.push_back(c);
				for (auto cp : chars){
					score_t s;
					fin >> s;
					set_score(c, cp, s);
				}
			}
		}
	}



	score_t sp_score(const chars_t& chars){
		score_t ret = 0;
		for (auto i = 0; i < chars.size(); i++){
			for (auto j = i + 1; j < chars.size(); j++){
				if (scoring_matrice[chars[i]][chars[j]] == INF)
					return INF;
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
