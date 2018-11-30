/*
 * main.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: ali
 */

#include <iostream>

#include "dp_aligner.h"
#include "astar_aligner.h"
using namespace std;

int main(){
	sequences_t seqs;
	string s;
	while (cin >> s){
		seqs.push_back(s);
	}

	scoring_function_t score;
	//dp_aligner_t aligner{score};
	astar_aligner_t aligner{score};

	const auto res = aligner.get_alignment(seqs);

	cout << "score: " << res.second << endl;
	for (const auto& s : res.first){
		cout << s << endl;
	}
}


