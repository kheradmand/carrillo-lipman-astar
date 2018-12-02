/*
 * main.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: ali
 */

#include <iostream>

#include "dp_aligner.h"
#include "astar_aligner.h"
#include "cl_star_aligner.h"
using namespace std;


void get_and_print_alignment(const sequences_t seqs, aligner_t& aligner){
	const auto res = aligner.get_alignment(seqs);

	cout << "score: " << res.second << endl;
	for (const auto& s : res.first){
		cout << s << endl;
	}
}

int main(){
	sequences_t seqs;
	string s;
	while (cin >> s){
		seqs.push_back(s);
	}

	scoring_function_t score;


	dp_aligner_t dp_aligner{score};
	astar_aligner_t astar_aligner{score};
	cl_star_aligner_t cl_star_aligner{score};

//	cout << "with dp_aligner:" << endl;
//	get_and_print_alignment(seqs, dp_aligner);

//	cout << "with astar_aligner:" << endl;
//	get_and_print_alignment(seqs, astar_aligner);

	cout << "with cl_star_aligner:" << endl;
	get_and_print_alignment(seqs, cl_star_aligner);


}


