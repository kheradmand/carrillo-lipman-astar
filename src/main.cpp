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

int main(int argc, char* argv[]){
	string aliger_name = "astar";
	string scoring_file = "";

	for (auto i = 1; i < argc; i++){
		string arg(argv[i]);
		if (arg == "-h" or arg == "--help"){
			cout << "usage: " << argv[0] << " [OPTIONS]\n" <<
			"\t OPTIONS:\n" <<
					"\t\t" << "-h,--help\t\t" << "print this help message\n" <<
					"\t\t" << "-a,--aligner <name>\t" << "set aligner type, one of\n\t\t\tastar: A* with Carrillo-Lipman (CL) heuristic (default)\n\t\t\tcl-cs: original CL method, z from center star alignment\n\t\t\tdp: naive dynamic programming\n" <<
					"\t\t" << "-s,--score <file>\t" << "uses the provided scoring file.\n\t\t\tdefault: 1 for mismatch and gap, 0 otherwise\n" <<
					endl;
			return 0;
		} else if (arg == "-a" or arg == "--aligner"){
			aliger_name = argv[++i];
		} else if (arg == "-s" or arg == "--score"){
			scoring_file = argv[++i];
		} else {
			cout << "invalid option";
			return 1;
		}
	}


	scoring_function_t* score;
	if (scoring_file == ""){
		score = new scoring_function_t();
	} else {
		score = new scoring_function_t(scoring_file);
	}

	aligner_t* aligner;
	if (aliger_name == "astar"){
		aligner = new astar_aligner_t{*score};
	} else if (aliger_name == "cl-cs"){
		aligner = new cl_star_aligner_t{*score};
	} else if (aliger_name == "dp"){
		aligner = new dp_aligner_t{*score};
	} else {
		cout << "invalid aligner name " << aliger_name << endl;
		return 1;
	}




	sequences_t seqs;
	string s;
	while (cin >> s){
		seqs.push_back(s);
	}
	get_and_print_alignment(seqs, *aligner);

	delete aligner;
	delete score;

	return 0;

}


