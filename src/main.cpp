/*
 * main.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: ali
 */

#include <iostream>

#include "naive_aligner.h"
using namespace std;

int main(){
	string s1,s2;
	cin >> s1 >> s2;

	scoring_function_t score;
	dp_aligner_t naive_aligner{score};

	const auto res = naive_aligner.get_alignment({s1,s2});

	cout << "score: " << res.second << endl;
	for (const auto& s : res.first){
		cout << s << endl;
	}



}


