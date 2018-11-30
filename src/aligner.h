/*
 * aligner.h
 *
 *  Created on: Nov 24, 2018
 *      Author: ali
 */

#ifndef SRC_ALIGNER_H_
#define SRC_ALIGNER_H_

#include "sequence.h"



class aligner_t {
protected:
	scoring_function_t& scoring;


public:

	aligner_t(scoring_function_t& _scoring) : scoring(_scoring) {}

	//virtual score_t get_alignment_score(const sequences_t& seqs, positions_t starting_pos) = 0;
	virtual std::pair<sequences_t,score_t> get_alignment(const sequences_t& seqs) = 0;

	virtual ~aligner_t(){}

};



#endif /* SRC_ALIGNER_H_ */
