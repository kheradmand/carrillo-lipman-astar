/*
 * sequence.h
 *
 *  Created on: Nov 24, 2018
 *      Author: ali
 */

#ifndef SRC_SEQUENCE_H_
#define SRC_SEQUENCE_H_

#include <string>
#include <vector>
#include <iostream>

typedef std::string sequence_t;
typedef std::vector<sequence_t> sequences_t;

typedef std::vector<std::size_t> positions_t;
typedef std::vector<char> chars_t;

std::ostream& operator<<(std::ostream& out, const positions_t pos){
	out << "(";
	for (auto i = 0; i < pos.size(); i++){
		out << pos[i] << ((i == pos.size() - 1) ? "" : ",");
	}
	out << ")";
	return out;
}

std::ostream& operator<<(std::ostream& out, const chars_t chars){
	out << "{";
	for (auto i = 0; i < chars.size(); i++){
		out << chars[i] << ((i == chars.size() - 1) ? "" : ",");
	}
	out << "}";
	return out;
}


#endif /* SRC_SEQUENCE_H_ */
