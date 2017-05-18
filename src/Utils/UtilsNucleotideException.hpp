/*
 * UtilsNucleotideException.hpp
 *
 *  Created on: 15 июня 2016 г.
 *      Author: thepunchy
 */

#ifndef UTILS_UTILSNUCLEOTIDEEXCEPTION_HPP_
#define UTILS_UTILSNUCLEOTIDEEXCEPTION_HPP_

#include <string>

class UtilsNucleotideException : public std::exception {
private:
	std::string msg;
public:
	UtilsNucleotideException(std::string msg = "Unhandled exception") : msg(msg) {};
	virtual const char *what() const throw()
	{
		return msg.c_str();
	}
	virtual ~UtilsNucleotideException() {};
};



#endif /* UTILS_UTILSNUCLEOTIDEEXCEPTION_HPP_ */
