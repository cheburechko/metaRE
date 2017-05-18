/*
 * UtilsKmerTooBigException.hpp
 *
 *  Created on: 15 июня 2016 г.
 *      Author: thepunchy
 */

#ifndef UTILS_UTILSKMERTOOBIGEXCEPTION_HPP_
#define UTILS_UTILSKMERTOOBIGEXCEPTION_HPP_

#include <string>

class UtilsKmerTooBigException : public std::exception {
private:
	std::string msg;
public:
	UtilsKmerTooBigException(std::string msg = "Unhandled exception") : msg(msg) {};
	virtual const char *what() const throw()
	{
		return msg.c_str();
	}
	virtual ~UtilsKmerTooBigException() {};
};


#endif /* UTILS_UTILSKMERTOOBIGEXCEPTION_HPP_ */
