/*
 * PatternWrongSizeException.h
 *
 *  Created on: 15 июня 2016 г.
 *      Author: thepunchy
 */

#ifndef PATTERN_PATTERNWRONGSIZEEXCEPTION_HPP_
#define PATTERN_PATTERNWRONGSIZEEXCEPTION_HPP_

class PatternWrongSizeException : public std::exception {
private:
	std::string msg;
public:
	PatternWrongSizeException(std::string msg = "Unhandled exception") : msg(msg) {};
	virtual const char *what() const throw()
	{
		return msg.c_str();
	}
	virtual ~PatternWrongSizeException() {};
};

#endif /* PATTERN_PATTERNWRONGSIZEEXCEPTION_HPP_ */
