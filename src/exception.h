/**
 * This file is part of the MultiDimSparseArray library
 *
 * @license  BSD-3
 * @author   Abhilekh Agarwal
 */


#ifndef __SPARSEMATRIX_EXCEPTIONS_H__
#define __SPARSEMATRIX_EXCEPTIONS_H__

#include <exception>
#include <string>
#include <stdexcept>

class Exception : public std::exception {
public:
    explicit Exception(const std::string & msg) : exception(), message(msg)
    {}
    virtual ~Exception(void) throw ()
    {}
    inline std::string getMessage(void) const {
        return this->message;
    }
protected:
    std::string message;
};

class InvalidStateException : public Exception
{
public:
    InvalidStateException(const std::string & msg) : Exception(msg)
    {}
};


class InvalidDimensionsException : public Exception
{
public:
    InvalidDimensionsException(const std::string & msg) : Exception(msg)
    {}
};


class InvalidCoordinatesException : public Exception {
public:
    InvalidCoordinatesException(const std::string & msg) : Exception(msg)
    {}
};

#endif
