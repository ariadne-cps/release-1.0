/***************************************************************************
 *            utility/macros.h
 *
 *  Copyright 2013-14  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file utility/macros.h
 *  \brief
 */



#ifndef ARIADNE_MACROS_H
#define ARIADNE_MACROS_H

#include <sstream>
#include <stdexcept>
#include "exceptions.h"

typedef std::string StringType;
typedef std::stringstream StringStream;

#define ARIADNE_USING_CONSTRUCTORS(Class,Base) \
    template<class T,typename std::enable_if<std::is_convertible<T,Base>::value,int>::type=0> \
    Class(const T& t) : Base(t) { } \
    template<class T,typename std::enable_if<std::is_constructible<T,Base>::value and not std::is_convertible<T,Base>::value,int>::type=0> \
    explicit Class(const T& t) : Base(t) { } \
    template<class ...Args> Class(Args&&... args) : Base(std::forward<Args>(args)...) { } \

#define ARIADNE_THROW(except,func,msg)          \
    { \
        StringStream ss; \
        ss << #except " in " << func << ": " << msg;    \
        throw except(ss.str()); \
    } \

#define ARIADNE_ASSERT(expression) \
    { \
        bool result = static_cast<bool>(expression); \
        if(!result) { \
            ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<__FUNCTION__,"Assertion `" << #expression << "' failed.\n"); \
        } \
    } \


#ifndef NDEBUG
#define ARIADNE_DEBUG_ASSERT_MSG(expression,error) \
    { \
        bool result = static_cast<bool>(expression); \
        if(!result) { \
            ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<ARIADNE_PRETTY_FUNCTION,"Assertion `" << #expression << "' failed.\n"<<"  "<<error<<"\n"); \
        } \
    } \

#else
#define ARIADNE_DEBUG_ASSERT_MSG(expression,error) \
    { }
#endif


#ifndef NDEBUG
#define ARIADNE_DEBUG_ASSERT(expression) \
    { \
        bool result = static_cast<bool>(expression); \
        if(!result) { \
            ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<__FUNCTION__,"Assertion `" << #expression << "' failed.\n"); \
        } \
    } \

#else
#define ARIADNE_DEBUG_ASSERT(expression) \
    { }
#endif


#define ARIADNE_PRECONDITION_MSG(expression,error)             \
    { \
        bool result = static_cast<bool>(expression); \
        if(!result) { \
            ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<ARIADNE_PRETTY_FUNCTION,"Precondition `" << #expression << "' failed.\n"<<"  "<<error<<"\n"); \
        } \
    } \

#define ARIADNE_PRECONDITION(expression)             \
    { \
        bool result = static_cast<bool>(expression); \
        if(!result) { \
            ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<ARIADNE_PRETTY_FUNCTION,"Precondition `" << #expression << "' failed.\n"); \
        } \
    } \

#ifndef NDEBUG
#define ARIADNE_DEBUG_PRECONDITION(expression) \
    { \
        bool result = static_cast<bool>(expression); \
        if(!result) { \
            ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<__FUNCTION__,"Precondition `" << #expression << "' failed.\n"); \
        } \
    } \

#else
#define ARIADNE_DEBUG_PRECONDITION(expression) \
    { }
#endif

#define ARIADNE_FAIL_MSG(error)             \
    { \
        ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<ARIADNE_PRETTY_FUNCTION,"ErrorTag "<<error<<"\n"); \
    } \

#define ARIADNE_ASSERT_MSG(expression,error)             \
    { \
        bool result = static_cast<bool>(expression); \
        if(!result) { \
            ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<ARIADNE_PRETTY_FUNCTION,"Assertion `" << #expression << "' failed.\n"<<"  "<<error<<"\n"); \
        } \
    } \

#define ARIADNE_ASSERT_EQUAL(expression1,expression2)    \
    { \
        bool result = static_cast<bool>((expression1) == (expression2));       \
        if(!result) { \
            ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<ARIADNE_PRETTY_FUNCTION,"Assertion `" << #expression1 << "==" << #expression2 << "' failed.\n"<<"  "<<expression1<<" != "<<expression2<<"\n"); \
        } \
    } \

#define ARIADNE_NOT_IMPLEMENTED                 \
    throw std::runtime_error(StringType("Not implemented: ")+ARIADNE_PRETTY_FUNCTION);

#define ARIADNE_DEPRECATED(fn,msg)          \
    static bool first_time=true; \
    if(first_time) { \
        first_time=false; \
        std::cerr << "DEPRECATED: Function " << #fn << " is deprecated. " << #msg << std::endl; \
    } \

#define ARIADNE_WARN(msg)          \
    {                                                                \
        std::cerr << "WARNING: " << msg << "" << std::endl;                \
    }
                                                                  \
#define ARIADNE_ERROR(msg)          \
    {                                                                \
        std::cerr << "ERROR: " << msg << "" << std::endl;                \
    }
                                                                  \
#if defined(linux) || defined(__linux) || defined(__linux__)
#define ARIADNE_PRETTY_FUNCTION __PRETTY_FUNCTION__
#elif defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#define ARIADNE_PRETTY_FUNCTION __FUNCTION__
#elif defined(darwin) || defined(__darwin) || defined(__darwin__)
#define ARIADNE_PRETTY_FUNCTION __PRETTY_FUNCTION__
#else
#define ARIADNE_PRETTY_FUNCTION ""
#endif


#endif // ARIADNE_MACROS_H
