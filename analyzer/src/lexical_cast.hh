// -*- C++ -*-

// Author: Tomonori Takahashi

#ifndef HDDAQ__LEXICAL_CAST_H
#define HDDAQ__LEXICAL_CAST_H

///////////////////////////////////////////////////////////////////////////////
//
// cast from/to "(const) char*" has not been supported.
// plaese use std::string instead of (const) char*, etc.
//
// If you are interested in the full implemetation of "lexical_cast",
// please see Boost C++ library. 
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <limits>
#include <string>
#include <sstream>
#include <map>

namespace hddaq
{

//______________________________________________________________________________
  template <typename Target, typename Source>
  struct lexical_cast_impl
  {
    static Target 
    lexical_cast(Source arg)
    {
      std::stringstream stream;
      if (std::numeric_limits<Target>::is_specialized)
	stream.precision(std::numeric_limits<Target>::digits10 + 1);
      else if(std::numeric_limits<Source>::is_specialized)
	stream.precision(std::numeric_limits<Source>::digits10 + 1);
      
      Target result;
      stream << arg;
      stream >> result;
      return result;
    };

  };

//______________________________________________________________________________
  template <typename Target, typename Source>
  inline
  Target 
  lexical_cast(Source arg)
  {
    return lexical_cast_impl<Target, Source>::lexical_cast(arg);
  };

//______________________________________________________________________________
  template <typename Target>
  inline
  Target
  b2i(const std::string& arg)
  {
    Target target = 0;
    int n_bits = arg.size()-1;
    for (int i=n_bits; i>0; --i)
      if ('1'==arg[i])
	target |= (1U<<(n_bits-i));
    return target;
  }

//______________________________________________________________________________
// template specialization : "string -> unsigned long long"
//______________________________________________________________________________
  template <>
  struct lexical_cast_impl<unsigned long long, std::string>
  {
    inline
    static unsigned long long
    lexical_cast(const std::string& arg)
    {
      unsigned long long result;
      std::stringstream stream;
      stream << std::setbase(0) << arg;
      stream >> result;
      return result;
    };

  };


//______________________________________________________________________________
// template specialization : "string -> unsigned int"
//______________________________________________________________________________
  template <>
  struct lexical_cast_impl<unsigned int, std::string>
  {
    inline
    static unsigned int 
    lexical_cast(const std::string& arg)
    {
      unsigned int result;
      std::stringstream stream;
      stream << std::setbase(0) << arg;
      stream >> result;
      return result;
    };

  };


//______________________________________________________________________________
// template specialization : "string -> int"
//______________________________________________________________________________
  template <>
  struct lexical_cast_impl<int, std::string>
  {
    inline
    static int 
    lexical_cast(const std::string& arg)
    {
      int result;
      std::stringstream stream;
      std::string::size_type minus_sign = arg.find("-");
      //      std::string::size_type white_space = arg.find_last_of(" ", minus_sign+1);
      if ((std::string::npos != minus_sign))
// 	  &&
// 	  (std::string::npos != white_space))
	{
	  //	  std::string tmp = arg.substr(white_space+1);
	  std::string tmp = arg.substr(minus_sign+1);
	  result 
	    = lexical_cast_impl<unsigned int, std::string>::lexical_cast(tmp) * (-1);
	}
      else
	result 
	  = lexical_cast_impl<unsigned int, std::string>::lexical_cast(arg);
      return result;
    };

  };

//______________________________________________________________________________
// template specialization : "string -> unsigned short"
//______________________________________________________________________________
  template <>
  struct lexical_cast_impl<unsigned short, std::string>
  {
    inline
    static unsigned short 
    lexical_cast(const std::string& arg)
    {
      unsigned short result;
      std::stringstream stream;
      stream << std::setbase(0) << arg;
      stream >> result;
      return result;
    };

  };

//______________________________________________________________________________
// template specialization : "string -> short"
//______________________________________________________________________________
  template <>
  struct lexical_cast_impl<short, std::string>
  {
    inline
    static short 
    lexical_cast(const std::string& arg)
    {
      short result;
      std::stringstream stream;
      std::string::size_type minus_sign = arg.find("-");
      //      std::string::size_type white_space = arg.find_last_of(" ", minus_sign+1);
      if ((std::string::npos != minus_sign))
// 	  &&
// 	  (std::string::npos != white_space))
	{
	  std::string tmp = arg.substr(minus_sign+1);
	  result
	    = lexical_cast_impl<unsigned short, std::string>::lexical_cast(tmp) * (-1);
	}
      else
	result 
	  = lexical_cast_impl<unsigned short, std::string>::lexical_cast(arg);

      return result;
    };

  };

//______________________________________________________________________________
// template specialization : "string -> bool"
//______________________________________________________________________________
  template <>
  struct lexical_cast_impl<bool, std::string>
  {
    inline
    static bool
    lexical_cast(const std::string& arg)
    {
      bool result = false;
      if (arg=="1" || arg=="true")
	result = true;
      else if (arg=="0" || arg=="false")
	result = false;
      return result;
    };
  };

//______________________________________________________________________________
// easy way
//______________________________________________________________________________

//______________________________________________________________________________
  inline
  unsigned long long
  a2ull(const std::string& arg)
  {
    return lexical_cast<unsigned long long>(arg);
  }

//______________________________________________________________________________
  inline
  unsigned int
  a2ui(const std::string& arg)
  {
    return lexical_cast<unsigned int>(arg);
  };
  
//______________________________________________________________________________
  inline
  int
  a2i(const std::string& arg)
  {
    return lexical_cast<int>(arg);
  };

//______________________________________________________________________________
  inline
  long long
  a2ll(const std::string& arg)
  {
    return lexical_cast<long long>(arg);
  };

//______________________________________________________________________________
  inline
  short
  a2s(const std::string& arg)
  {
    return lexical_cast<short>(arg);
  };

//______________________________________________________________________________
  inline
  unsigned short
  a2us(const std::string& arg)
  {
    return lexical_cast<unsigned short>(arg);
  };

//______________________________________________________________________________
  inline
  double
  a2d(const std::string& arg)
  {
    return lexical_cast<double>(arg);
  };


//______________________________________________________________________________
// int, ... -> string functions
//______________________________________________________________________________
  inline
  std::string
  ui2a(unsigned long long arg)
  {
    return lexical_cast<std::string>(arg);
  };

//______________________________________________________________________________
  inline
  std::string
  i2a(long long arg)
  {
    return lexical_cast<std::string>(arg);
  };

//______________________________________________________________________________
  inline
  std::string
  d2a(double arg)
  {
    return lexical_cast<std::string>(arg);
  };


  namespace string_operator
  {
    //  additional operators:
    //__________________________________________________________________________
    //  "string + int --> string", ....
    //__________________________________________________________________________
    inline
    std::string
    operator+(const std::string& lhs,
	      const unsigned long long rhs)
    {
      return lhs + lexical_cast<std::string>(rhs);
    }

    //__________________________________________________________________________
    inline
    std::string
    operator+(const std::string& lhs,
	      const long long rhs)
    {
      return lhs + lexical_cast<std::string>(rhs);
    }

    //__________________________________________________________________________
    inline
    std::string
    operator+(const std::string& lhs,
	      const double rhs)
    {
      return lhs + lexical_cast<std::string>(rhs);
    }

    //__________________________________________________________________________
    //  "int + string --> string", ....
    //__________________________________________________________________________
    inline
    std::string
    operator+(const unsigned long long lhs,
	      const std::string& rhs)
    {
      return lexical_cast<std::string>(lhs) + rhs;
    }

    //__________________________________________________________________________
    inline
    std::string
    operator+(const long long lhs,
	      const std::string& rhs)
    {
      return lexical_cast<std::string>(lhs) + rhs;
    }

    //__________________________________________________________________________
    inline
    std::string
    operator+(const double lhs,
	      const std::string& rhs)
    {
      return lexical_cast<std::string>(lhs) + rhs;
    }

    //__________________________________________________________________________
    //  "string += int", ....
    //__________________________________________________________________________
    template <typename T>
    inline
    void
    operator+=(std::string& lhs,
	       const unsigned long long rhs)
    {
      lhs += lexical_cast<std::string>(rhs);
      return;
    }

    //__________________________________________________________________________
    template <typename T>
    inline
    void
    operator+=(std::string& lhs,
	       const long long rhs)
    {
      lhs += lexical_cast<std::string>(rhs);
      return;
    }

    //__________________________________________________________________________
    template <typename T>
    inline
    void
    operator+=(std::string& lhs,
	       const double rhs)
    {
      lhs += lexical_cast<std::string>(rhs);
      return;
    }

  }

}
#endif
