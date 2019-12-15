// Copyright (C) 1996-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MATH_POWER_HPP
#define MATH_POWER_HPP

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>

namespace math {

using boost::is_arithmetic;
using boost::enable_if;
using boost::disable_if;

namespace detail {

template<typename T>
struct power_traits {
  typedef T power_type;
};
  
template<typename T>
struct power_traits<std::complex<T> > {
  typedef T power_type;
};

} // end namespace detail

//
// function power2 and p2
//

#ifndef BOOST_NO_SFINAE

template<typename T>
typename detail::power_traits<T>::power_type
power2(T t, typename enable_if<is_arithmetic<T> >::type* = 0) { return t * t; }

template<typename T>
typename detail::power_traits<T>::power_type
power2(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0) { return t * t; }

template<typename T>
typename detail::power_traits<T>::power_type
p2(T t, typename enable_if<is_arithmetic<T> >::type* = 0) { return power2(t); }

template<typename T>
typename detail::power_traits<T>::power_type
p2(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0) { return power2(t); }

#else
  
template<typename T>
typename detail::power_traits<T>::power_type
power2(T const& t) { return t * t; }

template<typename T>
typename detail::power_traits<T>::power_type
p2(T const& t) { return power2(t); }

#endif

template<typename T>
typename detail::power_traits<std::complex<T> >::power_type
power2(std::complex<T> const& t) { return power2(real(t)) + power2(imag(t)); }

//
// function power3 and p3
//

#ifndef BOOST_NO_SFINAE

template<typename T>
typename detail::power_traits<T>::power_type
power3(T t, typename enable_if<is_arithmetic<T> >::type* = 0) { return t * t * t; }

template<typename T>
typename detail::power_traits<T>::power_type
power3(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0) { return t * t * t; }

template<typename T>
typename detail::power_traits<T>::power_type
p3(T t, typename enable_if<is_arithmetic<T> >::type* = 0) { return power3(t); }

template<typename T>
typename detail::power_traits<T>::power_type
p3(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0) { return power3(t); }

#else
  
template<typename T>
typename detail::power_traits<T>::power_type
power3(T const& t) { return t * t * t; }

template<typename T>
typename detail::power_traits<T>::power_type
p3(T const& t) { return power3(t); }

#endif

//
// function power4 and p4
//

#ifndef BOOST_NO_SFINAE

template<typename T>
typename detail::power_traits<T>::power_type
power4(T t, typename enable_if<is_arithmetic<T> >::type* = 0) { return power2(power2(t)); }

template<typename T>
typename detail::power_traits<T>::power_type
power4(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0) { return power2(power2(t)); }

template<typename T>
typename detail::power_traits<T>::power_type
p4(T t, typename enable_if<is_arithmetic<T> >::type* = 0) { return power4(t); }

template<typename T>
typename detail::power_traits<T>::power_type
p4(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0) { return power4(t); }

#else
  
template<typename T>
typename detail::power_traits<T>::power_type
power4(T const& t) { return power2(power2(t)); }

template<typename T>
typename detail::power_traits<T>::power_type
p4(T const& t) { return power4(t); }

#endif

//
// function power6 and p6
//

#ifndef BOOST_NO_SFINAE

template<typename T>
typename detail::power_traits<T>::power_type
power6(T t, typename enable_if<is_arithmetic<T> >::type* = 0) { return power3(power2(t)); }

template<typename T>
typename detail::power_traits<T>::power_type
power6(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0) { return power3(power2(t)); }

template<typename T>
typename detail::power_traits<T>::power_type
p6(T t, typename enable_if<is_arithmetic<T> >::type* = 0) { return power6(t); }

template<typename T>
typename detail::power_traits<T>::power_type
p6(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0) { return power6(t); }

#else
  
template<typename T>
typename detail::power_traits<T>::power_type
power6(T const& t) { return power3(power2(t)); }

template<typename T>
typename detail::power_traits<T>::power_type
p6(T const& t) { return power6(t); }

#endif

} // end namespace math

#endif // MATH_POWER_HPP
