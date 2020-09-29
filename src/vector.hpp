#ifndef VECTOR_HPP
#define VECTOR_HPP

/**
 * \file Vector.hpp
 * \author Matthieu Schaller
 * \brief Defines a 3 dimensionnal vector class
 * \date 16 December 2011
 *
 * Usual operators are overloaded and common mathematical functions are available
 *
 * \todo (Optimization) Implement meta-prog operators (Blitz++).
 * \todo (Optimization) Implement SSE operations.
 * \todo Write an optimized version of the inverseNorm function
 */


#include <cmath>
#include <sstream>


/**
 * \class Vector<>
 * \brief A 3D vector template class
 *
 * All the usual operators are overloaded. \n
 * Components are in the public part of the class.
 */
template<typename T>
struct Vector
{

  T x;  /**< X component*/
  T y;  /**< Y component*/
  T z;  /**< Z component*/

  /**
   *  \brief Builds a Vector with the 3 components specified
   */
  explicit Vector(const T& x_ = static_cast<T>(0),const T& y_ = static_cast<T>(0),const T& z_ = static_cast<T>(0))
    :x(x_),y(y_),z(z_) {}


  static const size_t dimension = 3; /**< Dimension of the vector*/

  template<typename S>
  Vector<T>& operator=(const Vector<S>& v)
  {
    x=v.x;
    y=v.y;
    z=v.z;
    return *this;
  }

  /**
   *  \brief Adds up two vectors
   */
  Vector<T>& operator+=(const Vector<T>& v)
  {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }

  /**
   *  \brief Adds up a vector with a number
   */
  Vector<T>& operator+=(const T& alpha)
  {
    x += alpha;
    y += alpha;
    z += alpha;
    return *this;
  }

  /**
   *  \brief Substracts two vectors
   */
  Vector<T>& operator-=(const Vector<T>& v)
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  /**
   *  \brief Substracts a vector by a number
   */
  Vector<T>& operator-=(const T& alpha)
  {
    x -= alpha;
    y -= alpha;
    z -= alpha;
    return *this;
  }

  /**
   *  \brief Multiplies a vector by a number
   */
  Vector<T>& operator*=(const T& alpha)
  {
    x *= alpha;
    y *= alpha;
    z *= alpha;
    return *this;
  }
  /**
   *  \brief Divides a vector by a number
   */
  Vector<T>& operator/=(const T& alpha)
  {
    x /= alpha;
    y /= alpha;
    z /= alpha;
    return *this;
  }

  /**
   * \brief Increments all components by 1 and returns the new vector
   */
  Vector<T>& operator++ ()
  {
    ++x;
    ++y;
    ++z;
    return *this;
  }

  /**
   *  \brief Increments all components by 1 and returns the old vector
   */
  Vector<T> operator++ (int)
  {
    Vector<T> v(*this);
    ++x;
    ++y;
    ++z;
    return v;
  }

  /**
   *  \brief Reduces all components by 1 and returns the new vector
   */
  Vector<T>& operator-- ()
  {
    --x;
    --y;
    --z;
    return *this;
  }

  /**
   *  \brief Reduces all components by 1 and return the old vector
   */
  Vector<T> operator-- (int)
  {
    Vector<T> v(*this);
    --x;
    --y;
    --z;
    return v;
  }

  /**
   * \brief Access to the specified component
   */
  T& operator[](size_t i)
  {
    return *(reinterpret_cast<T*>(this) + i);         // Dirty trick. Don't do this at home
  }

  /**
   * \brief Access to the specified component
   */
  const T& operator[](size_t i) const
  {
    return *(reinterpret_cast<const T*>(this) + i);         // Dirty trick. Don't do this at home
  }

};

/**
 *  \brief Compares two vectors component by component
 */
template<typename T>
inline bool operator==(const Vector<T>& v,const Vector<T>& w)
{
  return (v.x == w.x && v.y == w.y && v.z == w.z);
}

/**
 *  \brief Compares two vectors component by component
 */
template<typename T>
inline bool operator!=(const Vector<T>& v,const Vector<T>& w)
{
  return (v.x != w.x || v.y != w.y || v.z != w.z);
}

/**
 *  \brief Compares two vectors component by component
 *  Vectors are sorted according to their x component then to their y component and finally z
 *  This is useful for some operations implying STL containers and algorithms
 */
template<typename T>
inline bool operator<(const Vector<T>& v,const Vector<T>& w)
{
  return ((v.x < w.x) ||
	  (v.x == w.x && v.y < w.y) ||
	  (v.x == w.x && v.y == w.y && v.z < w.z));
}

/**
 *  \brief Compares two vectors component by component
 */
template<typename T>
inline bool operator>=(const Vector<T>& v,const Vector<T>& w)
{
  return !(v < w);
}

/**
 *  \brief Compares two vectors component by component
 */
template<typename T>
inline bool operator>(const Vector<T>& v,const Vector<T>& w)
{
  return ((v.x > w.x) ||
	  (v.x == w.x && v.y > w.y) ||
	  (v.x == w.x && v.y == w.y && v.z > w.z));
}

/**
 *  \brief Compares two vectors component by component
 */
template<typename T>
inline bool operator<=(const Vector<T>& v,const Vector<T>& w)
{
  return !(v > w);
}

/**
 *  \brief Returns the opposite of the Vector
 */
template<typename T>
inline Vector<T> operator-(const Vector<T>& v)
{
  return Vector<T>(-v.x, -v.y, -v.z);
}

/**
 *  \brief Applies the + operator to each of the Vector's component
 */
template<typename T>
inline Vector<T> operator+(const Vector<T>& v)
{
  return Vector<T>(+v.x, +v.y, +v.z);
}

/**
 *  \brief Computes the sum of two vectors
 */
template<typename T>
inline Vector<T> operator+(const Vector<T>& v,const Vector<T>& w)
{
  return Vector<T>(v) += w;
}

/**
 *  \brief Computes the difference between two vectors
 */
template<typename T>
inline Vector<T> operator-(const Vector<T>& v,const Vector<T>& w)
{
  return Vector<T>(v) -= w;
}

/**
 *  \brief Computes the product of a vector and a number
 */
template<typename T>
inline Vector<T> operator*(const Vector<T>& v,const T& alpha)
{
  return Vector<T>(v) *= alpha;
}

/**
 *  \brief Computes the product of a number and a vector
 */
template<typename T>
inline Vector<T> operator*(const T& alpha,const Vector<T>& v)
{
  return Vector<T>(v) *= alpha;
}

/**
 *  \brief Computes the ratio of a vector and a number
 */
template<typename T>
inline Vector<T> operator/(const Vector<T>& v,const T& alpha)
{
  return Vector<T>(v) /= alpha;
}

/**
 * \brief Computes the dot product of v and w
 */
template<typename T>
inline T operator*(const Vector<T>& v,const Vector<T>& w)
{
  return v.x*w.x + v.y*w.y + v.z*w.z;
}

/**
 * \brief Computes the dot product of v and w
 */
template<typename T>
inline T dotProduct(const Vector<T>& v,const Vector<T>& w)
{
  return v.x*w.x + v.y*w.y + v.z*w.z;
}

/**
 * \brief Computes the cross product of v and w
 */
template<typename T>
inline Vector<T> crossProduct(const Vector<T>& v ,const Vector<T>& w)
{
  return Vector<T>(v.y * w.z - v.z * w.y,
		   v.z * w.x - v.x * w.z,
		   v.x * w.y - v.y * w.x);
}

/**
 * \brief Computes the triple product of u,v and w
 */
template<typename T>
T tripleProduct(const Vector<T>& u,const Vector<T>& v ,const Vector<T>& w)
{
  return u.x*v.y*w.z - u.x*v.z*w.y + u.y*v.z*w.x - u.y*v.x*w.z + u.z*v.x*w.y - u.z*v.y*w.x;
}

/**
 * \brief Stream insertion operator
 *
 * Vectors are written as follow : (x,y,z)
 */
template<typename T, typename charT, class traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Vector<T>& v)
{
  std::basic_ostringstream<charT, traits> s;
  s.flags(os.flags());
  s.imbue(os.getloc());
  s.precision(os.precision());
  s << v.x << " " << v.y << " " << v.z;
  return os << s.str();
}


/**
 * \brief Stream extraction operator
 *
 * Vectors of the follwing form will be accepted:\n
 * (x,y,z) \n
 * (x) \n
 * x \n
 * The latter two producing a vector whose components are (x,x,x)
 * Any other input will raise the stream failbit.
 *
 */
template<typename T, typename charT, class traits>
std::basic_istream<charT, traits>&
operator>>(std::basic_istream<charT, traits>& is, Vector<T>& v)
{
  T x, y, z;
  charT temp;
  is >> temp;
  if (temp == '(')
    {
      is >> x >> temp;
      if (temp == ',')
        {
	  is >> y >> temp;
	  if (temp == ',')
            {
	      is >> z >> temp;
	      if(temp == ')')
		v = Vector<T>(x,y,z);
	      else
		is.setstate(std::ios_base::failbit);
            }
	  is.setstate(std::ios_base::failbit);
        }
      else if (temp == ')')
	v = Vector<T>(x,x,x);
      else
	is.setstate(std::ios_base::failbit);
    }
  else
    {
      is.putback(temp);
      is >> x;
      v = Vector<T>(x,x,x);
    }
  return is;
}

/**
 * \brief Computes the norm of the vector
 */
template <typename T>
T norm(const Vector<T>& v)
{
  return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

/**
 * \brief Computes the square norm of the vector
 */
template <typename T>
inline T norm2(const Vector<T>& v)
{
  return v.x*v.x + v.y*v.y + v.z*v.z;
}

/**
 * \brief Computes the Manhattan norm of the vector
 */
template <typename T>
inline T normManhattan(const Vector<T>& v)
{
  return std::abs(v.x) + std::abs(v.y) + std::abs(v.z);
}

/**
 * \brief Computes inverse of the norm of the vector
 */
template<typename T>
inline T inverseNorm(const Vector<T>& v)
{
  return 1. / norm(v);
}

/**
 * \brief Computes cosinus of the angle between two vectors
 */
template <typename T>
T cosAngle(const Vector<T>& v, const Vector<T>& w)
{
  return (v.x*w.x + v.y*w.y + v.z*w.z) / std::sqrt(norm2(v)*norm2(w));
}

/**
 * \brief Computes the angle between two vectors
 */
template <typename T>
T angle(const Vector<T>& v, const Vector<T>& w)
{
  return std::acos((v.x*w.x + v.y*w.y + v.z*w.z) / std::sqrt(norm2(v)*norm2(w)));
}

/**
 * \brief Computes the projection of the vector b on the vector a
 */
template<typename T>
Vector<T> proj(const Vector<T>& a,const Vector<T>& b)
{
  return Vector<T>(a)*=(a*b/norm2(a));
}

/**
 * \brief Returns a unitary vector in the direction of a
 */
template<typename T>
Vector<T> normalize(const Vector<T>& a)
{
  return Vector<T>(a)/norm(a);
}

/**
 * \brief Returns a vector whose components are the absolute value of the vector
 * given as an argument.
 */
template<typename T>
Vector<T> abs(const Vector<T>& a)
{
  return Vector<T>(std::abs(a.x), std::abs(a.y), std::abs(a.z));
}

/**
 * \brief Returns true if the vector is (0,0,0)
 */
template<typename T>
inline bool isNull(const Vector<T>& v)
{
  return (v.x == static_cast<T>(0) && v.y == static_cast<T>(0) && v.z == static_cast<T>(0));
}

/**
 * \brief Returns true if the vector is different from (0,0,0)
 */
template<typename T>
inline bool isNotNull(const Vector<T>& v)
{
  return (v.x != static_cast<T>(0) || v.y != static_cast<T>(0) || v.z != static_cast<T>(0));
}

/**
 * \brief Casts a vector of type T in a vector of type S. This is done by performing a static_cast on each component
 */
template<typename S,typename T>
Vector<S> vector_cast(const Vector<T>& v)
{
  return Vector<S>(static_cast<S>(v.x),static_cast<S>(v.y),static_cast<S>(v.z));
}



#endif  // VECTOR_HPP
