#ifndef TUPLE_HPP
#define TUPLE_HPP

#include <ostream>

template<size_t N, typename T>
struct Tuple
{
  const T& operator[](size_t i) const {return m_data[i];}
  T& operator[](size_t i) {return m_data[i];}

  template<typename S>
  void operator*=(const S& other)
  {
    for(size_t i(0); i<N; ++i)
      m_data[i] *= other;
  }

  size_t size() const
  {
    return N;
  }

  T m_data[N];

};

template<size_t N, typename T>
std::ostream& operator<<(std::ostream& flux, const Tuple<N, T>& tuple)
{
  for(int i(0); i<N-1; ++i)
    flux << tuple[i] << " ";
  flux << tuple[N-1];
  return flux;
}

template<size_t N, typename T>
inline T sum(const Tuple<N, T>& tuple)
{
  T total(0);
  for(size_t i(0); i<N; ++i)
    total += tuple[i];
  return total;
}

#endif //TUPLE_HPP
