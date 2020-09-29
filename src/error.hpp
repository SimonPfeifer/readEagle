#ifndef ERROR_HPP
#define ERROR_HPP

#include <iostream>
#include <cstdlib>
 
#define error(text)  do{std::cerr << "ERROR:" << __FILE__ << ':' << __func__ << "()" << __LINE__ << ": " << text << '\n';  std::abort();}while(0)

#define debug(text) do{std::cerr << "DEBUG:" << __FILE__ << ':' << __func__ << "()" << __LINE__ << ": " << text << '\n';}while(0)

#define debug_variable(text) do{std::cerr << "DEBUG:" << __FILE__ << ':' << __func__ << "()" << __LINE__ << ": " << #text << "=" << text << '\n';} while(0)

#endif //ERROR_HPP
