// This file is part of spar which is released under the Boost Software
// License, Version 1.0. See accompanying file LICENSE or copy at
// https://www.boost.org/LICENSE_1_0.txt

#ifndef SPAR_BENCHMARKS_TIMER_HPP
#define SPAR_BENCHMARKS_TIMER_HPP
#pragma once

#include <chrono>
#include <climits>


class timer
{
  public:
    timer();
    
    void start(const bool reset_first=false);
    void stop();
    void reset();
    double elapsed() const {return _elapsed;};
  
  private:
    std::chrono::high_resolution_clock::time_point query_clock() const;
    std::chrono::high_resolution_clock::time_point _start;
    double _elapsed;
};



inline timer::timer()
{
  reset();
}



inline void timer::start(const bool reset_first)
{
  if (reset_first)
    this->reset();
  
  _start = query_clock();
}


inline void timer::stop()
{
  std::chrono::duration<double> elapsed = query_clock() - _start;
  _elapsed += elapsed.count();
}



inline void timer::reset()
{
  _elapsed = 0.0;
}



inline std::chrono::high_resolution_clock::time_point timer::query_clock() const
{
  return std::chrono::high_resolution_clock::now();
}


#endif
