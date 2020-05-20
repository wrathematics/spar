#ifndef TIMER_HPP
#define TIMER_HPP


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
