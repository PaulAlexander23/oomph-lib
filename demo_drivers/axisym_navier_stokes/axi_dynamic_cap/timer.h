#ifndef TIMER_H
#define TIMER_H

#include <chrono>

namespace oomph
{
  class Timer
  {
  private:
    /// Start time
    /// The time when the timer was started
    std::chrono::high_resolution_clock::time_point Start_time;

  public:
    /// Constructor
    /// Set the start time to the current time
    /// when the object is created
    Timer()
    {
      Start_time = std::chrono::high_resolution_clock::now();
    }

    /// Destructor
    ~Timer() {}

    /// The time in seconds since the timer was created can be retrieved by
    /// calling the time_elapsed() function
    std::chrono::duration_cast<std::chrono::seconds> time_elapsed()
    {
      return std::chrono::high_resolution_clock::now() - Start_time;
    }
  };
}; // namespace oomph
#endif
