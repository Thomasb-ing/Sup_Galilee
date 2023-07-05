#include <time.h>
# if not defined(WIN32) && not defined(__USE_POSIX199309)
#   include <sys/time.h>
# endif
# include "Chronometer.hpp"

Chronometer::Chronometer() : m_time(0)
{}
// ------------------------------------------------------------------------
double
Chronometer::click()
{
# ifdef WIN32
  clock_t chrono;
  chrono = clock();
  double t = ((double)chrono)/CLOCKS_PER_SEC;
# elif defined(CLOCK_MONOTONIC_RAW)
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW , &tp);
  double t = tp.tv_sec+1.E-9*tp.tv_nsec;
# else
  struct timeval tv;
  gettimeofday(&tv,NULL);
  double t = tv.tv_sec+1.E-6*tv.tv_usec;
# endif
  double dt = t - m_time;
  m_time = t;
  return dt;
}
