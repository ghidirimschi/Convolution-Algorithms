#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>

struct Timer
{
  struct timeval start;
  struct timeval stop;
};

typedef struct Timer Timer;

void startTimer(Timer *timer);
void stopTimer(Timer *timer);
double seconds(Timer timer);
int milliseconds(Timer timer);

#endif
