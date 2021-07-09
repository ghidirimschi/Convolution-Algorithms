#include <stdlib.h>
#include "timer.h"

#define SECONDS      1000000
#define MILLISECONDS    1000

void startTimer(Timer *timer) {
  gettimeofday(&(timer->start), NULL);
}

void stopTimer (Timer *timer)  {
  gettimeofday(&(timer->stop), NULL);
}

double seconds(Timer timer) {
  double tm = timer.stop.tv_sec*SECONDS + 
              timer.stop.tv_usec - timer.start.tv_sec*SECONDS -
              timer.start.tv_usec;
  return tm/SECONDS;
}

int milliseconds(Timer timer) {
  return (int)(seconds(timer)*1000);
}
