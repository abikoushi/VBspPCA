#include "lr.h"

double lr_default(const double & t,
                  const double & delay,
                  const double & forgetting){
  return pow(t+delay, -forgetting);
}