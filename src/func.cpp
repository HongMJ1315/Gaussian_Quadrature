#include "func.h"

#define FUNCSTR "(sin(4 * pi * x) + 1) * (cos(4 * pi * y) + 1)"
double F(double x, double y){
    return (sin(4 * PI * x) + 1) * (cos(4 * PI * y) + 1);
}