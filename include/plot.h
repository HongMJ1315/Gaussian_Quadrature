#ifndef PLOT_H
#define PLOT_H

#include <vector>
#include <iostream>
#include <fstream>
#include "func.h"

void generate_data(double x_min, double x_max, double y_min, double y_max, double x_step, double y_step, const std::string& filename);

#endif

