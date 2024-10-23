#ifndef TABULATION_H
#define TABULATION_H

/*-------------------------
tabulation.h
建立用於計算高斯積分的表格
-------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include "func.h"


double legendre(int n, double x);
double legendre_derivative(int n, double x);
std::vector<double> legendre_roots(int n, int max_iter, double tol);
std::vector<double> compute_weights(int n, const std::vector<double>& roots);




#endif