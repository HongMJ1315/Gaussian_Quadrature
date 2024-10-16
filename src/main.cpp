#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "tabulation.h"

double F(double x, double y){
    return (x * x * x) / 9.0 + (y * y * y) / 10.0 + (y * y) / 5.0;
}

double gauss_quadrature_2D(int n, double (*F)(double, double), double integralMinX, double integralMaxX, double integralMinY, double integralMaxY){
    std::vector<double> roots = legendre_roots(n, 100, 1e-10);
    std::vector<double> weights = compute_weights(n, roots);

    double integral = 0.0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            double w = weights[i] * weights[j];
            double fx = F(roots[i] * (integralMaxX - integralMinX) / 2.0 + (integralMaxX + integralMinX) / 2.0, roots[j] * (integralMaxY - integralMinY) / 2.0 + (integralMaxY + integralMinY) / 2.0);
            integral += w * fx;            
        }
    }
    integral *= (integralMaxX - integralMinX) / 2.0 * (integralMaxY - integralMinY) / 2.0;
    std::cout << "(" << integralMinX << ", " << integralMaxX << "), (" << integralMinY << ", " << integralMaxY << "), Integral: " << integral << std::endl;

    return integral;
}

double gauss_quadrature_2D_grid(int n, double (*F)(double, double), double integralMinX, double integralMaxX, double integralMinY, double integralMaxY, double gridSize){
    double integral = 0.0;
    for(double x = integralMinX; x < integralMaxX; x += gridSize){
        for(double y = integralMinY; y < integralMaxY; y += gridSize){
            double gridIntegral = gauss_quadrature_2D(n, F, x, x + gridSize, y, y + gridSize);
            integral += gridIntegral;
            // std::cout << "X: (" << x << ", " << x + gridSize << "), Y: (" << y << ", " << y + gridSize << "), Integral: " << gridIntegral << std::endl; 
        }
    }
    std::cout << "Total Integral: " << integral << std::endl;
    return integral;
}

double integralMinX = -3.0;
double integralMaxX = 3.0;
double integralMinY = -3.0;
double integralMaxY = 3.0;

int main(){
    int n = 64;
    std::vector<double> roots = legendre_roots(n, 100, 1e-10);
    std::vector<double> weights = compute_weights(n, roots);

    double ans = gauss_quadrature_2D(n, F, integralMinX, integralMaxX, integralMinY, integralMaxY);
    gauss_quadrature_2D_grid(n, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 3);
    gauss_quadrature_2D_grid(n, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 2);
    gauss_quadrature_2D_grid(n, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 1.5);
    std::cout << ans << std::endl;

    return 0;
}
