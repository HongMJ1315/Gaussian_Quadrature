#include "tabulation.h"

double legendre(int n, double x) {
    if (n == 0)
        return 1.0;
    else if (n == 1)
        return x;
    else {
        double Pn_1 = x;       
        double Pn_2 = 1.0;    
        double Pn;
        for (int k = 2; k <= n; ++k) {
            Pn = ((2.0 * k - 1.0) * x * Pn_1 - (k - 1.0) * Pn_2) / k;
            Pn_2 = Pn_1;
            Pn_1 = Pn;
        }
        return Pn;
    }
}

double legendre_derivative(int n, double x) {
    return n / (x * x - 1.0) * (x * legendre(n, x) - legendre(n - 1, x));
}

std::vector<double> legendre_roots(int n, int max_iter, double tol) {
    std::vector<double> roots(n);
    int m = (n + 1) / 2;  

    for (int i = 0; i < m; ++i) {
        double x = std::cos(M_PI * (n - i - 0.25) / (n + 0.5));
        double x_prev;

        for (int iter = 0; iter < max_iter; ++iter) {
            double Pn = legendre(n, x);
            double Pn_prime = legendre_derivative(n, x);
            x_prev = x;
            x = x_prev - Pn / Pn_prime;

            if (std::abs(x - x_prev) < tol)
                break;
        }
        roots[i] = x;                      
        roots[n - i - 1] = -x;            
    }
    return roots;
}


std::vector<double> compute_weights(int n, const std::vector<double>& roots) {
    std::vector<double> weights(n);
    for (int i = 0; i < n; ++i) {
        double x = roots[i];
        double Pn_prime = legendre_derivative(n, x);
        weights[i] = 2.0 / ((1.0 - x * x) * Pn_prime * Pn_prime);
    }
    return weights;
}