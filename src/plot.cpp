#include "plot.h"  


void generate_data(double x_min, double x_max, double y_min, double y_max, double x_step, double y_step, const std::string& filename) {
    std::ofstream outfile(filename);
    for (double x = x_min; x <= x_max; x += x_step) {
        for (double y = y_min; y <= y_max; y += y_step) {
            double z = F(x, y);
            outfile << x << " " << y << " " << z << std::endl;
        }
        outfile << std::endl;  // 分隔不同的 x 值块
    }
    outfile.close();
}

