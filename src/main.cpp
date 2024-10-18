#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include "tabulation.h"

double gauss_quadrature_2D(int n, double (*F)(double, double), double integralMinX, double integralMaxX, double integralMinY, double integralMaxY){
    std::vector<double> roots = legendre_roots(n, 100, 1e-10);
    std::vector<double> weights = compute_weights(n, roots);

    double integral = 0.0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            double w = weights[i] * weights[j];
            double fx = F(
                roots[i] * (integralMaxX - integralMinX) / 2.0 + (integralMaxX + integralMinX) / 2.0,
                roots[j] * (integralMaxY - integralMinY) / 2.0 + (integralMaxY + integralMinY) / 2.0
            );
            integral += w * fx;            
        }
    }
    integral *= (integralMaxX - integralMinX) / 2.0 * (integralMaxY - integralMinY) / 2.0;
    return integral;
}

double gauss_quadrature_2D_grid(int n, double (*F)(double, double), double integralMinX, double integralMaxX, double integralMinY, double integralMaxY, double gridSize){
    double integral = 0.0;
    for(double x = integralMinX; x < integralMaxX; x += gridSize){
        for(double y = integralMinY; y < integralMaxY; y += gridSize){
            integral += gauss_quadrature_2D(n, F, x, x + gridSize, y, y + gridSize);
        }
    }
    return integral;
}


double integralMinX = -3.0;
double integralMaxX = 3.0;
double integralMinY = -3.0;
double integralMaxY = 3.0;

std::vector<std::pair<double, double> > variableN, variableGridSize;
std::vector<std::pair<int, std::vector<std::pair<int, double > > > > totalIntegralResult;
int main(){
    // for(int i = 2; i <= 128; i *= 2){
    //     variableN.push_back(std::make_pair(i, gauss_quadrature_2D(i, F, integralMinX, integralMaxX, integralMinY, integralMaxY)));
    // }
    // for(int i = 1; i <= 128; i *= 2){
    //     variableGridSize.push_back(std::make_pair(i, gauss_quadrature_2D_grid(32, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 6.0 / i)));
    // }


    std::cout << gauss_quadrature_2D_grid(2, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 6.0 / 6) << std::endl;

    for(int i = 2; i <= 128; i *= 2){
        std::vector<std::pair<int, double> > row;
        for(int j = 1; j <= 128; j *= 2){
            row.push_back(std::make_pair(j, gauss_quadrature_2D_grid(i, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 6.0 / j)));
        }
        totalIntegralResult.push_back(std::make_pair(i, row));
    }

    // 將 (x, y, z) 座標數據寫入文件
    std::ofstream outFile("xy_z_coordinates.dat");
    for (const auto& pair : totalIntegralResult) {
        int x = pair.first;
        for (const auto& innerPair : pair.second) {
            int y = innerPair.first;
            double z = innerPair.second;
            outFile << x << " " << y << " " << z << std::endl;  // 輸出 (x, y, z)
        }
    }
    outFile.close();

 // 調用 gnuplot 繪製 (x, y, z) 柱狀圖
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        // 使用 wxt 終端彈出視窗
        fprintf(gnuplotPipe, "set terminal wxt size 800,600\n");
        fprintf(gnuplotPipe, "set title '3D Bar Chart (x, y, z)'\n");
        fprintf(gnuplotPipe, "set xlabel 'X (Sample Points)'\n");
        fprintf(gnuplotPipe, "set ylabel 'Y (Grid Size)'\n");
        fprintf(gnuplotPipe, "set zlabel 'Z (Integration Result)'\n");

        // 設置x和y軸為以2為底的對數刻度
        fprintf(gnuplotPipe, "set logscale x 2\n");
        fprintf(gnuplotPipe, "set logscale y 2\n");

        // 自定義x軸和y軸的刻度標籤，顯示2^n數據
        fprintf(gnuplotPipe, "set xtics ('2' 2, '4' 4, '8' 8, '16' 16, '32' 32, '64' 64, '128' 128, '256' 256)\n");
        fprintf(gnuplotPipe, "set ytics ('1' 1, '2' 2, '4' 4, '8' 8, '16' 16, '32' 32, '64' 64, '128' 128, '256' 256)\n");

        // 使用 boxes 繪製 3D 柱狀圖，並根據 z 值設置顏色
        fprintf(gnuplotPipe, "set style fill solid\n"); // 設置柱狀圖為實心
        fprintf(gnuplotPipe, "set palette defined (0 'blue', 1 'green', 2 'yellow', 3 'red')\n"); // 設置顏色範圍
        fprintf(gnuplotPipe, "set style data boxes\n"); // 使用 boxes 繪製 3D 柱狀圖
        fprintf(gnuplotPipe, "splot 'xy_z_coordinates.dat' using 1:2:3 with boxes lc palette notitle\n");

        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }


    /*
     // 將數據寫入文件
    std::ofstream outFileN("variableN.dat");
    for (const auto& pair : variableN) {
        outFileN << pair.first << " " << pair.second << std::endl;
    }
    outFileN.close();

    std::ofstream outFileGrid("variableGridSize.dat");
    for (const auto& pair : variableGridSize) {
        outFileGrid << pair.first << " " << pair.second << std::endl;
    }
    outFileGrid.close();

    // 調用 gnuplot 繪製 variableN 圖
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n");
        fprintf(gnuplotPipe, "set output 'variableN.png'\n");
        fprintf(gnuplotPipe, "set title 'Integration Result vs. N'\n");
        fprintf(gnuplotPipe, "set xlabel 'N (Number of Sample Points)'\n");
        fprintf(gnuplotPipe, "set ylabel 'Integration Result'\n");
        fprintf(gnuplotPipe, "set style data histograms\n");
        fprintf(gnuplotPipe, "set style histogram cluster gap 1\n");
        fprintf(gnuplotPipe, "set style fill solid border -1\n");
        fprintf(gnuplotPipe, "set boxwidth 0.9\n");
        fprintf(gnuplotPipe, "set yrange [0:*]\n"); // 設置Y軸範圍從0開始
        fprintf(gnuplotPipe, "plot 'variableN.dat' using 2:xtic(1) title 'Result'\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }

    // 調用 gnuplot 繪製 variableGridSize 圖
    gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n");
        fprintf(gnuplotPipe, "set output 'variableGridSize.png'\n");
        fprintf(gnuplotPipe, "set title 'Integration Result vs. Grid Size'\n");
        fprintf(gnuplotPipe, "set xlabel 'Grid Size (1/grid)'\n");
        fprintf(gnuplotPipe, "set ylabel 'Integration Result'\n");
        fprintf(gnuplotPipe, "set style data histograms\n");
        fprintf(gnuplotPipe, "set style histogram cluster gap 1\n");
        fprintf(gnuplotPipe, "set style fill solid border -1\n");
        fprintf(gnuplotPipe, "set boxwidth 0.9\n");
        fprintf(gnuplotPipe, "set yrange [0:*]\n"); // 設置Y軸範圍從0開始
        fprintf(gnuplotPipe, "plot 'variableGridSize.dat' using 2:xtic(1) title 'Result'\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }
    */


    return 0;
}
