#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "tabulation.h"

#define FORMATE(x) std::scientific << std::setprecision(x)

struct Data{
    int n;
    double integral;
    time_t time;
    double error;
};

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
            // std::cerr << x << " " << x + gridSize << " " << y << " " << y + gridSize << std::endl;
            integral += gauss_quadrature_2D(n, F, x, std::min(3.0, x + gridSize), y, std::min(3.0,y + gridSize));
        }
    }
    return integral;
}

double dabs(double x){
    return x < 0 ? -x : x;
}


double integralMinX = -3.0;
double integralMaxX = 3.0;
double integralMinY = -3.0;
double integralMaxY = 3.0;

std::vector<std::pair<double, double> > variableN, variableGridSize;
std::vector<std::pair<int, std::vector<std::pair<int, double > > > > totalIntegralResult;
std::vector<Data> resultDiffP;
std::vector<Data> resultDiffH;
int main(){
    double currentIntegral = gauss_quadrature_2D_grid(15, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 6.0 / 16);
    std::cout << std::setprecision(15) << "Current Integral = " << currentIntegral << std::endl;

    for(int i = 2; i <= 4; i++){
        for(int j = 2; j <= 4; j++){
            std::cout << "N = " << i << ", Grid Size = " << j * j << ", Result = " << gauss_quadrature_2D_grid(i, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 6.0 / j) << std::endl;

        }
    }    

    /*-------------------------
    P-Refinement
    -------------------------*/
    std::ofstream pRef("P-Refinement.txt");
    for(int i = 2; i <= 33; i++){
        time_t start = time(NULL);
        double res = gauss_quadrature_2D(i, F, integralMinX, integralMaxX, integralMinY, integralMaxY);
        resultDiffP.push_back({i, res, time(NULL) - start, dabs((res - currentIntegral) / currentIntegral) * 100});
        // pRef << "N = " << i << ", Result = " << std::setprecision(10) << res << std::endl;
    }
    for(int i = 0; i < resultDiffP.size() / 2; i++){
        pRef << FORMATE(10) << resultDiffP[i].n << " & " << resultDiffP[i].integral  << " & " << resultDiffP[i].error << "\\% & " 
        << resultDiffP[resultDiffP.size() / 2 + i].n << " & " << resultDiffP[resultDiffP.size() / 2 + i].integral << " & "  << resultDiffP[resultDiffP.size() / 2+ i].error << "\\% \\\\" << std::endl; 
        pRef << "\\hline" << std::endl;
    }
    if(resultDiffP.size() % 2 == 1){
        pRef << FORMATE(10)  << resultDiffP[resultDiffP.size() / 2].n << " & " << resultDiffP[resultDiffP.size() / 2].integral << " & " << resultDiffP[resultDiffP.size() / 2].time << " & " << resultDiffP[resultDiffP.size() / 2].error << "\\% & " << " & " << " & " << " & " << " \\\\" << std::endl;
        pRef << "\\hline" << std::endl;
    }
    pRef.close();

    /*-------------------------
    H-Refinement
    -------------------------*/
    std::ofstream hRef("H-Refinement.txt");
    for(int i = 1; i <= 32; i++){
        time_t start = time(NULL);
        double res = gauss_quadrature_2D_grid(8, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 6.0 / i);
        resultDiffH.push_back({i, res, time(NULL) - start, dabs(dabs((res - currentIntegral) / currentIntegral) * 100)});
    }
    for(int i = 0; i < resultDiffH.size() / 2; i++){
        hRef << FORMATE(10)  << resultDiffH[i].n << " & " << resultDiffH[i].integral  << " & " << resultDiffH[i].error << "\\% & " 
        << resultDiffH[resultDiffH.size() / 2 + i].n << " & " << resultDiffH[resultDiffH.size() / 2 + i].integral << " & "  << resultDiffH[resultDiffH.size() / 2+ i].error << "\\% \\\\" << std::endl; 
        hRef << "\\hline" << std::endl;
    }
    if(resultDiffH.size() % 2 == 1){
        hRef << FORMATE(10) << resultDiffH[resultDiffH.size() / 2].n << " & " << resultDiffH[resultDiffH.size() / 2].integral << " & " << resultDiffH[resultDiffH.size() / 2].time << " & " << resultDiffH[resultDiffH.size() / 2].error << "\\% & " << " & " << " & " << " & " << " \\\\" << std::endl;
        hRef << "\\hline" << std::endl;
    }
    hRef.close();

    /*-------------------------
    HP-Refinement
    -------------------------*/
    for(int i = 2; i <= 11; i++){
        std::vector<std::pair<int, double> > row;
        for(int j = 3; j <= 10; j++){
            row.push_back(std::make_pair(j, dabs((gauss_quadrature_2D_grid(i, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 6.0 / j) - currentIntegral) / currentIntegral * 100)));
        }
        totalIntegralResult.push_back(std::make_pair(i, row));
    }
    std::ofstream hpRef("HP-Refinement.txt");
    for(int i = 0; i < totalIntegralResult[0].second.size(); i++)
        hpRef << " & " << totalIntegralResult[0].second[i].first;
    hpRef << "\\\\" << std::endl;
    hpRef << "\\hline" << std::endl;
    for(int i = 0; i < totalIntegralResult.size(); i++){
        hpRef << totalIntegralResult[i].first;
        for(int j = 0; j < totalIntegralResult[i].second.size(); j++){
            hpRef << " & " << std::fixed <<  std::setprecision(2) << totalIntegralResult[i].second[j].second << "\\%";
        }
        hpRef << "\\\\" << std::endl;
        hpRef << "\\hline" << std::endl;
    }
    totalIntegralResult.clear();
    for(int i = 2; i <= 20; i++){
        std::vector<std::pair<int, double> > row;
        for(int j = 1; j <= 20; j++){
            row.push_back(std::make_pair(j, dabs((gauss_quadrature_2D_grid(i, F, integralMinX, integralMaxX, integralMinY, integralMaxY, 6.0 / j) - currentIntegral) / currentIntegral * 100)));
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
            outFile << FORMATE(20) << x << " " << y << " " << z << std::endl;  // 輸出 (x, y, z)
        }
    }
    outFile.close();

 // 調用 gnuplot 繪製 (x, y, z) 柱狀圖
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set terminal wxt size 800,600\n");
        // 设置绘图参数
        fprintf(gnuplotPipe, "set xlabel 'Sample Point'\n");
        fprintf(gnuplotPipe, "set ylabel 'Grid Size'\n");
        fprintf(gnuplotPipe, "set zlabel 'Error Value'\n");
        fprintf(gnuplotPipe, "set title 'Error Value 3D Bar Chart'\n");
        // 设置颜色映射
        fprintf(gnuplotPipe, "set palette defined (0 'blue', 1 'green', 2 'yellow', 3 'red')\n");
        // 设置视角
        fprintf(gnuplotPipe, "set view 60, 30, 1, 1\n");
        // 绘图
        fprintf(gnuplotPipe, "splot 'xy_z_coordinates.dat' using 1:2:3 with impulses lc palette lw 5 notitle\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }

    // 使用 gnuplot 繪製 F(x, y) 函數圖(使用FUNCSTR宏定義的函數) 將結果存儲在文件中
    gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set terminal qt size 800,600\n");
        fprintf(gnuplotPipe, "set title 'F(x, y) = %s'\n", FUNCSTR);
        fprintf(gnuplotPipe, "set xlabel 'X'\n");
        fprintf(gnuplotPipe, "set ylabel 'Y'\n");
        fprintf(gnuplotPipe, "set pm3d\n");
        fprintf(gnuplotPipe, "set hidden3d\n");
        fprintf(gnuplotPipe, "set isosamples 50\n");
        fprintf(gnuplotPipe, "splot [-3:3][-3:3] %s\n", FUNCSTR);
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    }

    return 0;
}
