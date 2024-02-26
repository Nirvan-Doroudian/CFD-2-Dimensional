#include <iostream>
#include <fstream>
#include <vector>

int main() {

    double L = 2.0;    // Length of board
    double W = 1.0;    // Width of board
    int Nx = 50;    // Number of grid points in the x direction
    int Ny = 25;    // Number of grid points in the y direction
    double dx = L / Nx;  // Delta-x
    double dy = W / Ny;  // Delta-y
    int nt = 500;           // Number of temporal steps
    double dt = 0.01;          // Time step size
    double Cp = 0.897;         // Cp (Aluminum)
    double Rho = 2700;         // Density (Aluminum)
    double C = Cp * Rho;        // For Convenience
    int a = 10;
    int iplus = 0;
    // Initialize Temp. array
    std::vector<std::vector<double>> Temperature(Nx, std::vector<double>(Ny, 20.0));

    // Apply boundary conditions
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            if (0 <= i && i <= 20 && j == Ny - 1) {
                Temperature[i][j] = 3 * i + 60;     // chip #1
            } else if (20 <= i && i <= 50 && j == Ny - 1) {
                Temperature[i][j] = -3 * i + 180;   // Top Edge B.C
            } else if (0 <= i && i <= 10 && j == 0) {
                Temperature[i][j] = 7 * i + 40;     // chip #2
            } else if (10 <= i && i <= 20 && j == 0) {
                Temperature[i][j] = -3 * i + 140;   // Bottom Edge B.C 1
            } else if (20 <= i && i <= 40 && j == 0) {
                iplus++;
                Temperature[i][j] = 2 * iplus + 80;     // chip #3
            } else if (40 <= i && i <= 50 && j == 0) {
                Temperature[i][j] = -9 * i + 480;   // Bottom Edge B.C 2
            } else if (0 <= j && j <= 25 && i == 0) {
                Temperature[i][j] = 0.4 * j + 40;   // Left Edge B.C
            } else if (0 <= j && j <= 5 && i == Nx - 1) {
                Temperature[i][j] = -j + 30;        // Right Edge B.C 1
            } else if (6 <= j && j <= 20 && i == Nx - 1) {
                Temperature[i][j] = 20;             // Cooling Fin
            } else if (21 <= j && j <= 25 && i == Nx - 1) {
                Temperature[i][j] = j - 20;        // Right Edge B.C 2
            }
        }
    }
    for (int k = 0; k < nt; k++) {
        std::vector<std::vector<double>> CTemperature = Temperature;

        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                double q = a * (1 + CTemperature[i][j] / 50);
                Temperature[i][j] = CTemperature[i][j] + (q * dt / C) * ((CTemperature[i+1][j] - 2 * CTemperature[i][j] + CTemperature[i-1][j]) / (dx * dx) + (CTemperature[i][j+1] - 2 * CTemperature[i][j] + CTemperature[i][j-1]) / (dy * dy));
            }
        }
    }

    // Export
    std::ofstream file("Profile.txt");
    file << "ZONE";
    for (int i = 0; i < Nx; i++) {
        for (int j = 12; j < 13; j++) {
            file << i << " " << Temperature[i][12] << "\n";
        }
    }
    file.close();

    return 0;
}