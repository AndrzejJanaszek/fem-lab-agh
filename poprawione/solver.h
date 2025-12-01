#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

// #############################
// wygenerowane chatemGPT
// #############################

// Rozwiązuje układ A x = b metodą Gaussa z częściowym pivotowaniem
std::vector<double> solveLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = (int)A.size();
    if (n == 0 || (int)A[0].size() != n || (int)b.size() != n) {
        throw std::runtime_error("Zle wymiary macierzy A lub wektora b");
    }

    // Eliminacja w przód
    for (int col = 0; col < n; ++col) {
        // 1. Znajdź największy element w kolumnie (pivot)
        int pivotRow = col;
        double pivotVal = std::fabs(A[col][col]);
        for (int row = col + 1; row < n; ++row) {
            if (std::fabs(A[row][col]) > pivotVal) {
                pivotVal = std::fabs(A[row][col]);
                pivotRow = row;
            }
        }

        if (std::fabs(pivotVal) < 1e-12) {
            throw std::runtime_error("Macierz osobliwa lub prawie osobliwa – brak jednoznacznego rozwiazania");
        }

        // 2. Zamień bieżący wiersz z wierszem pivotu
        if (pivotRow != col) {
            std::swap(A[pivotRow], A[col]);
            std::swap(b[pivotRow], b[col]);
        }

        // 3. Wyzeruj elementy poniżej pivota
        for (int row = col + 1; row < n; ++row) {
            double factor = A[row][col] / A[col][col];
            A[row][col] = 0.0; // już wiemy, że się wyzeruje

            for (int k = col + 1; k < n; ++k) {
                A[row][k] -= factor * A[col][k];
            }
            b[row] -= factor * b[col];
        }
    }

    // Podstawianie wsteczne
    std::vector<double> x(n);
    for (int row = n - 1; row >= 0; --row) {
        double sum = b[row];
        for (int col = row + 1; col < n; ++col) {
            sum -= A[row][col] * x[col];
        }
        x[row] = sum / A[row][row];
    }

    return x;
}
