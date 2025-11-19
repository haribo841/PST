#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

// Function to calculate Erlang B using the summation definition (Eq 2 and 3)
// We use an iterative approach to calculate terms to avoid factorial overflow.
double calculateErlangB_Definition(double A, int C) {
    if (A < 0) return 0.0;

    double sum = 1.0;   // Represents the 0-th term (A^0 / 0!)
    double term = 1.0;  // Holds the value of the current term (A^i / i!)

    for (int i = 1; i <= C; ++i) {
        term = term * (A / i); // Calculate next term based on previous
        sum += term;           // Add to denominator sum
    }

    // Blocking probability E = P(C) = (A^C / C!) / Sum
    return term / sum;
}

double calculateErlangB_Recursive(double A, int C) {
    double E = 1.0; // E0
    for (int n = 1; n <= C; ++n)
        E = (A * E) / (n + A * E);
    return E; // EC
}

int main() {
    // Variables for input data
    double a_min, a_max, a_step;
    int C, m;

    // User Interface
    cout << "--- Erlang B Calculator (Exercise 2.1) ---" << endl;

    cout << "Enter system capacity (C): ";
    cin >> C;

    cout << "Enter min traffic per unit (a_min): ";
    cin >> a_min;

    cout << "Enter max traffic per unit (a_max): ";
    cin >> a_max;

    cout << "Enter step size (a_step): ";
    cin >> a_step;

    //'m' (number of classes), though Eq 4 uses single stream logic.
    // We ask for it to satisfy the input requirements.
    cout << "Enter number of traffic classes (m): ";
    cin >> m;

    // Prepare output file
    string filename = "wyniki_erlang.txt";
    ofstream outFile(filename);

    if (!outFile.is_open()) {
        cerr << "Error: Could not open file for writing!" << endl;
        return 1;
    }

    // Write Header as requested
    // Format: Capacity of system
    outFile << "Capacity (C): " << C << endl;
    outFile << "a_per_unit\tPb(E)" << endl;

    cout << "\nCalculating and writing to " << filename << "..." << endl;

    // Loop through traffic values
    // We use a small epsilon for float comparison safety in loops
    for (double a = a_min; a <= a_max + 0.00001; a += a_step) {

        // Calculate Total Offered Traffic A
        // A = a * C
        double A = a * C;

        // Calculate Blocking Probability
        //double E = calculateErlangB_Definition(A, C);
        double E = calculateErlangB_Recursive(A, C);
        // Write to file
        // 1st column: traffic per unit (a)
        // 2nd column: Blocking Probability (E)
        outFile << fixed << setprecision(6) << a << "\t" << E << endl;

        // Optional: Output to console to see progress
        cout << "a=" << a << ", A=" << A << " -> E=" << E << endl;
    }

    outFile.close();
    cout << "Calculations complete." << endl;

    return 0;
}