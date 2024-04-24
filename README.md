# Power System Analysis Project - Load Analysis using Newton-Raphson Method

## Project Overview

This project aims to perform a load analysis of a power system with 4 buses using the Newton-Raphson iterative method. The objective is to calculate the voltage magnitudes, angles, power flows, and losses in the system. The analysis involves creating the admittance matrix, forming the mismatches vector, constructing the Jacobian matrix, updating the voltage values iteratively, calculating the slack bus power, and finally determining the load flow and losses of the system.

## Steps of Solution
1. Creation of the Admittance Matrix: Construct the admittance matrix based on the system topology and impedance data.

2. Formation of Mismatches Vector: Calculate the active power (P) and reactive power (Q) mismatches by comparing the calculated and scheduled power values. Calculate ΔP and ΔQ.

3. Formation of the Jacobian Matrix: Build the Jacobian matrix using the partial derivatives of the power flow equations with respect to voltage magnitudes and angles.

4. Solve the Jacobian Matrix: Solve the inverse of the Jacobian matrix using the mismatches vector to obtain the updated vector of voltage magnitudes and angles.

5. Calculate Slack Bus Power: Calculate the power injected at the slack bus based on the updated voltage values.

6. Repeat Iteratively: Repeat the above steps, updating the voltage values until a specific tolerance value is reached, indicating convergence.

7. Calculate Load Flow and Losses: Once convergence is achieved, calculate the load flow and losses in the power system.

## Repository Contents

The repository includes the following files:

1. **Calculate.m**: This file contains the calculation functions for determining the voltage magnitude, angle, tolerance, and updating the values. It also calculates the slack bus power based on the updated voltage values.

2. **Ybus.m**: This file implements the Ybus function, which calculates the admittance matrix of the power system.

3. **display_analysis.m**: This function displays the analysis outputs in a user-friendly manner, presenting the voltage magnitudes, angles, power flows, and system losses.

4. **load_flow.m**: This file contains the load_flow function, which performs the load flow calculations for the buses and determines the losses of the power system.

5. **main.m**: The main file for running the code correctly. It orchestrates the execution of the load analysis using the Newton-Raphson method. It calls the necessary functions and displays the results using display_analysis.m.

## Getting Started

To run the power system load analysis, follow these steps:

1. Place all the project files in the same directory.

2. Modify the input data and parameters in the main.m file, if required (but note that most of the code is hardcoded to solve a specific problem so beware of your changes).

4. Run the main.m file to execute the load analysis using the Newton-Raphson method.

5. Examine the results obtained, including voltage magnitudes, angles, power flows, and system losses, displayed in a user-friendly manner.
