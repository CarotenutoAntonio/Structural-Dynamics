# Structural Dynamics with MATLAB

This repository contains MATLAB scripts developed for the Structural Dynamics course (A.A. 2024-2025) at the University of Naples “Federico II.” The projects focus on modal analysis, Frequency Response Functions (FRFs), and modal parameter extraction techniques.

## Projects / Main Scripts

### 1. Exercise 1: Modal Analysis - Analytical vs. FEM
   (Scripts: `autovalori_ROD.m`, `autovalori_bar.m`)

*   **Objective:** Determine and compare the natural frequencies of structures (rod and beam) using both analytical solutions and the Finite Element Method (FEM).
*   **Methodology (Rod - `autovalori_ROD.m`):**
    *   **Analytical Solution:** Calculates natural frequencies for a simply-supported rod using the formula $\omega_n = n (\pi/L) \sqrt{E/\rho}$.
    *   **FEM Implementation:**
        *   Discretizes the rod into N 1D elements (2 DOFs per element: longitudinal translation).
        *   Assembles global stiffness ($K$) and consistent mass ($M$) matrices.
        *   Applies boundary conditions for a simply-supported rod (ends constrained).
        *   Solves the generalized eigenvalue problem $K \Phi = \lambda M \Phi$ to find numerical natural frequencies.
        *   Compares FEM results with analytical solutions and evaluates percentage error as a function of the number of elements (DOFs).
*   **Methodology (Beam - `autovalori_bar.m`):**
    *   **Analytical Solution:** Calculates natural frequencies for a simply-supported beam using the formula $\omega_n = (n \pi/L)^2 \sqrt{EI/\rho A}$.
    *   **FEM Implementation:**
        *   Discretizes the beam into N 1D Euler-Bernoulli beam elements (4 DOFs per element: transverse displacement and rotation at each node).
        *   Assembles global stiffness ($K$) and consistent mass ($M$) matrices.
        *   Applies boundary conditions for a simply-supported beam (ends constrained for transverse displacement).
        *   Solves the generalized eigenvalue problem $K \Phi = \lambda M \Phi$.
        *   Compares FEM results with analytical solutions and evaluates percentage error.
*   **Key MATLAB Functions:** `linspace`, `eig`, `sort`, `diag`, plotting functions.

### 2. Exercise 2: Modal Characteristics and FRF
   (Script: `esercizio_3msse_6molle.m`)

*   **Objective:** Obtain modal characteristics (eigenvalues, eigenvectors) and Frequency Response Functions (FRFs) for a 3-mass, 6-spring system (as per Ewins) under different damping conditions.
*   **System Properties:** Uses defined mass ($m_1, m_2, m_3$) and stiffness ($k_1$ to $k_6$) values.
*   **Damping Conditions Studied:**
    1.  **Undamped System:**
    2.  **Proportional Hysteretic Damping:** Damping matrix $D_1 = \beta K$.
    3.  **Non-Proportional Hysteretic Damping:** Damping matrix $D_2$ with a non-zero term $D_2(1,1) = 0.3k$.
*   **Methodology:**
    *   Constructs global mass ($M$) and stiffness ($K$) matrices.
    *   For each damping case, solves the complex eigenvalue problem $(K + jD) \Phi = \lambda M \Phi$.
    *   Eigenvectors are mass-normalized.
    *   Calculates FRFs (receptance $\alpha_{jk}(\omega)$) using two methods:
        1.  Direct inversion: $[\alpha(\omega)] = (K + jD - \omega^2 M)^{-1}$.
        2.  Modal summation: $\alpha_{jk}(\omega) = \sum_{r=1}^{N} (\phi_j^{(r)} \phi_k^{(r)}) / (\lambda_r - \omega^2)$.
    *   Analyzes and plots eigenvalues, eigenvectors (magnitude and phase), and FRFs (direct and cross-receptance), discussing the influence of damping and frequency step on accuracy.
*   **Key MATLAB Functions:** `eig` (for complex systems), `sort`, `diag`, `inv`, `linspace`, plotting functions.

### 3. Exercise 3: Modal Parameter Extraction with Circle-Fit Technique
   (Script: `Circle_fit.m`)

*   **Objective:** Extract modal parameters (natural frequencies, damping ratios, modal constants) from a numerically generated FRF of a simply-supported beam using the Single-Degree-of-Freedom (SDOF) Circle-Fit technique.
*   **Methodology:**
    *   **FEM Model & FRF Generation:**
        *   Constructs an FEM model of a simply-supported beam (N=14 elements) with proportional hysteretic damping ($\beta = 0.05$).
        *   Determines modal properties (frequencies, damping, modes) and reconstructs a receptance function ($\alpha_{3,3}(\omega)$) to simulate experimental data.
    *   **Circle-Fit Procedure (applied to the "experimental" FRF):**
        1.  **Peak-Picking:** Identifies resonance peaks in the FRF magnitude (`findpeaks_basic` custom function).
        2.  **Frequency Band Selection:** Selects a $\pm10\%$ frequency window around each identified resonance peak.
        3.  **Circle Fitting:** Fits a circle to the real and imaginary parts of the FRF data within the selected window using a least squares method to find the circle's center ($x_c, y_c$) and radius ($R$).
        4.  **Resonance Frequency Estimation:** Estimates natural frequency using three methods: max imaginary part, zero real part, and max FRF magnitude (using `fminbnd`, `fzero`).
        5.  **Damping Ratio ($\eta_r$) Estimation:** Calculated from the frequencies corresponding to the horizontal diameter of the fitted circle.
        6.  **Modal Constant (${}_rA_{jk}$) Estimation:** Determined from the circle's radius, estimated natural frequency, and damping ratio.
    *   **Validation:** Compares the extracted modal parameters with the known "true" values from the FEM database.
*   **Key MATLAB Functions:** `eig`, `sort`, `diag`, `linspace`, `findpeaks_basic` (custom), `interp1` (for spline interpolation), `fminbnd`, `fzero`, plotting functions.

## Prerequisites

*   MATLAB (developed with a standard version, no specific toolboxes explicitly required by these scripts beyond standard MATLAB functionality).

## How to Run

For each `.m` script:
1.  Open the script in MATLAB.
2.  Ensure any required data files (e.g., `.mat` files from Exercise 2 if referenced, though these scripts appear self-contained) are in the MATLAB path or the same directory.
3.  Review and adjust parameters within the script if desired (e.g., number of elements, damping coefficients).
4.  Run the script.
5.  Figures and console outputs will be generated as described in the code and the associated report.

## Author

*   Antonio Carotenuto

## Acknowledgements

*   The course material and concepts from the Structural Dynamics course at the University of Naples “Federico II.”
*   Reference to "Modal Testing: Theory, Practice and Application" by D. J. Ewins for the 3-mass 6-spring system in Exercise 2.

