// -----------------------------------------------
// Synchrotron Polarization Header File ----------
// -----------------------------------------------

// ------------------ Libraries ------------------
#include <stdio.h>             // Standard I/O functions (printf, etc.)
#include <stdlib.h>            // Standard library (memory allocation, exit, etc.)
#include <math.h>              // Math functions (pow, sqrt, etc.)
#include <gsl/gsl_sf_bessel.h> // GSL Bessel special functions
#include <gsl/gsl_errno.h>     // GSL error handling
#include "parameters.h"         // User-defined constants

// ------------------ Numerical Limits ------------------
#define UPPER_LIMIT   1e3      // Integration upper limit
#define MAX_DATA      1001     // Max electron data points
#define MAX_SPECTRUM  1001     // Max spectrum data points
#define MAX_LOOKUP    10001    // Max lookup table size

// ------------------ Resolution Parameters ------------------
#define N_INT_STEPS   500      // Integration steps
#define N_F_STEPS     3000     // F-function steps
#define NUM_POINTS    1000     // Number of grid points
#define X_MIN_LOG     -4       // Minimum log10(x)
#define X_MAX_LOG      4       // Maximum log10(x)

// ------------------ Synchrotron Kernel Constants ------------------
#define F1 (M_PI * pow(2.0, 5.0 / 3.0) / (sqrt(3.0) * tgamma(1.0 / 3.0)))
#define F2 (sqrt(M_PI / 2.0))
#define G1 (F1 / 2.0)
#define G2 (F2)

// ------------------ Polarization Evaluation Modes ------------------
#define POL_BESSEL      0      // Use GSL Bessel functions
#define POL_LOOKUP      1      // Use precomputed lookup tables
#define POL_ANALYTICAL  2      // Use analytical expressions for F and G

int polarization_mode = POL_LOOKUP;   // Default mode

// ------------------ Synchrotron Self-Absorption Modes ------------------
#define SSA_OFF            0   // Do not apply SSA weighting
#define SSA_ESCAPE_FACTOR  1   // Apply (1 - exp(-tau)) / tau

// int ssa_mode = SSA_OFF;        // Default mode
int ssa_mode = SSA_ESCAPE_FACTOR; 
// #define POL_UNSCALED 0
// #define POL_SCALED   1

// int polarization_variant = POL_SCALED; //

#define LOOKUP_BESSEL 0
#define LOOKUP_ANALYTICAL 1
int lookup_mode = LOOKUP_ANALYTICAL;

// ------------------ Global Variables ------------------

double gamma_arr[MAX_DATA];    // Array to hold gamma values
double elecdist_arr[MAX_DATA]; // Array to hold electron distribution values
int data_size = 0;             // Size of the data loaded (i.e., number of valid points read from the file)

double nu_arr[MAX_SPECTRUM];   // Array to hold gamma values
double flux_arr[MAX_SPECTRUM]; // Array to hold electron distribution values
int spectrum_size = 0;

double ftotal_nu[MAX_SPECTRUM];  // Frequency grid for f_total
double ftotal_arr[MAX_SPECTRUM]; // Flux values
int ftotal_size = 0;

double tau_nu[MAX_SPECTRUM]; // Frequency grid for tau
double tau_arr[MAX_SPECTRUM]; // SSA optical depth values
int tau_size = 0;

double x_table[MAX_LOOKUP]; // Array to hold x values
double F_table[MAX_LOOKUP]; // Array to hold F values
double G_table[MAX_LOOKUP]; // Array to hold G values
int lookup_size = 0;
