// -----------------------------------------------
// Synchrotron Polarization Calculation ----------
// -----------------------------------------------
/*
Synchrotron Polarization Calculator

Compilation:
  gcc code.c -lgsl -lgslcblas -lm

Required libraries:
  - GSL        : GNU Scientific Library
  - libm       : Standard math library

Usage (SSA off):
  ./code <Electron-Distribution.dat> <Spectrum.dat> <FTotal.dat> <Output.dat> [LookupTable.dat]

Usage (SSA on):
  ./code <Electron-Distribution.dat> <Spectrum.dat> <FTotal.dat> <Output.dat> [LookupTable.dat] <nu-tau.dat>

Required input files:
  1. Electron-Distribution.dat
     Electron energy distribution file: gamma, N(gamma)

  2. Spectrum.dat
     Synchrotron spectrum file: nu, nuFnu

  3. FTotal.dat
     Total flux spectrum / full blazar SED file: nu, nuFnu

Output file:
  4. Output.dat
     Calculated polarization degree file: nu, PD

Optional file:
  5. LookupTable.dat
     Precomputed lookup table for kernel values
*/
// ------------------ Libraries ------------------

#include "polarization.h"

// ------------------ Prototypes ------------------

double calculate_K(double nu);
double calculate_N_gamma(int n_intervals);
double calculate_deltaB();
double calculate_eta();
double syn_F(double K_gamma2);
double syn_G(double K_gamma2);
double interpolate(double gamma_query);
double lookup_F(double K_gamma2);
double lookup_G(double K_gamma2);
double integrand_G(double gamma, void *params);
double integrand_F(double gamma, void *params);
void load_data(const char *filename);
void load_lookup_table(const char *filename);
void load_spectrum(const char *filename);
double log_integrate(double (*func)(double, void *), double gamma_min, double gamma_max, void *params, int n_steps);
double calculate_polarization(double K_scale, double nu);
double interpolate_nu(double nu_query, double *nu_arr, double *flux_arr, int size);
void load_ftotal(const char *filename);
void load_tau(const char *filename);
double calculate_ssa_weight(double tau);
double calculate_x(void);

// ------------------ Main code ------------------

int main(int argc, char *argv[])
{
    if (polarization_mode == POL_BESSEL)
    {
        gsl_set_error_handler_off();
    }

    // Required argument check
    int required_argc = 5;
    if (polarization_mode == POL_LOOKUP)
    {
        required_argc++;
    }
    if (ssa_mode != SSA_OFF)
    {
        required_argc++;
    }

    if (argc < required_argc)
    {
        printf("Usage: %s <Electron-Distribution.dat> <Spectrum.dat> <FTotal.dat> <Output.dat> [LookupTable.dat] [nu-tau.dat]\n", argv[0]);
        return 1;
    }

    // Load electron distribution
    load_data(argv[1]);
    if (data_size == 0)
    {
        printf("Error: No data loaded from %s.\n", argv[1]);
        return 1;
    }

    // Load SSA spectrum
    load_spectrum(argv[2]);
    if (spectrum_size == 0)
    {
        printf("Error: No spectrum data loaded from %s.\n", argv[2]);
        return 1;
    }

    // Load final SED (total flux spectrum)
    load_ftotal(argv[3]);
    if (ftotal_size == 0)
    {
        printf("Error: No final SED data loaded from %s.\n", argv[3]);
        return 1;
    }

    // Open output file
    FILE *outfile = fopen(argv[4], "w");
    if (!outfile)
    {
        printf("Error: Cannot create output file %s.\n", argv[4]);
        return 1;
    }

    // Load lookup table if needed
    if (polarization_mode == POL_LOOKUP)
    {
        load_lookup_table(argv[5]);
    }

    if (ssa_mode != SSA_OFF)
    {
        int tau_arg = (polarization_mode == POL_LOOKUP) ? 6 : 5;
        load_tau(argv[tau_arg]);
        if (tau_size < 2)
        {
            printf("Error: No usable tau data loaded from %s.\n", argv[tau_arg]);
            fclose(outfile);
            return 1;
        }
    }

    double deltaB = calculate_deltaB();
    double eta = calculate_eta();
    double N_gamma = calculate_N_gamma(1000);

    // Main loop: compute and write polarization
    for (jj = 0; jj < spectrum_size; jj++)
    {
        double nu = nu_arr[jj];

        double f_ssa = flux_arr[jj];
        double f_total = interpolate_nu(nu, ftotal_nu, ftotal_arr, ftotal_size);

        double synF_ratio = 0.0;
        if (f_total > 0.0 && isfinite(f_total) && isfinite(f_ssa))
        {
            synF_ratio = f_ssa / f_total;
            if (!isfinite(synF_ratio) || synF_ratio < 0.0)
            {
                synF_ratio = 0.0;
            }
            if (synF_ratio > 1.0)
            {
                synF_ratio = 1.0;
            }
        }

        double K_scale = calculate_K(nu);
        double polarization = calculate_polarization(K_scale, nu);
        double eta = calculate_eta();
        double ssa_weight = 1.0;
        if (ssa_mode == SSA_ESCAPE_FACTOR)
        {
            if (nu >= tau_nu[0] && nu <= tau_nu[tau_size - 1])
            {
                double tau = interpolate_nu(nu, tau_nu, tau_arr, tau_size);
                ssa_weight = calculate_ssa_weight(tau);
            }
        }
        double P_final = polarization * ssa_weight * synF_ratio * (1 - eta) * ff;

        if (isfinite(P_final))
        {
            if (P_final < 0.0)
            {
                P_final = 0.0;
            }
            if (P_final > 1.0)
            {
                P_final = 1.0;
            }
            printf("nu = %.3e Hz | synF_ratio = %.3f | P_final = %.6f\n", nu, synF_ratio, P_final);
            fprintf(outfile, "%.6e %.6f\n", nu, P_final);
        }
    }

    // Summary
    printf("Integrated N_gamma = %e\n", N_gamma);
    printf("Calculated deltaB  = %e G\n", deltaB);
    printf("Calculated eta     = %e\n", eta);

    fclose(outfile);
    return 0;
}

// ------------------ Functions ------------------

double calculate_ssa_weight(double tau)
{
    if (!isfinite(tau) || tau <= 1.0e-10)
    {
        return 1.0;
    }

    double weight = (1.0 - exp(-tau)) / tau;
    if (!isfinite(weight) || weight < 0.0)
    {
        return 1.0;
    }
    if (weight > 1.0)
    {
        return 1.0;
    }

    return weight;
}

double calculate_K(double nu)
{
    return (4 * M_PI * nu * me * c) / (3 * e_charge * B);
}

double syn_F(double K_gamma2)
{
    if (K_gamma2 <= 0.0)
        return NAN;

    const int n_intervals = N_F_STEPS;
    const double upper = UPPER_LIMIT;

    const double uu_min = log(K_gamma2);
    const double uu_max = log(upper);
    const double duu = (uu_max - uu_min) / n_intervals;

    double sum = 0.0;
    for (ii = 0; ii < n_intervals; ++ii)
    {
        double uu = uu_min + (ii + 0.5) * duu;
        double yy = exp(uu);
        if (yy < 700.0)
        {
            double kval = gsl_sf_bessel_Knu(5.0 / 3.0, yy);

            if (kval >= 1e-300)
            {
                sum += kval * yy;
            }
        }
    }

    return K_gamma2 * duu * sum;
}

double syn_G(double K_gamma2)
{
    return K_gamma2 * gsl_sf_bessel_Knu(2.0 / 3.0, K_gamma2);
}

double lookup_F(double K_gamma2)
{
    if (K_gamma2 <= x_table[0])
        return F_table[0];
    if (K_gamma2 >= x_table[lookup_size - 1])
        return F_table[lookup_size - 1];

    for (ii = 0; ii < lookup_size - 1; ii++)
    {
        if (K_gamma2 >= x_table[ii] && K_gamma2 <= x_table[ii + 1])
        {
            double t = (K_gamma2 - x_table[ii]) / (x_table[ii + 1] - x_table[ii]);
            return F_table[ii] + t * (F_table[ii + 1] - F_table[ii]);
        }
    }
    return 0.0;
}

double lookup_G(double K_gamma2)
{
    if (K_gamma2 <= x_table[0])
        return G_table[0];
    if (K_gamma2 >= x_table[lookup_size - 1])
        return G_table[lookup_size - 1];

    for (ii = 0; ii < lookup_size - 1; ii++)
    {
        if (K_gamma2 >= x_table[ii] && K_gamma2 <= x_table[ii + 1])
        {
            double t = (K_gamma2 - x_table[ii]) / (x_table[ii + 1] - x_table[ii]);
            return G_table[ii] + t * (G_table[ii + 1] - G_table[ii]);
        }
    }
    return 0.0;
}

double analytical_F(double K_gamma2)
{
    if (K_gamma2 <= 0.0)
        return NAN;
    else if (K_gamma2 < 1.0e-4)
        return F1 * pow(K_gamma2, 1.0 / 3.0);
    else if (K_gamma2 > 1.0e4)
        return F2 * exp(-K_gamma2) * sqrt(K_gamma2);
    else
        return syn_F(K_gamma2);
}

double analytical_G(double K_gamma2)
{
    if (K_gamma2 <= 0.0)
        return NAN;
    else if (K_gamma2 < 1.0e-1)
        return G1 * pow(K_gamma2, 1.0 / 3.0);
    else if (K_gamma2 > 10.0)
        return G2 * exp(-K_gamma2) * sqrt(K_gamma2);
    else
        return K_gamma2 * gsl_sf_bessel_Knu(2.0 / 3.0, K_gamma2);
}

// ------------------ Data Loading ------------------

void load_data(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Error: Cannot open file %s\n", filename);
        data_size = 0;
        return;
    }
    ii = 0;
    while (ii < MAX_DATA && fscanf(file, "%lf %lf", &gamma_arr[ii], &elecdist_arr[ii]) == 2)
    {
        ii++;
    }
    data_size = ii;
    fclose(file);
}

void load_spectrum(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Error: Cannot open file %s\n", filename);
        spectrum_size = 0;
        return;
    }
    ii = 0;
    while (ii < MAX_SPECTRUM && fscanf(file, "%lf %lf", &nu_arr[ii], &flux_arr[ii]) == 2)
    {
        ii++;
    }
    spectrum_size = ii;
    fclose(file);
}

void load_lookup_table(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Error: Cannot open lookup table %s\n", filename);
        exit(1);
    }

    ii = 0;
    while (ii < MAX_LOOKUP && fscanf(file, "%lf %lf %lf", &x_table[ii], &F_table[ii], &G_table[ii]) == 3)
        ii++;

    lookup_size = ii;
    fclose(file);

    if (lookup_size < 2)
    {
        printf("Error: Lookup table too small or improperly formatted.\n");
        exit(1);
    }
}

double interpolate(double gamma_query)
{
    if (data_size < 2)
        return 0.0;

    for (ii = 0; ii < data_size - 1; ii++)
    {
        if (gamma_query >= gamma_arr[ii] && gamma_query <= gamma_arr[ii + 1])
        {
            double slope = (elecdist_arr[ii + 1] - elecdist_arr[ii]) / (gamma_arr[ii + 1] - gamma_arr[ii]);
            return elecdist_arr[ii] + slope * (gamma_query - gamma_arr[ii]);
        }
    }

    return (gamma_query < gamma_arr[0]) ? elecdist_arr[0] : elecdist_arr[data_size - 1];
}

double integrand_F(double gamma, void *params)
{
    double K_scale = *(double *)params;
    double K_gamma2 = K_scale / pow(gamma, 2.0);

    switch (polarization_mode)
    {
    case POL_BESSEL:
        return syn_F(K_gamma2) * interpolate(gamma);
    case POL_LOOKUP:
        return lookup_F(K_gamma2) * interpolate(gamma);
    case POL_ANALYTICAL:
        return analytical_F(K_gamma2) * interpolate(gamma);
    default:
        fprintf(stderr, "Error: Unknown polarization_mode (%d) in integrand_F\n", polarization_mode);
        return 0.0;
    }
}

double integrand_G(double gamma, void *params)
{
    double K_scale = *(double *)params;
    double K_gamma2 = K_scale / pow(gamma, 2.0);

    switch (polarization_mode)
    {
    case POL_BESSEL:
        return syn_G(K_gamma2) * interpolate(gamma);
    case POL_LOOKUP:
        return lookup_G(K_gamma2) * interpolate(gamma);
    case POL_ANALYTICAL:
        return analytical_G(K_gamma2) * interpolate(gamma);
    default:
        fprintf(stderr, "Error: Unknown polarization_mode (%d)", polarization_mode);
        return 0.0;
    }
}

double log_integrate(double (*func)(double, void *), double gamma_min, double gamma_max, void *params, int n_steps)
{
    if (gamma_min <= 0.0 || gamma_max <= 0.0)
    {
        fprintf(stderr, "Error: gamma_min and gamma_max must be positive for log integration.\n");
        return NAN;
    }

    double uu_min = log(gamma_min);
    double uu_max = log(gamma_max);
    double duu = (uu_max - uu_min) / n_steps;
    double sum = 0.0;

    for (nn = 0; nn < n_steps; nn++)
    {
        double uu = uu_min + (nn + 0.5) * duu;
        double gamma = exp(uu);
        double weight = gamma;
        sum += func(gamma, params) * weight;
    }

    return sum * duu;
}

double calculate_polarization(double K_scale, double nu)
{
    double gamma_min = gamma_arr[0];
    double gamma_max = gamma_arr[data_size - 1];
    int n_intervals = N_INT_STEPS;

    double pol_num = log_integrate(integrand_G, gamma_min, gamma_max, &K_scale, n_intervals);
    double pol_denom = log_integrate(integrand_F, gamma_min, gamma_max, &K_scale, n_intervals);

    return (pol_denom != 0.0) ? (pol_num / pol_denom) : 0.0;
}

double integrand_Ngamma(double gamma, void *params)
{
    // Linear interpolation using global arrays
    if (gamma <= gamma_arr[0])
        return elecdist_arr[0];
    if (gamma >= gamma_arr[data_size - 1])
        return elecdist_arr[data_size - 1];

    for (int i = 0; i < data_size - 1; i++)
    {
        if (gamma >= gamma_arr[i] && gamma < gamma_arr[i + 1])
        {
            double x0 = gamma_arr[i];
            double x1 = gamma_arr[i + 1];
            double y0 = elecdist_arr[i];
            double y1 = elecdist_arr[i + 1];
            return y0 + (gamma - x0) * (y1 - y0) / (x1 - x0);
        }
    }

    return 0.0; // fallback
}

double calculate_N_gamma(int n_intervals)
{
    if (data_size < 2)
    {
        fprintf(stderr, "Error: Not enough data points to integrate.\n");
        return 0.0;
    }

    double gamma_min = gamma_arr[0];
    double gamma_max = gamma_arr[data_size - 1];

    return log_integrate(integrand_Ngamma, gamma_min, gamma_max, NULL, n_intervals);
}

double calculate_deltaB()
{
    double N_gamma = calculate_N_gamma(1000); // Assumed comoving density

    double R_blob = (t_var * deltaD * c) / (1.0 + zz);
    double w0 = 2.0 * M_PI / R_blob;

    double numerator = 24.0 * D0 * N_gamma * c * me;
    double denominator = M_PI * pow(R_blob, 3.0) * w0;

    double delta_B = sqrt(numerator / denominator);
    return delta_B;
}

double calculate_eta()
{
    double deltaB = calculate_deltaB();
    double deltaB_squared = deltaB * deltaB;
    double B_squared = B * B;

    double eta = deltaB_squared / (B_squared + deltaB_squared);
    return eta;
}
void load_ftotal(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Error: Cannot open file %s\n", filename);
        ftotal_size = 0;
        return;
    }

    int i = 0;
    while (i < MAX_SPECTRUM && fscanf(file, "%lf %lf", &ftotal_nu[i], &ftotal_arr[i]) == 2)
    {
        i++;
    }

    ftotal_size = i;
    fclose(file);

    if (ftotal_size < 2)
    {
        printf("Error: Not enough data points loaded from %s.\n", filename);
    }
}
void load_tau(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Error: Cannot open file %s\n", filename);
        tau_size = 0;
        return;
    }

    int i = 0;
    while (i < MAX_SPECTRUM && fscanf(file, "%lf %lf", &tau_nu[i], &tau_arr[i]) == 2)
    {
        i++;
    }

    tau_size = i;
    fclose(file);

    if (tau_size < 2)
    {
        printf("Error: Not enough data points loaded from %s.\n", filename);
    }
}
double interpolate_nu(double nu_query, double *nu_arr, double *flux_arr, int size)
{
    if (size < 2)
        return 0.0;

    // Extrapolate if out of bounds
    if (nu_query <= nu_arr[0])
        return flux_arr[0];
    if (nu_query >= nu_arr[size - 1])
        return flux_arr[size - 1];

    for (int i = 0; i < size - 1; i++)
    {
        if (nu_query >= nu_arr[i] && nu_query <= nu_arr[i + 1])
        {
            double slope = (flux_arr[i + 1] - flux_arr[i]) / (nu_arr[i + 1] - nu_arr[i]);
            return flux_arr[i] + slope * (nu_query - nu_arr[i]);
        }
    }

    return 0.0; // fallback
}
