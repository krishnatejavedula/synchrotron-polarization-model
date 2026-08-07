// -----------------------------------------------
// Synchrotron Polarization Constants File -------
// -----------------------------------------------

// ------------------ Physical Constants (CGS-Gaussian) ------------------
#define c_light  2.9979e10         // Speed of light [cm/s]
const double me       = 9.11e-28;  // Electron mass [g]
const double e_charge = 4.803e-10; // Electron charge [statC]

// Physical/model constants used in the synchrotron polarization calculation
const double zz      = 0.859;      // Source redshift
const double t_var   = 70.20e4;    // Variability timescale [s]
const double B       = 0.6;        // Magnetic field strength [G]
const double deltaD  = 26.0;       // Doppler factor
const double D0      = 3.0e-6;     // Diffusion coefficient
const double ft      = 1;          // Turbulent scaling factor