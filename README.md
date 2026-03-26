# Synchrotron Polarization

## Background and Motivation

The synchrotron emission from relativistic electrons is intrinsically linearly polarized due to the anisotropic nature of their acceleration in magnetic fields. The frequency-dependent polarization degree is given by:

$$
\Pi(\omega) = \frac{\int_{0}^{\infty} \ G(x) N(\gamma) \ dx}{\int_{0}^{\infty} F(x)\ N(\gamma)\ dx}
$$

where:

* $N(\gamma)$: electron energy distribution
* $F(x), G(x)$: synchrotron kernel functions

Historically, the synchrotron polarization degree has been computed assuming a power-law electron distribution. In this limit, the polarization reduces to (see *Rybicki & Lightman, 1985*):

$$
\Pi = \frac{p + 1}{p + {7}/{3}}
$$

where (p) is the electron spectral index.

While widely used, this expression is only valid for idealized power-law distributions and uniform magnetic fields. In realistic astrophysical environments—such as blazar jets—the electron distribution is shaped by acceleration, cooling, and transport processes, and can deviate significantly from a pure power law.

This code instead computes the polarization using a **general electron energy distribution**, enabling a fully numerical and physically consistent treatment.

## Parameters related to turbulence

The following parameters (defined in `parameters.h`) are required:

```c
const double z       = 0.859;      // Source redshift
const double t_var   = 70.20e4;    // Variability timescale [s]
const double B       = 0.6;        // Magnetic field strength [G]
const double delta_D = 26.0;       // Doppler factor
const double D0      = 3.0e-6;     // Diffusion coefficient
```

These parameters control the emitting region, magnetic field properties, relativistic boosting, and turbulence (via diffusion), facilitate self-consistency between particle transport and polarization calculations.

## Method Summary

* Log-space integration of synchrotron functions $F(x), G(x)$
* Linear interpolation of electron distribution $N(\gamma)$
* Optional lookup tables for speed
* Asymptotic approximations for extreme regimes

## Compilation

```bash
gcc polarization.c -lgsl -lgslcblas -lm
```

**Libraries required:**

* GSL (GNU Scientific Library)
* libm (math library)


## Usage

```bash
./polarizarion <Electron-Distribution.dat> <Spectrum.dat> <FTotal.dat> <Output.dat> [LookupTable.dat]
```

### Input Files

* **Electron-Distribution.dat**
  Columns: gamma, N(gamma)

* **Spectrum.dat**
  Columns: nu, nuFnu (synchrotron spectrum)

* **FTotal.dat**
  Columns: nu, nuFnu (total SED)

* **LookupTable.dat (optional)**
  Precomputed values of synchrotron kernel functions

*(All files: two columns, no headers)*


### Output

* **Output.dat**
  Columns: nu, PD (polarization degree, 0–1)