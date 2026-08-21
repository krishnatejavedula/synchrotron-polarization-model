# Synchrotron Polarization

## Background and Motivation

The synchrotron emission from relativistic electrons is intrinsically linearly polarized. The frequency-dependent polarization degree is given by:

$$
\Pi(\omega) = \frac{\int_{0}^{\infty} \ G(x) N(\gamma) \ dx}{\int_{0}^{\infty} F(x)\ N(\gamma)\ dx}
$$

where:

* $N(\gamma)$: electron energy distribution
* $F(x), G(x)$: synchrotron kernel functions

Frequently this expression is simplified by assuming a power law particle distribution. This code instead computes the polarization using a **general electron energy distribution**. This aproach is a fully numerical implementation that can physically connect polarization analysis to a complementary spectral analysis.

## Parameters related to turbulence

The following parameters (defined in `parameters.h`) are required:

```c
const double zz      = 0.859;      // Source redshift
const double t_var   = 70.20e4;    // Variability timescale [s]
const double B       = 0.6;        // Magnetic field strength [G]
const double deltaD  = 26.0;       // Doppler factor
const double D0      = 3.0e-6;     // Diffusion coefficient
const double ft      = 1;          // Turbulent scaling factor
```

These parameters can be set to values consistent with a spectral analysis and pertain to the emitting region, magnetic field properties, relativistic boosting, and turbulence, which can facilitate self-consistency between e.g. a particle transport model and this polarization calculation.

The library file `parameters.h` also defines the physical constants used in the calculation (e.g. speed of light, electron mass, electron charge, in CGS-Gaussian units). These are fixed and not meant to be changed.

## Method Summary

* Log-space integration of synchrotron functions $F(x), G(x)$
* Linear interpolation of input electron distribution $N(\gamma)$
* Optional lookup tables of the Bessel functions for computational efficiency
* Asymptotic approximations for extreme regimes of the Bessel functions (helps prevent errors and also for computational efficiency)

## Configuration Switches

Settings available in `polarization.h` (defaults shown):

* `polarization_mode` — `POL_BESSEL` / `POL_LOOKUP` (default) / `POL_ANALYTICAL`: how F(x), G(x) are evaluated.
* `lookup_mode` — `LOOKUP_BESSEL` / `LOOKUP_ANALYTICAL` (default): how the lookup table itself is generated.
* `scaling_mode` — `POL_UNSCALED` / `POL_SCALED` (default): master switch for the four factors below, applied once at the start of `main()`. `POL_UNSCALED` forces `ssa_mode`, `flux_weight_mode`, `field_disorder_mode`, and `turbulence_mode` all off, overriding whatever they're set to. `POL_SCALED` applies no override — the four switches take effect as set individually, so you can enable just one or two factors (e.g. SSA only) by leaving `scaling_mode = POL_SCALED` and setting the others to `_OFF`.
* `ssa_mode` — `SSA_ON` / `SSA_OFF` (default): exclude SSA weighting.
* `flux_weight_mode` — `FLUX_WEIGHT_OFF` / `FLUX_WEIGHT_ON` (default): include the synchrotron/total flux ratio.
* `field_disorder_mode` — `FIELD_DISORDER_OFF` / `FIELD_DISORDER_ON` (default): include the $(1- \eta)$ field disorder factor.
* `turbulence_mode` — `TURBULENCE_OFF` / `TURBULENCE_ON` (default): include the `ft` turbulence scaling factor.

The polarization degree value written to `Output.dat` is:

$$
\Pi_{\text{final}} = F_{SSA} \ F_{Syn} \ F_{\eta} \ F_T \ \Pi(\omega)
$$

All switches default to on (fully scaled — `scaling_mode = POL_SCALED`). Setting `scaling_mode = POL_UNSCALED` reduces `P_final` to the intrinsic Π(ω) by turning the other four switches off automatically without having to switch them off individually.

## Compilation in Command Line

The standard compilation anticipates the user has the GNU Scientific Library pre-installed. It is possible to compile the C code without GNU Tools, but you would need to comment out the GNU libraries in 'polarization.h' and use the analytical Bessel functions. 

```bash
gcc polarization.c -lgsl -lgslcblas -lm -o polarization
```

**Libraries required:**

* GSL (GNU Scientific Library)
* libm (math library)


## Usage

```bash
./polarization <Electron-Distribution.dat> <Spectrum.dat> <FTotal.dat> <Output.dat> <LookupTable.dat> <nu-tau.dat>
```

`LookupTable.dat` is required unless `polarization_mode` is not set to `POL_LOOKUP`. `nu-tau.dat` is required unless `ssa_mode` is set to `SSA_OFF`.

### Input Files

* **Electron-Distribution.dat**
  Columns: <gamma> <N(gamma)>

* **Spectrum.dat**
  Columns: <nu> <nuFnu (synchrotron spectrum)>

* **FTotal.dat**
  Columns: <nu> <nuFnu (total SED)>

* **LookupTable.dat**
  Precomputed values of synchrotron kernel functions - DO NOT EDIT

* **nu-tau.dat**
  Columns: <nu> <tau (SSA optical depth)>

*(All files: two space-separated columns, no headers)*


### Output

* **Output.dat**
  Columns: <nu> <PD (polarization degree, 0–1)>
