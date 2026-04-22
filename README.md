# 31Cl β-delayed Proton-Detection Efficiency Simulation

Monte-Carlo simulation and analysis code for the proton-detection
efficiency of the **Proton Detector**, the gaseous component of the
**GADGET** (Gaseous Detector with Germanium Tagging) system at the
National Superconducting Cyclotron Laboratory (NSCL), Michigan State
University. Written in the context of a β-delayed proton-decay
measurement of <sup>31</sup>Cl. The project consists of three programs
that run sequentially:

1. **`RecoverBeamSpot`** fits the experimental beam-spot position and
   width from the measured central-pad and outer-quadrant-pad
   intensities.
2. **`GetParticleDistributions`** uses those beam-spot parameters (plus
   variations for systematic uncertainty) to sample the input
   distributions needed by the simulation.
3. **`SGMC`** runs the per-proton Monte-Carlo scan and produces the
   final efficiency summary.

## Physics overview

<sup>31</sup>Cl undergoes β<sup>+</sup> decay (T<sub>1/2</sub> = 190 ms)
to excited states of <sup>31</sup>S, which can then de-excite by
emitting a proton. Precise measurement of these β-delayed proton
branches provides an experimental constraint on the <sup>30</sup>P(p,
γ)<sup>31</sup>S reaction rate, which is important for understanding
nucleosynthesis in classical novae — particularly the abundance of
<sup>30</sup>Si observed in pre-solar stardust grains thought to have
condensed from nova ejecta.

The measurement is performed with GADGET, which pairs the Proton
Detector — a custom-built cylindrical gas-filled ionization chamber
with a Micromegas segmented-pad readout — with the Segmented Germanium
Array (SeGA) for prompt γ-ray detection. The <sup>31</sup>Cl beam is
implanted onto a thin Kapton window at the center of the Proton
Detector's active volume, where it decays in situ; the emitted protons
ionize P10 counting gas (90% Ar, 10% CH<sub>4</sub>) as they traverse
the chamber, and the liberated electrons drift to the Micromegas plane
where they are amplified and collected. A central circular pad measures
the signal of interest, four outer quadrant pads and an annular guard
ring provide veto and positional information, and SeGA tags coincident
γ-rays to pin down which excited state of <sup>31</sup>S emitted the
proton. For details of the apparatus see Friedman et al., Nucl.
Instrum. Methods A **940**, 93 (2019).

Detection efficiency depends on:

- **Beam-spot shape** on the Kapton window — a Gaussian profile whose
  centroid and 1/e<sup>2</sup> radius are fit from the measured pad
  intensities (see `RecoverBeamSpot`) and which drifts outward during
  the implantation half-lives.
- **Electron diffusion** during drift to the Micromegas, which blurs
  the apparent track and can push ionization charge outside the active
  region.
- **Stopping power** of the gas (dE/dx curve), which determines how
  many ionization electrons are produced along the track.
- **Veto threshold** on the outer guard ring, which discriminates
  decays that deposit too much charge outside the active volume.

This code evaluates the efficiency across a grid of variations in each
of these inputs so that the systematic uncertainty on the extracted
branching ratios can be read off directly from the spread in results.

## Repository layout

```
.
├── Makefile
├── README.md
├── .gitignore
├── RecoverBeamSpot.cpp           # stage-0: fit beam-spot x_0, y_0, R
├── GetParticleDistributions.cpp  # stage-1: sample input distributions
├── SGMC.cpp                      # stage-2: Monte-Carlo efficiency scan
├── RecoverBeamSpot.cfg           # stage-0 config
├── GetParticleDistributions.cfg  # stage-1 config
├── simulations.ctg               # stage-2 config
├── data/                         # read-only inputs
│   ├── tripleCoincidences.root   # reference beam-profile histogram
│   ├── padCounts.txt             # (optional) measured pad intensities
│   ├── protonStoppingPower808TorrP10gas.txt
│   ├── stoppingPower96percent.txt
│   └── stoppingPower104percent.txt
└── results/                      # generated outputs (gitignored)
    ├── ParticleDistributions.root
    └── efficiencySummary.txt
```

## Requirements

- **ROOT 6.24 or newer.** Earlier versions lack the
  `TF2::GetRandom2(x, y, TRandom*)` and `TH2::GetRandom2(x, y, TRandom*)`
  overloads that the parallel sampling relies on. Check with
  `root-config --version`.
- **A C++14 compiler** (`g++` 5.0+ or `clang++` 3.4+).
- **GNU make**. On macOS, the default `make` works; on BSD you may
  need `gmake`.
- **pthreads** (standard on Linux/macOS; the Makefile passes
  `-pthread`).

## Building

From the project root:

```bash
make                           # builds all three executables
make RecoverBeamSpot           # stage-0 only
make GetParticleDistributions  # stage-1 only
make SGMC                      # stage-2 only
make clean                     # removes binaries
```

Make sure ROOT's environment is set up first — if `root-config --libs`
prints nothing, source ROOT's setup script:

```bash
source /path/to/root/bin/thisroot.sh
```

## Running

The three stages run sequentially: stage 0 extracts beam-spot
parameters from data, stage 1 generates the input histograms, stage 2
runs the MC.

```bash
./RecoverBeamSpot          RecoverBeamSpot.cfg
./GetParticleDistributions GetParticleDistributions.cfg
./SGMC                     simulations.config
```

If no config path is passed, each program falls back to its default
filename (`RecoverBeamSpot.cfg`, `GetParticleDistributions.cfg`, and
`SGMC.config` respectively).

### Stage 0: `RecoverBeamSpot`

Fits a 2D Gaussian beam spot to the five measured pad intensities from
the Proton Detector — one central-pad count and four outer-quadrant
counts. The code integrates the assumed beam profile over each pad and
uses MINUIT (SIMPLEX → MIGRAD) to minimize the χ<sup>2</sup> between
calculated and measured pad-intensity ratios.

Outputs the best-fit `(x_0, y_0, R)` to stdout. Set `drawBeam = 1` in
the config to pop a TCanvas showing the fitted beam spot overlaid on
the pad geometry.

The extracted values feed into stage 1's `x_0`, `y_0`, and `R`
parameters (along with any ±1σ variations you want to sweep for
systematic uncertainty).

### Stage 1: `GetParticleDistributions`

Reads `data/tripleCoincidences.root`, samples the longitudinal
drift-time distribution, builds N Gaussian beam-spot histograms, and
builds M × 9 electron-cloud histograms (one per De value per drift time
from 0 to 8 µs). Writes everything to
`results/ParticleDistributions.root`.

Output histograms:

- `DriftTimes` — 1D reconstruction of the drift-time reference
  distribution
- `BeamSpot_k` for k = 0..N-1 — 2D beam spot for preset k
  (`numHalfLives`, `x_0`, `y_0`, `R` in the config are length-N parallel
  lists)
- `EC_De{k}_t{j}us` for k = 0..M-1, j = 0..8 — electron cloud family
  for De[k] at drift time j µs (j = 0 is internally treated as 0.1 µs
  to avoid a degenerate σ = 0)

### Stage 2: `SGMC`

Reads `results/ParticleDistributions.root` and runs the per-proton MC.
For each energy in `labEnergies`, for each combination of
`(beamSpotHist, stoppingFile, De-family)`, the code tracks `nProtons`
protons through the active volume and counts how many pass the
veto-threshold cut. The resulting efficiencies define the scan used to
extract the central value and systematic band.

Output: `results/efficiencySummary.txt` — human-readable table with
one block per energy, listing every variant's efficiency plus the
median/min/max across variants.

## Configuration

All three programs use the same simple format: `key = value` per line,
`#` starts a comment, blank lines are ignored. List-valued keys use
comma-separated values in stage 1 and whitespace-separated values in
stages 0 and 2 (historical accident; kept for backwards compatibility
with existing config files).

Required and optional keys for each stage are documented at the top of
the respective `.cpp` file, with defaults given in parentheses.

### Typical workflow

1. Edit `RecoverBeamSpot.cfg` to use your experiment's pad intensities
   (either set `padsExperimental` inline or point `intensityFile` at a
   text file with five numbers).
2. Run `./RecoverBeamSpot` and record the fitted `(x_0, y_0, R)`.
3. Edit `GetParticleDistributions.cfg` to use those values as the
   central preset in the `x_0`, `y_0`, `R` lists, and add ±1σ shifts
   as additional presets for systematic uncertainty.
4. Run `./GetParticleDistributions` to produce
   `ParticleDistributions.root`.
5. Edit `simulations.config` to use the list of beam-spot histograms
   and De values from the previous step, plus your stopping-power
   files, energies, and thresholds.
6. Run `./SGMC` to produce the efficiency summary.

### Tuning for your workload

The two knobs that matter most for runtime:

- **`nProtons` / `nElectrons`** (stage 1) and **`nProtons`** (stage 2)
  control per-histogram and per-variant statistics. Scale these
  together with the number of (stoppingFile × beamSpot × De × energy)
  variants to budget total CPU time.
- **`nThreads`** (stages 1 and 2) sets the worker-pool size. `0`
  auto-detects the core count. Stage 1 parallelizes over
  histogram-filling jobs (beam-spot jobs are ~100× more expensive than
  electron-cloud jobs, so speedup caps at ~N when you have N beam
  spots). Stage 2 parallelizes over proton energies (so speedup caps
  at the number of energies).

Stage 0 is not parallelized — it takes a few seconds at most, and
MINUIT has no natural axis to parallelize over.

## Parallelization notes

Stages 1 and 2 use a `std::thread` pool with atomic job-index dispatch.
ROOT's thread-safety caveats required several specific workarounds:

- **`TF2` formulas are pre-compiled on the main thread** and
  copy-constructed into each worker, because formula construction
  invokes the Cling JIT which has process-global state.
- **ROOT class dictionaries are warmed on the main thread** before
  workers are launched, so no worker triggers a lazy Cling dictionary
  load.
- **Input histograms are cloned serially on the main thread** into
  per-thread workspaces. `TH1::Clone` walks ROOT's global directory
  lists and is not thread-safe.
- **Each worker owns its own `TRandom3`** and all `GetRandom*` calls
  route through it explicitly, so `gRandom` is never touched.
- **Output file writes are serialized** on the main thread after all
  workers have joined. `TFile` is not thread-safe.

If you port to a newer ROOT that enables `ROOT::EnableThreadSafety()`
cleanly, some of the above can be simplified.

## Reproducibility

Each simulation job is seeded deterministically from
`std::time(nullptr) + constant * jobIndex`, so re-running with the
same clock time produces bit-identical output regardless of how jobs
migrate between threads. Different invocations differ because
`std::time` changes; hardcode a fixed seed in `seedBase` if you need
cross-run reproducibility.

## References

For the GADGET detector system:

> M. Friedman, D. Pérez-Loureiro, T. Budner, E. Pollacco, C. Wrede, et
> al., *GADGET: a Gaseous Detector with Germanium Tagging*, Nucl.
> Instrum. Methods A **940**, 93 (2019).
> [doi:10.1016/j.nima.2019.05.100](https://doi.org/10.1016/j.nima.2019.05.100)

For the GADGET II upgrade (a TPC-based successor to the present
Proton Detector):

> R. Mahajan, T. Wheeler, E. Pollacco, C. Wrede, et al., *Time
> Projection Chamber for GADGET II*, Phys. Rev. C **110**, 035807
> (2024).

## Citation

If you use this code in published work, please cite:

> *[Thesis / paper citation to be added upon publication.]*
