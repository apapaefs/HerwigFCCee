# HerwigFCCee

[ThePEG](https://thepeg.hepforge.org/) and [Herwig 7](https://herwig.hepforge.org/) development code with three additions for e+e- studies at LEP and FCC-ee:

- **Beam energy spread:** `ThePEG::GaussianBeamLuminosity` gives the two beams independent, truncated Gaussian energy spectra.
- **Decay cylinder:** with `DecayHandler:LimitCylinder On`, a particle decays only if its sampled laboratory decay point is inside a cylinder set by `CylinderRadius` and `CylinderHalfLength`.
- **Vertex smearing:** `ThePEG::GaussianVertexDistribution` smears the interaction point in position and time.

All three are off by default. The code also fixes two faults in `ShowerTree::clear()` that affect events whose shower is retried, and it adds a faster way to build Herwig. The readers of Les Houches event files, ThePEG's and the FxFx one, are also fixed for files that lack some of the processes their init block declares. They used to read past the end of their per-process arrays, which could give NaN or wrong cross sections.

This is not an official Herwig release.

## Contents

| Directory | Hepforge Mercurial repository | Base revision |
|---|---|---|
| `thepeghg/` | `thepeghg` | `2236:8912658c8679` |
| `herwighg/` | `herwighg` | `10868:040287d931af` |

The first commit is the unmodified upstream code, so comparing a later commit with it shows the local changes. The features and fixes follow as separate commits. The commits of the geometry series were not built one by one; only the complete set of changes has been validated.

`herwighg/` also includes two later upstream changesets, `10869:0d0396b157c1` and `10870:85e0e8fd839d`, which fix the LHE example cards in `herwighg/src`. They are separate commits, so a comparison with the first commit shows them too.

## Build

You need GNU Autotools with libtool, `gengetopt`, C++17 and Fortran compilers, GSL, Boost, LHAPDF, HepMC3 and FastJet. EvtGen and Pythia 8 are optional.

Build and install ThePEG first, then build Herwig against it. `configure` is not tracked, so run `autoreconf` in each tree first. If your libtool version differs from the one tracked in ThePEG, this also updates `thepeghg/Config/ltmain.sh` and the libtool macros in `thepeghg/m4`, and those files then show as modified.

In these commands, `PREFIX` is the install location, `DEPS` is the prefix that holds LHAPDF, FastJet, HepMC3, EvtGen and Pythia 8, and `BOOST_ROOT` and `GSL_ROOT` are the prefixes of Boost and GSL. `LDFLAGS` gives the plugins a run-time path to the libraries in `DEPS`.

On macOS with EvtGen, the Herwig install step adds `$PREFIX/lib` to the EvtGen plugin as another run-time path. The plugin needs `-headerpad_max_install_names` to leave room for it, and `DEPS` and `PREFIX` must be different directories, because `install_name_tool` refuses a duplicate path.

```bash
export LDFLAGS="-Wl,-rpath,$DEPS/lib"
# On macOS only:
export LDFLAGS="$LDFLAGS -Wl,-headerpad_max_install_names"
```

```bash
(cd thepeghg && autoreconf -fi)
mkdir -p build/thepeg && cd build/thepeg
../../thepeghg/configure --prefix="$PREFIX" \
  --with-lhapdf="$DEPS" --with-fastjet="$DEPS" --with-hepmc="$DEPS" --with-hepmcversion=3 \
  --with-boost="$BOOST_ROOT" --with-gsl="$GSL_ROOT" \
  --without-rivet --without-javagui --disable-readline
make -j8 && make -C src check-TESTS && make install
cd ../..
```

```bash
(cd herwighg && autoreconf -fi)
mkdir -p build/herwig && cd build/herwig
../../herwighg/configure --prefix="$PREFIX" --with-thepeg="$PREFIX" \
  --with-fastjet="$DEPS" --with-boost="$BOOST_ROOT" --with-gsl="$GSL_ROOT" \
  --with-evtgen="$DEPS" --with-pythia="$DEPS" \
  --without-openloops --without-gosam --without-gosam-contrib --without-njet \
  --without-vbfnlo --without-madgraph
make -j8 && make -C src check-TESTS && make install
cd ../..
```

EvtGen and Pythia are optional: use `--without-evtgen --without-pythia` to leave them out.

To compile Herwig faster, run `make -j8 parallel-objects` before `make -j8`. The recursive build compiles one directory at a time and leaves most cores idle. This step first makes the generated sources of every directory, then compiles the objects of all directories together, and `make` is left with little more than the links. On an 18-core Mac at `-j18`, this cut the Herwig build from 266 s to 129 s, with the same objects. Running more compilations at once also takes more memory. `make install` is faster too, because it installs the headers one directory at a time. Both changes are based on a patch by Juan Miguel.

The code has been tested on macOS arm64 with GCC 16 and Python 3.14, using these settings:

```bash
export CC=gcc-16 CXX=g++-16 FC=gfortran-16
export CFLAGS="-O2"
export CXXFLAGS="-O2 -std=c++17 -Wno-deprecated-declarations -Wno-deprecated-copy"
export FCFLAGS="-O2 -ffixed-line-length-none -fallow-argument-mismatch -fno-range-check"
export FFLAGS="$FCFLAGS"
export PYTHON=python3.14
```

This code uses ThePEG ABI 31:0:0. Rebuild anything compiled against an earlier ThePEG, and recreate saved runs.

## Example

```bash
Herwig read FCC-ee.in
Herwig run FCC-ee.run -N 1000
```

`FCC-ee.in` enables all three additions: independent 0.221% beam energy spreads, FCC vertex and time widths, and a cylinder of radius 2250 mm and half-length 2500 mm. It generates s sbar at 91.188 GeV with lepton ISR, QCD and QED final-state radiation and native decays.

`herwighg/README` describes the other example cards, and `thepeghg/README` describes the ThePEG changes.

## Licence

ThePEG and Herwig are distributed under version 3 of the GNU General Public License; see `COPYING` in each directory. This repository contains modified versions, changed in October and November 2026, and the commit history records the changes.
