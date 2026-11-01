// -*- C++ -*-
//
// GaussianBeamLuminosity is part of ThePEG - Toolkit for HEP Event Generation.
// ThePEG is licensed under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_GaussianBeamLuminosity_H
#define ThePEG_GaussianBeamLuminosity_H

#include "LuminosityFunction.h"

namespace ThePEG {

/**
 * Independent Gaussian energy spectra for head-on, ultrarelativistic e-/e+
 * beams. Each spectrum is normalized within NSigma standard deviations of
 * its mean. RelativeSpreadA/B are sigma(E)/mean(E); a zero width leaves that
 * beam fixed and consumes no integration dimension.
 *
 * The beam energies neglect the electron mass. The event boost is computed
 * from the selected XComb's luminosity logarithms, including for unequal
 * beam means or widths. Maximum energies are derived from the spectrum.
 */
class GaussianBeamLuminosity : public LuminosityFunction {
public:
  /** Construct fixed 45.6 GeV beams; Gaussian widths default to zero. */
  GaussianBeamLuminosity();

  /** Accept an electron/positron pair in either beam order. */
  bool canHandle(const cPDPair &) const override;

  /** Number of beams with a nonzero width. */
  int nDim(const cPDPair &) const override;

  /** Upper and lower bounds on the collision energy before ISR. */
  Energy maximumCMEnergy() const override;
  Energy minimumCMEnergy() const override;

  /** Rapidity at the maximum beam energies. */
  double Y() const override;

  /** Laboratory rapidity at the supplied luminosity logarithms. */
  double beamRapidity(double lA, double lB) const override;

  /** Laboratory boost of the selected event. */
  LorentzRotation getBoost() const override;

  /** Map uniform variates to logarithmic energy losses and multiply the Jacobian. */
  pair<double,double> generateLL(const double *, double &) const override;

  /** Normalized luminosity density in the logarithmic energy variables. */
  double value(const cPDPair &, double, double) const override;

  /** Persistent output and input. */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /** Register the repository interfaces. */
  static void Init();

  /** Set and validate both means, relative widths and the Gaussian cutoff. */
  void configure(Energy a, Energy b, double da, double db, double cutoff);

  /** Upper endpoints of the individual beam spectra. */
  Energy maximumA() const { return meanA_*(1.0 + cutoff_*spreadA_); }
  Energy maximumB() const { return meanB_*(1.0 + cutoff_*spreadB_); }

protected:
  void doinit() override;
  IBPtr clone() const override { return new_ptr(*this); }
  IBPtr fullclone() const override { return new_ptr(*this); }

private:
  void validate() const;
  void syncMaxima();
  double sample(double u, double spread, double &jacobian) const;
  double density(double l, double spread) const;
  Energy meanA_, meanB_;
  double spreadA_, spreadB_, cutoff_;
};
}
#endif // ThePEG_GaussianBeamLuminosity_H
