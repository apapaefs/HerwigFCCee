// -*- C++ -*-
//
// GaussianBeamLuminosity is part of ThePEG - Toolkit for HEP Event Generation.
// ThePEG is licensed under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "GaussianBeamLuminosity.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include <gsl/gsl_cdf.h>
#include <cmath>

using namespace ThePEG;

GaussianBeamLuminosity::GaussianBeamLuminosity()
  : meanA_(45.6*GeV), meanB_(45.6*GeV), spreadA_(0), spreadB_(0), cutoff_(6) {}

bool GaussianBeamLuminosity::canHandle(const cPDPair &p) const {
  return p.first && p.second && abs(p.first->id()) == 11 &&
    p.first->id() == -p.second->id();
}

void GaussianBeamLuminosity::validate() const {
  if (!std::isfinite(meanA_/GeV) || !std::isfinite(meanB_/GeV) ||
      !std::isfinite(spreadA_) || !std::isfinite(spreadB_) ||
      !std::isfinite(cutoff_) || meanA_ <= ZERO || meanB_ <= ZERO ||
      spreadA_ < 0 || spreadB_ < 0 || cutoff_ <= 0 || cutoff_ > 12 ||
      cutoff_*spreadA_ >= 1 || cutoff_*spreadB_ >= 1)
    throw InitException() << "GaussianBeamLuminosity requires finite positive means, "
      "nonnegative spreads, 0 < NSigma <= 12 and NSigma*spread < 1."
      << Exception::runerror;
}

void GaussianBeamLuminosity::syncMaxima() {
  validate();
  beamEMaxA(maximumA());
  beamEMaxB(maximumB());
}

void GaussianBeamLuminosity::configure(Energy a, Energy b, double da, double db, double k) {
  meanA_ = a; meanB_ = b; spreadA_ = da; spreadB_ = db; cutoff_ = k;
  syncMaxima();
}

void GaussianBeamLuminosity::doinit() {
  syncMaxima();
  LuminosityFunction::doinit();
}

int GaussianBeamLuminosity::nDim(const cPDPair &) const {
  return int(spreadA_ > 0) + int(spreadB_ > 0);
}
Energy GaussianBeamLuminosity::maximumCMEnergy() const {
  return 2*sqrt(maximumA()*maximumB());
}
Energy GaussianBeamLuminosity::minimumCMEnergy() const {
  return 2*sqrt(meanA_*meanB_*(1-cutoff_*spreadA_)*(1-cutoff_*spreadB_));
}
double GaussianBeamLuminosity::Y() const { return 0.5*log(maximumA()/maximumB()); }
double GaussianBeamLuminosity::beamRapidity(double a, double b) const {
  return Y() + 0.5*(b-a);
}
LorentzRotation GaussianBeamLuminosity::getBoost() const {
  return LorentzRotation(0.0, 0.0, tanh(beamRapidity(lastP1(), lastP2())));
}

double GaussianBeamLuminosity::sample(double u, double d, double &jac) const {
  if (!(u >= 0 && u <= 1))
    throw Exception() << "GaussianBeamLuminosity: integration variate outside [0,1]."
                      << Exception::runerror;
  const double q = gsl_cdf_ugaussian_P(-cutoff_);
  const double norm = 1-2*q;
  // Evaluate the small tail on both sides: never form 1 - a tiny probability.
  const double p = q + norm*std::min(u, 1-u);
  const double z = std::max(-cutoff_, std::min(cutoff_,
    (u < 0.5 ? 1 : -1)*gsl_cdf_ugaussian_Pinv(p)));
  const double phi = exp(-0.5*z*z)/sqrt(2*Constants::pi);
  jac *= d*norm/((1+d*z)*phi);
  return log1p(cutoff_*d)-log1p(d*z);
}

pair<double,double> GaussianBeamLuminosity::generateLL(const double *r, double &jac) const {
  double a = 0, b = 0;
  if (spreadA_ > 0) a = sample(*r++, spreadA_, jac);
  if (spreadB_ > 0) b = sample(*r, spreadB_, jac);
  return make_pair(a,b);
}

double GaussianBeamLuminosity::density(double l, double d) const {
  if (d == 0) return l == 0 ? 1 : 0;
  const double lmax = log1p(cutoff_*d)-log1p(-cutoff_*d);
  if (!std::isfinite(l) || l < 0 || l > lmax) return 0;
  const double de = expm1(log1p(cutoff_*d)-l);
  const double z = de/d;
  const double norm = erf(cutoff_/sqrt(2.0));
  return (1+de)*exp(-0.5*z*z)/(sqrt(2*Constants::pi)*d*norm);
}
double GaussianBeamLuminosity::value(const cPDPair &, double a, double b) const {
  return density(a, spreadA_)*density(b, spreadB_);
}

void GaussianBeamLuminosity::persistentOutput(PersistentOStream &os) const {
  os << ounit(meanA_,GeV) << ounit(meanB_,GeV) << spreadA_ << spreadB_ << cutoff_;
}
void GaussianBeamLuminosity::persistentInput(PersistentIStream &is, int) {
  is >> iunit(meanA_,GeV) >> iunit(meanB_,GeV) >> spreadA_ >> spreadB_ >> cutoff_;
  syncMaxima();
}

DescribeClass<GaussianBeamLuminosity,LuminosityFunction> describeGaussianBeamLuminosity
  ("ThePEG::GaussianBeamLuminosity", "");

void GaussianBeamLuminosity::Init() {
  static ClassDocumentation<GaussianBeamLuminosity> doc
    ("Independent truncated Gaussian accelerator energy spectra for head-on e-/e+ beams. "
     "RelativeSpread is sigma(E)/mean(E); the means and widths determine the beam maxima.");
  static Parameter<GaussianBeamLuminosity,Energy> a("MeanEnergyA", "Mean +z beam energy.",
    &GaussianBeamLuminosity::meanA_, GeV, 45.6*GeV, ZERO, ZERO, false, false, Interface::lowerlim);
  static Parameter<GaussianBeamLuminosity,Energy> b("MeanEnergyB", "Mean -z beam energy.",
    &GaussianBeamLuminosity::meanB_, GeV, 45.6*GeV, ZERO, ZERO, false, false, Interface::lowerlim);
  static Parameter<GaussianBeamLuminosity,double> da("RelativeSpreadA", "Relative standard deviation of beam A.",
    &GaussianBeamLuminosity::spreadA_, 0., 0., 1., false, false, Interface::limited);
  static Parameter<GaussianBeamLuminosity,double> db("RelativeSpreadB", "Relative standard deviation of beam B.",
    &GaussianBeamLuminosity::spreadB_, 0., 0., 1., false, false, Interface::limited);
  static Parameter<GaussianBeamLuminosity,double> k("NSigma", "Gaussian cutoff in standard deviations.",
    &GaussianBeamLuminosity::cutoff_, 6., 0., 12., false, false, Interface::limited);
  // Shadow the inherited writable maxima: these quantities are derived.
  static Parameter<GaussianBeamLuminosity,Energy> ma("BeamEMaxA", "Derived maximum beam A energy.",
    nullptr, GeV, 45.6*GeV, ZERO, ZERO, false, true, Interface::nolimits,
    nullptr, &GaussianBeamLuminosity::maximumA);
  static Parameter<GaussianBeamLuminosity,Energy> mb("BeamEMaxB", "Derived maximum beam B energy.",
    nullptr, GeV, 45.6*GeV, ZERO, ZERO, false, true, Interface::nolimits,
    nullptr, &GaussianBeamLuminosity::maximumB);
}
