// -*- C++ -*-
#include "GaussianVertexDistribution.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <cmath>
using namespace ThePEG;
GaussianVertexDistribution::GaussianVertexDistribution():
  meanX_(ZERO),meanY_(ZERO),meanZ_(ZERO),meanCT_(ZERO),
  sigmaX_(ZERO),sigmaY_(ZERO),sigmaZ_(ZERO),sigmaCT_(ZERO),
  spatialCutoff_(5.),timeCutoff_(5.) {}
void GaussianVertexDistribution::configure(const std::array<Length,4> & m,
    const std::array<Length,4> & s,double space,double time) {
  meanX_=m[0];meanY_=m[1];meanZ_=m[2];meanCT_=m[3];
  sigmaX_=s[0];sigmaY_=s[1];sigmaZ_=s[2];sigmaCT_=s[3];
  spatialCutoff_=space;timeCutoff_=time;validate();
}
void GaussianVertexDistribution::validate() const {
  for (Length x : {meanX_,meanY_,meanZ_,meanCT_})
    if (!std::isfinite(x/mm))
      throw Exception() << "Nonfinite Gaussian vertex mean." << Exception::runerror;
  for (Length x : {sigmaX_,sigmaY_,sigmaZ_,sigmaCT_})
    if (!std::isfinite(x/mm) || x<ZERO)
      throw Exception() << "Gaussian vertex widths must be finite and nonnegative." << Exception::runerror;
  if (!std::isfinite(spatialCutoff_) || !std::isfinite(timeCutoff_) || spatialCutoff_<=0 || timeCutoff_<=0)
    throw Exception() << "Gaussian vertex cutoffs must be finite and positive." << Exception::runerror;
}
void GaussianVertexDistribution::doinit() { Interfaced::doinit(); validate(); }
LorentzPoint GaussianVertexDistribution::sample() const {
  double x=0,y=0,z=0,t=0;
  do {
    if (sigmaX_>ZERO) x=UseRandom::rndGauss();
    if (sigmaY_>ZERO) y=UseRandom::rndGauss();
    if (sigmaZ_>ZERO) z=UseRandom::rndGauss();
  } while (x*x+y*y+z*z>spatialCutoff_*spatialCutoff_);
  if (sigmaCT_>ZERO) do { t=UseRandom::rndGauss(); } while (std::abs(t)>timeCutoff_);
  return LorentzPoint(meanX_+x*sigmaX_,meanY_+y*sigmaY_,meanZ_+z*sigmaZ_,meanCT_+t*sigmaCT_);
}
void GaussianVertexDistribution::persistentOutput(PersistentOStream & os) const {
  os << ounit(meanX_,mm) << ounit(meanY_,mm) << ounit(meanZ_,mm) << ounit(meanCT_,mm)
     << ounit(sigmaX_,mm) << ounit(sigmaY_,mm) << ounit(sigmaZ_,mm) << ounit(sigmaCT_,mm)
     << spatialCutoff_ << timeCutoff_;
}
void GaussianVertexDistribution::persistentInput(PersistentIStream & is,int) {
  is >> iunit(meanX_,mm) >> iunit(meanY_,mm) >> iunit(meanZ_,mm) >> iunit(meanCT_,mm)
     >> iunit(sigmaX_,mm) >> iunit(sigmaY_,mm) >> iunit(sigmaZ_,mm) >> iunit(sigmaCT_,mm)
     >> spatialCutoff_ >> timeCutoff_; validate();
}
DescribeClass<GaussianVertexDistribution,Interfaced> describeGaussianVertexDistribution
  ("ThePEG::GaussianVertexDistribution","libThePEG.so");
void GaussianVertexDistribution::Init() {
  static ClassDocumentation<GaussianVertexDistribution> docs
    ("Gaussian laboratory interaction point. Spatial coordinates share an ellipsoid cutoff; "
     "ct has a separate cutoff. All dimensional parameters are in mm. Zero widths draw no random numbers.");
#define VERTEX_MEAN(name, member) \
  static Parameter<GaussianVertexDistribution,Length> mean##name \
    ("Mean" #name,"Mean laboratory " #name " in mm (CT means c times time).", \
     &GaussianVertexDistribution::member,mm,ZERO,ZERO,ZERO,false,false,Interface::nolimits)
  VERTEX_MEAN(X,meanX_); VERTEX_MEAN(Y,meanY_); VERTEX_MEAN(Z,meanZ_); VERTEX_MEAN(CT,meanCT_);
#undef VERTEX_MEAN
#define VERTEX_SIGMA(name, member) \
  static Parameter<GaussianVertexDistribution,Length> sigma##name \
    ("Sigma" #name,"Gaussian " #name " width before truncation, in mm.", \
     &GaussianVertexDistribution::member,mm,ZERO,ZERO,ZERO,false,false,Interface::lowerlim)
  VERTEX_SIGMA(X,sigmaX_); VERTEX_SIGMA(Y,sigmaY_); VERTEX_SIGMA(Z,sigmaZ_); VERTEX_SIGMA(CT,sigmaCT_);
#undef VERTEX_SIGMA
  static Parameter<GaussianVertexDistribution,double> space
    ("MaxDeviationVertex","Maximum spatial distance in standardized coordinates.",
     &GaussianVertexDistribution::spatialCutoff_,5.,0.,0.,false,false,Interface::lowerlim);
  static Parameter<GaussianVertexDistribution,double> time
    ("MaxDeviationTime","Maximum absolute standardized ct displacement.",
     &GaussianVertexDistribution::timeCutoff_,5.,0.,0.,false,false,Interface::lowerlim);
}
