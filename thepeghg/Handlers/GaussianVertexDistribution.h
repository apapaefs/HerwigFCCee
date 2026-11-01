// -*- C++ -*-
#ifndef ThePEG_GaussianVertexDistribution_H
#define ThePEG_GaussianVertexDistribution_H
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Vectors/LorentzVector.h"
#include <array>
namespace ThePEG {
/** Laboratory Gaussian coordinates with a spatial ellipsoid and separate ct cutoff. */
class GaussianVertexDistribution: public Interfaced {
public:
  GaussianVertexDistribution();
  LorentzPoint sample() const;
  void configure(const std::array<Length,4> &,const std::array<Length,4> &,
                 double spatialCutoff=5.,double timeCutoff=5.);
  void validate() const;
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &,int);
  static void Init();
protected:
  void doinit() override;
  IBPtr clone() const override { return new_ptr(*this); }
  IBPtr fullclone() const override { return new_ptr(*this); }
private:
  Length meanX_,meanY_,meanZ_,meanCT_,sigmaX_,sigmaY_,sigmaZ_,sigmaCT_;
  double spatialCutoff_,timeCutoff_;
};
}
#endif
