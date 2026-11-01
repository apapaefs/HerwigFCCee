// The source XComb supplies both luminosity coordinates, including their order.
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/GaussianBeamLuminosity.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/PDF/PartonBin.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/RemnantHandler.h"
#include "ThePEG/PDT/ParticleData.h"
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace ThePEG;
namespace {
void require(bool condition, const char * message) {
  if (!condition) throw std::runtime_error(message);
}

// Supply the real luminosity function without a repository or event generator.
struct TestHandler : EventHandler {
  TestHandler(const cPDPair & beams, LumiFnPtr luminosity) : EventHandler(false) {
    theIncoming = beams;
    lumiFn(luminosity);
  }
};

// Only the StandardXComb constructor uses this scaffold; fill() must not
// evaluate an amplitude or generate phase space.
struct TestME : MEBase {
  IBPtr clone() const override { return new_ptr(*this); }
  unsigned int orderInAlphaS() const override { return 0; }
  unsigned int orderInAlphaEW() const override { return 2; }
  double me2() const override { throw std::runtime_error("unexpected amplitude call"); }
  Energy2 scale() const override { return GeV2; }
  int nDim() const override { return 0; }
  bool generateKinematics(const double *) override { return false; }
  CrossSection dSigHatDR() const override { return ZERO; }
  void getDiagrams() const override {}
  Selector<const ColourLines *> colourGeometries(tcDiagPtr) const override {
    return Selector<const ColourLines *>();
  }
};

void checkAccepted(Energy meanA, Energy meanB, double spread, DPair logs) {
  auto electron = ParticleData::Create(11,"contract-electron");
  auto positron = ParticleData::Create(-11,"contract-positron");
  auto photon = ParticleData::Create(22,"contract-photon");
  const cPDPair types(electron,positron);
  auto luminosity = new_ptr(GaussianBeamLuminosity());
  luminosity->configure(meanA,meanB,spread,spread,6);
  auto handler = new_ptr(TestHandler(types,luminosity));
  auto cuts = new_ptr(Cuts());
  const Energy emax = luminosity->maximumCMEnergy();
  cuts->initialize(sqr(emax),luminosity->Y());
  const PBPair bins(new_ptr(PartonBin()),new_ptr(PartonBin()));
  XComb source(emax,types,handler,tPExtrPtr(),tCascHdlPtr(),bins,cuts);
  source.lastP1P2(logs);
  const Energy a = luminosity->maximumA()*exp(-logs.first);
  const Energy b = luminosity->maximumB()*exp(-logs.second);
  const Energy2 s = 4*a*b;
  const Energy e = sqrt(a*b);
  source.lastS(s);
  const PPair beams(electron->produceParticle(LorentzMomentum(ZERO,ZERO,e,e)),
                    positron->produceParticle(LorentzMomentum(ZERO,ZERO,-e,e)));
  const PPair partons(electron->produceParticle(LorentzMomentum(ZERO,ZERO,.8*e,.8*e)),
                      positron->produceParticle(LorentzMomentum(ZERO,ZERO,-.6*e,.6*e)));
  const vector<Lorentz5Momentum> momenta{partons.first->momentum(),partons.second->momentum()};
  const DVector randoms{.123,.456};
  auto me = new_ptr(TestME());
  auto extractor = new_ptr(PartonExtractor());
  StandardXComb::DiagramVector diagrams;
  diagrams.push_back(new_ptr((Tree2toNDiagram(2),electron,positron,1,photon,3,electron,3,positron,-1)));
  StandardXComb target(emax,types,handler,tSubHdlPtr(),extractor,tCascHdlPtr(),
                      bins,cuts,me,diagrams,false);
  target.fill(beams,partons,momenta,randoms,source);
  require(target.lastP1()==logs.first && target.lastP2()==logs.second,
          "source coordinates were not copied in order");
  require(abs(target.lastS()/s-1)<1e-12 && abs(cuts->currentS()/s-1)<1e-12,
          "sampled collision energy was not restored");
  require(abs(cuts->currentBeamY()-.5*log(a/b))<1e-12,
          "beam rapidity was not restored from both source coordinates");
  require(abs(target.lastX1()-.8)<1e-12 && abs(target.lastX2()-.6)<1e-12 &&
          abs(target.lastSHat()/s-.48)<1e-12 && abs(target.lastY()-.5*log(.8/.6))<1e-12,
          "ISR momentum fractions or subprocess kinematics changed");
  require(target.lastParticles()==beams && target.lastPartons()==partons &&
          target.lastRandomNumbers()==randoms && target.meMomenta()==momenta,
          "fill did not copy the supplied phase-space state");
}
}

int main() {
  try {
    // Same invariant mass but opposite log ordering must give different Y.
    checkAccepted(60*GeV,35*GeV,.01,DPair(.01,.02));
    checkAccepted(60*GeV,35*GeV,.01,DPair(.02,.01));
    checkAccepted(45.6*GeV,45.6*GeV,0,DPair(0,0));
    checkAccepted(60*GeV,35*GeV,0,DPair(0,0));

    const PBPair bins(new_ptr(PartonBin()),new_ptr(PartonBin()));
    XComb source(100*GeV,cPDPair(),tEHPtr(),tPExtrPtr(),tCascHdlPtr(),bins,tCutsPtr());
    const DPair correct(.1,.2);
    source.lastP1P2(correct);
    const Energy2 s=10000*exp(-.3)*GeV2;
    source.lastS(s);
    auto data=ParticleData::Create(11,"contract-electron");
    const Energy e=sqrt(s)/2;
    const PPair beams(data->produceParticle(LorentzMomentum(ZERO,ZERO,e,e)),
                      data->produceParticle(LorentzMomentum(ZERO,ZERO,-e,e)));
    auto rejected=[&](const PPair & particles) {
      StandardXComb target;
      try { target.fill(particles,beams,vector<Lorentz5Momentum>(),DVector(),source); }
      catch(const Exception & error) {
        error.handle();
        require(!target.lastParticles().first && target.lastS()==ZERO,
                "invalid input changed the target before rejection");
        return true;
      }
      return false;
    };
    const double nan=std::numeric_limits<double>::quiet_NaN();
    for(const auto & logs : {DPair(nan,.2),DPair(.1,nan),DPair(.1,std::numeric_limits<double>::infinity())}) {
      source.lastP1P2(logs);
      require(rejected(beams),"nonfinite source coordinates accepted");
    }
    source.lastP1P2(correct);
    require(rejected(PPair()),"missing beams accepted");
    PPair wrong(beams);
    wrong.first=data->produceParticle(LorentzMomentum(ZERO,ZERO,1.01*e,1.01*e));
    require(rejected(wrong),"wrong collision mass accepted");
    source.lastS(1.01*s);
    require(rejected(beams),"inconsistent parent energy accepted");
    source.lastS(nan*GeV2);
    require(rejected(beams),"nonfinite parent energy accepted");
    std::cout << "PASS: parent-sourced luminosity, both log orderings, zero spread, "
                 "ISR fractions and cut context; invalid contexts rejected before mutation\n";
  } catch(const std::exception & error) {
    std::cerr << error.what() << '\n';return 1;
  }
}
