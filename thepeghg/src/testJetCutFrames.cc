#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Cuts/JetRegion.h"
#include "ThePEG/Cuts/JetPairRegion.h"
#include "ThePEG/Repository/Repository.h"
#include <iostream>
#include <stdexcept>

using namespace ThePEG;
namespace {
void require(bool ok, const char * message) {
  if (!ok) throw std::runtime_error(message);
}
void command(const string & text) {
  require(Repository::exec(text,std::cerr).find("Error")==string::npos,
          "could not configure cut fixture");
}
LorentzMomentum momentum(double y, double px=10.) {
  return LorentzMomentum(px*GeV,ZERO,10*sinh(y)*GeV,10*cosh(y)*GeV);
}
}

int main() {
  try {
    auto cuts=new_ptr(Cuts()); cuts->initialize(40000*GeV2,0.);
    auto jet=new_ptr(JetRegion()), first=new_ptr(JetRegion()), second=new_ptr(JetRegion());
    auto pair=new_ptr(JetPairRegion());
    Repository::Register(jet,"/FrameJet");
    Repository::Register(first,"/FrameFirst");
    Repository::Register(second,"/FrameSecond");
    Repository::Register(pair,"/FramePair");
    command("do /FrameJet:YRange -1 1");
    command("set /FramePair:FirstRegion /FrameFirst");
    command("set /FramePair:SecondRegion /FrameSecond");
    command("set /FramePair:OppositeHemispheres Yes");
    cuts->initSubProcess(10000*GeV2,0.);
    require(jet->matches(cuts,1,momentum(.99)),"unboosted edge must pass");
    jet->reset(); cuts->beamKinematics(40000*GeV2,.02);
    require(!jet->matches(cuts,1,momentum(.99)),"beam boost must cross jet rapidity edge");
    first->matches(cuts,1,momentum(-.01));
    second->matches(cuts,2,momentum(.1,-10));
    require(!pair->matches(cuts),"beam boost must change detector hemispheres");
    // Public no-parent mode already receives the complete explicit shift.
    jet->reset(); require(!jet->matches(tcCutsPtr(),1,momentum(.99),.02),"explicit no-parent shift");

    unsigned points=0;
    for(double beamY : {0.,.5*log(46.1/45.1),.5*log(30./920.),-.04,.04})
      for(double yhat : {-.3,0.,.2}) {
        cuts->beamKinematics(40000*GeV2,beamY);
        require(cuts->initSubProcess(10000*GeV2,yhat),"subprocess fixture");
        LorentzRotation boost(0.,0.,tanh(beamY+yhat));
        for(int i=0;i<401;++i) {
          auto p=momentum(-2.003+i*.01001), q=momentum(.207-i*.0021,-10);
          auto lab=p, labq=q;lab.transform(boost);labq.transform(boost);
          jet->reset();
          require(jet->matches(cuts,1,p)==(abs(lab.rapidity())<1.),
                  "JetRegion disagrees with explicit laboratory selection");
          first->reset();second->reset();
          first->matches(cuts,1,p);second->matches(cuts,2,q);
          require(pair->matches(cuts)==(lab.rapidity()*labq.rapidity()<0.),
                  "JetPairRegion disagrees with explicit laboratory hemispheres");
          require(abs((p.rapidity()-q.rapidity())-(lab.rapidity()-labq.rapidity()))<1e-12,
                  "rapidity separation must remain boost invariant");
          ++points;
        }
      }
    // Fixed asymmetric beams must also work without sampled beam context.
    cuts->initialize(40000*GeV2,.02);cuts->initSubProcess(10000*GeV2,0.);
    jet->reset();require(!jet->matches(cuts,1,momentum(.99)),"fixed asymmetric beam rapidity");
    std::cout << "PASS: jet rapidity/hemisphere edges and " << points
              << " pointwise laboratory closure comparisons\n";
  } catch(const std::exception & error) {
    std::cerr << error.what() << '\n';return 1;
  }
}
