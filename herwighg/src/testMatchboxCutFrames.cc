#include "Herwig/MatrixElement/Matchbox/Cuts/IdentifiedParticleCut.h"
#include "Herwig/MatrixElement/Matchbox/Cuts/PairRapidityCut.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/Repository.h"
#include <iostream>
#include <stdexcept>

using namespace ThePEG;
using namespace Herwig;
namespace {
void require(bool ok, const char * message) {
  if (!ok) throw std::runtime_error(message);
}
void command(const string & text) {
  require(Repository::exec(text,std::cerr).find("Error")==string::npos,
          "could not configure cut fixture");
}
LorentzMomentum momentum(double y, double py) {
  const double mt=sqrt(25.+py*py);
  return LorentzMomentum(5*GeV,py*GeV,mt*sinh(y)*GeV,mt*cosh(y)*GeV);
}
}

int main() {
  try {
    auto cuts=new_ptr(Cuts());cuts->initialize(40000*GeV2,0.);
    auto single=new_ptr(IdentifiedParticleCut());auto pair=new_ptr(PairRapidityCut());
    auto type=ParticleData::Create(13,"frame-muon");
    command("create ThePEG::Matcher<Any> /FrameAny");
    Repository::Register(single,"/FrameSingle");Repository::Register(pair,"/FramePair");
    command("set /FrameSingle:Matcher /FrameAny");
    command("do /FrameSingle:YRange -1 1");
    command("set /FramePair:FirstMatcher /FrameAny");
    command("set /FramePair:SecondMatcher /FrameAny");
    command("do /FramePair:YRange -1 1");
    cuts->initSubProcess(10000*GeV2,0.);
    require(single->passCuts(cuts,type,momentum(.99,3)),"unboosted single edge");
    require(pair->passCuts(cuts,type,type,momentum(.99,3),momentum(.99,-3)),"unboosted pair edge");
    cuts->beamKinematics(40000*GeV2,.02);
    require(!single->passCuts(cuts,type,momentum(.99,3)),"beam boost must cross single rapidity edge");
    require(!pair->passCuts(cuts,type,type,momentum(.99,3),momentum(.99,-3)),"beam boost must cross pair rapidity edge");
    unsigned points=0;
    for(bool pseudo : {false,true}) {
      command(string("set /FramePair:Pseudo ")+(pseudo?"Yes":"No"));
      for(double beamY : {0.,.5*log(46.1/45.1),.5*log(30./920.),-.04,.04})
        for(double yhat : {-.3,0.,.2}) {
          cuts->beamKinematics(40000*GeV2,beamY);
          require(cuts->initSubProcess(10000*GeV2,yhat),"subprocess fixture");
          LorentzRotation boost(0.,0.,tanh(beamY+yhat));
          for(int i=0;i<401;++i) {
            auto p=momentum(-2.0037+i*.01001,3), q=momentum(-1.8273+i*.0091,-3);
            auto lab=p, pairLab=p+q;lab.transform(boost);pairLab.transform(boost);
            require(single->passCuts(cuts,type,p)==(abs(lab.rapidity())<=1.),
                    "single cut disagrees with explicit laboratory selection");
            require(pair->passCuts(cuts,type,type,p,q)==
                    (abs(pseudo?pairLab.eta():pairLab.rapidity())<=1.),
                    "pair cut disagrees with explicit laboratory selection");
            ++points;
          }
        }
    }
    cuts->initialize(40000*GeV2,.02);cuts->initSubProcess(10000*GeV2,0.);
    require(!single->passCuts(cuts,type,momentum(.99,3)),"fixed asymmetric single cut");
    std::cout << "PASS: single/pair rapidity edges, massive-pair pseudorapidity and "
              << points << " pointwise laboratory closure comparisons\n";
  } catch(const std::exception & error) {
    std::cerr << error.what() << '\n';return 1;
  }
}
