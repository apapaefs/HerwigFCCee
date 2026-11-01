// Deterministic checks for the Gaussian luminosity and sampled beam context.
#include "ThePEG/Handlers/GaussianBeamLuminosity.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/PDF/PDFCuts.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include <sstream>
#include <gsl/gsl_integration.h>
#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace ThePEG;
void check(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
int main() {
  try {
    GaussianBeamLuminosity g;
    const cPDPair types;
    double maxError = 0;
    for (auto d : std::vector<std::pair<double,double>>{{0,0},{.00221,0},{0,.00221},{.00221,.00221},{.001,.003},{.01,.01}}) {
      for (double k : {6.,7.}) {
        g.configure(46.1*GeV,45.1*GeV,d.first,d.second,k);
        check(g.nDim(types)==int(d.first>0)+int(d.second>0), "sampling dimension");
        for (double u : {0.,1e-12,.001,.1,.5,.9,.999,1-1e-12,1.}) {
          const double r[] = {u,1-u};
          double jac=1;
          const auto l=g.generateLL(r,jac);
          const double err=std::abs(jac*g.value(types,l.first,l.second)-1);
          maxError=std::max(maxError,err);
          if(err>=1e-10) std::cerr << "d="<<d.first<<","<<d.second<<" k="<<k<<" u="<<u<<" ell="<<l.first<<","<<l.second<<" J="<<jac<<" L="<<g.value(types,l.first,l.second)<<" error="<<err<<'\n';
          check(err<1e-10,"density times Jacobian");
          const Energy ea=g.maximumA()*exp(-l.first), eb=g.maximumB()*exp(-l.second);
          const Energy2 s=sqr(g.maximumCMEnergy())*exp(-l.first-l.second);
          check(std::abs(s/(4*ea*eb)-1)<1e-12,"collision energy");
          check(std::abs(g.beamRapidity(l.first,l.second)-.5*log(ea/eb))<1e-12,"rapidity");
          auto xc=new_ptr(XComb());
          xc->lastP1P2(l);
          g.select(xc);
          // A subsequent rejected integration point must not change the accepted boost.
          const double rejected[] = {.2,.7}; double j=1;
          g.generateLL(rejected,j);
          LorentzMomentum p(ZERO,ZERO,ZERO,sqrt(s));
          p.transform(g.getBoost());
          check(std::abs(p.e()/(ea+eb)-1)<1e-12,"selected-event boost energy");
          check(std::abs((p.z()-(ea-eb))/(ea+eb))<1e-12,"selected-event boost momentum");
        }
        check(g.value(types,-.1,0)==0,"negative logarithm support");
        check(g.value(types,100,100)==0,"upper logarithm support");
      }
    }
    Cuts cuts;
    cuts.initialize(10000*GeV2,0.1);
    cuts.beamSMin(6400*GeV2);
    check(cuts.beamSLower()==6400*GeV2,"run lower support");
    cuts.beamKinematics(8100*GeV2,-.2);
    check(cuts.currentS()==8100*GeV2 && cuts.currentBeamY()==-.2,"current beam context");
    cuts.initEvent();
    check(cuts.currentS()==10000*GeV2 && cuts.currentBeamY()==.1,"reset beam context");
    auto boundCuts=new_ptr(Cuts(20*GeV));
    Repository::Register(boundCuts,"/BeamSpreadTestCuts");
    for (std::string setting:{"MHatMax 90","YHatMin -0.01","YHatMax 0.01"}) {
      const auto response=Repository::exec("set /BeamSpreadTestCuts:"+setting,std::cerr);
      check(response.find("Error")==std::string::npos,"set bounded cut interface");
    }
    boundCuts->initialize(10000*GeV2,0);
    boundCuts->beamSMin(85*85*GeV2);
    PDFCuts pdfBounds(*boundCuts,true,100*GeV);
    check(std::abs(pdfBounds.lMin())<1e-14,"upper mass cut must not exclude x=1 at low beam energy");
    boundCuts->beamKinematics(86*86*GeV2,.07);
    check(boundCuts->initSubProcess(86*86*GeV2,0),"sampled energy fractions under upper mass cut");
    auto original=new_ptr(GaussianBeamLuminosity());
    original->configure(46.1*GeV,45.1*GeV,.001,.003,7);
    std::stringstream buffer;
    { PersistentOStream stream(buffer); stream << original << boundCuts; }
    Ptr<GaussianBeamLuminosity>::pointer loaded;
    Ptr<Cuts>::pointer loadedCuts;
    { PersistentIStream stream(buffer); stream >> loaded >> loadedCuts; }
    check(bool(loaded) && bool(loadedCuts),"persistence types");
    const double variates[]={.173,.827}; double ja=1,jb=1;
    auto la=original->generateLL(variates,ja),lb=loaded->generateLL(variates,jb);
    check(la==lb && ja==jb,"save/load reproducibility");
    check(loaded->beamEMaxA()==loaded->maximumA(),"persistent derived maximum");
    check(loadedCuts->beamSLower()==85*85*GeV2 && loadedCuts->currentS()==10000*GeV2,
          "saved support and transient collision context");
    g.configure(45.6*GeV,45.6*GeV,.00221,0,6);
    gsl_function integrand;
    integrand.params=&g;
    integrand.function=[](double ell,void *ptr) {
      return static_cast<GaussianBeamLuminosity*>(ptr)->value(cPDPair(),ell,0);
    };
    auto workspace=gsl_integration_workspace_alloc(1000);
    double integral=0,error=0;
    const double limit=log1p(6*.00221)-log1p(-6*.00221);
    const int status=gsl_integration_qag(&integrand,0,limit,1e-12,1e-12,1000,6,workspace,&integral,&error);
    gsl_integration_workspace_free(workspace);
    check(status==0 && std::abs(integral-1)<1e-10,"independent logarithmic-density normalization");
    for(double u:{.2,.4,.6,.8}) {
      const double h=1e-4;
      auto f=[&](double x) { double j=1;return g.generateLL(&x,j).first; };
      double jac=1;g.generateLL(&u,jac);
      const double derivative=(f(u-2*h)-8*f(u-h)+8*f(u+h)-f(u+2*h))/(12*h);
      check(std::abs(derivative/jac+1)<1e-10,"numerical mapping derivative");
    }
    for(auto invalid:std::vector<std::array<double,5>>{{{-1,45.6,.01,.01,6}},
        {{45.6,45.6,-.01,0,6}},{{45.6,45.6,.2,0,6}},{{45.6,45.6,0,0,0}},
        {{45.6,45.6,std::numeric_limits<double>::quiet_NaN(),0,6}}}) {
      bool failed=false;
      try {g.configure(invalid[0]*GeV,invalid[1]*GeV,invalid[2],invalid[3],invalid[4]);}
      catch(const Exception& error) {error.handle();failed=true;}
      check(failed,"invalid parameter rejection");
    }
    std::cout << "PASS: mapping, support, zero widths, selected-event boost, cut context; max(JL-1)=" << maxError << '\n';
    return 0;
  } catch (const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
}
