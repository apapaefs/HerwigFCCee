// Laboratory geometry, deterministic flights, continuation and persistence.
#include "ThePEG/Handlers/DecayHandler.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/GaussianVertexDistribution.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Handlers/SamplerBase.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/StandardRandom.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <random>
using namespace ThePEG;
namespace {
void check(bool x,const char *s) {if(!x) throw std::runtime_error(s);}
double distance(const LorentzPoint &a,const LorentzPoint &b) {
  return (abs(a.x()-b.x())+abs(a.y()-b.y())+abs(a.z()-b.z())+abs(a.t()-b.t()))/mm;
}
template<class T> std::string saved(const T &p) {
  std::stringstream b; {PersistentOStream s(b);s<<p;} return b.str();
}
class ResumeHandler:public StandardEventHandler {
public:
  void prepare(tEGPtr g,tEventPtr e) {
    setGenerator(g);currentEvent(e);currentCollision(e->primaryCollision());
    currentStep(e->primaryCollision()->finalStep());
  }
  void pause() {
    auto h=new_ptr(Hint());h->stop(true);performStep(tStepHdlPtr(),h);
  }
};
DescribeNoPIOClass<ResumeHandler,StandardEventHandler> describeResume("GeometryTest::ResumeHandler","");

// Construct the reference with the stable boost constructor. In particular,
// boost() itself loses small-beta spin precision and cannot be the reference.
void checkSmallBoostPersistence() {
  double maximumVector=0.,maximumSpin=0.;unsigned int cases=0;
  auto roundTrip=[&](const LorentzRotation &reference) {
    for(bool complete : {false,true}) {
      auto event=new_ptr(Event(PPair()));event->geometry(reference,LorentzPoint());
      if(complete) event->finishGeometry();
      std::stringstream buffer;{PersistentOStream out(buffer);out<<event;}
      EventPtr loaded;{PersistentIStream in(buffer);in>>loaded;}
      check(loaded->hasGeometry() && loaded->geometryComplete()==complete,
            "small-boost persistence changed completion state");
      const auto &actual=loaded->laboratoryTransform();
      double vectorError=0.,spinError=0.;
      for(int i=0;i<4;++i) {
        LorentzMomentum a(i==0?GeV:ZERO,i==1?GeV:ZERO,i==2?GeV:ZERO,i==3?GeV:ZERO),b=a;
        a.transform(reference);b.transform(actual);
        const double deviation=std::max({abs(a.x()-b.x())/GeV,
                     abs(a.y()-b.y())/GeV,abs(a.z()-b.z())/GeV,abs(a.t()-b.t())/GeV});
        vectorError=std::max(vectorError,deviation);
        for(int j=0;j<4;++j)
          spinError=std::max(spinError,abs(reference.half()(i,j)-actual.half()(i,j)));
      }
      maximumVector=std::max(maximumVector,vectorError);
      maximumSpin=std::max(maximumSpin,spinError);++cases;
      if(vectorError>=1e-14 || spinError>=1e-14) {
        std::ostringstream error;
        error<<"small-boost persistence: vector error "<<vectorError
             <<", spin-1/2 error "<<spinError<<" (tolerance 1e-14)";
        throw std::runtime_error(error.str());
      }
    }
  };
  // This case failed by about 5e-9 when the loader used boost().
  roundTrip(LorentzRotation(0.,0.,1.01e-8));
  for(double beta : {0.,1e-12,1e-9,1.01e-8,1e-7,1e-6,1e-4,1e-3,.1,.5,.9})
    for(double sign : {-1.,1.})
      for(const Axis &axis : {Axis(0.,0.,1.),Axis(1./3.,2./3.,-2./3.)})
        for(double angle : {0.,1e-12,Constants::pi,2*Constants::pi,
                            2*Constants::pi+1e-8,3*Constants::pi}) {
          LorentzRotation rotation;rotation.setRotate(angle,Axis(1.,2.,-3.));
          roundTrip(LorentzRotation(sign*beta*axis)*rotation);
        }
  // A local deterministic engine leaves the event generator's RNG untouched.
  std::mt19937_64 random(781028);std::uniform_real_distribution<double> flat(0.,1.);
  auto axis=[&]() {
    const double z=2*flat(random)-1,phi=Constants::twopi*flat(random),r=sqrt(1-z*z);
    return Axis(r*cos(phi),r*sin(phi),z);
  };
  for(int i=0;i<2048;++i) {
    const double beta=pow(10.,-9+log10(.9e9)*flat(random));
    const Axis direction=axis(),rotationAxis=axis();
    LorentzRotation rotation;rotation.setRotate(2*Constants::twopi*flat(random),rotationAxis);
    roundTrip(LorentzRotation(beta*direction)*rotation);
  }
  std::cout<<"PASS: "<<cases<<" small-boost persistence round trips; max vector error "
           <<maximumVector<<", max spin-1/2 error "<<maximumSpin<<" (tolerance 1e-14)\n";
}
}
int main() {
try {
  checkSmallBoostPersistence();
  auto rng=new_ptr(StandardRandom());rng->setSeed(78123);UseRandom random(rng);
  auto data=ParticleData::Create(991001,"geometry-test");
  data->stable(false);
  auto particle=[&](double px,double pz,double ct) {
    auto p=data->produceParticle(LorentzMomentum(px*GeV,ZERO,pz*GeV,sqrt(1+px*px+pz*pz)*GeV));
    p->properLifeLength(ct*mm);return p;
  };
  auto p=particle(3,4,2);
  check(distance(p->lifeLength(),LorentzPoint(6*mm,ZERO,8*mm,2*sqrt(26.)*mm))<1e-10,"flight construction");
  const auto before=saved(rng);
  p->set5Momentum(Lorentz5Momentum(GeV,ZERO,3*GeV,4*GeV)); // m=sqrt(6) GeV
  p->synchronizeLifeLength();
  check(abs(p->lifeLength().tau()/mm-2)<1e-10,"proper lifetime changed");
  check(abs(p->lifeLength().x()/mm-2/sqrt(6.))<1e-10,"stale flight direction");
  check(saved(rng)==before,"flight reconstruction drew random numbers");
  auto photon=particle(0,0,0);photon->set5Momentum(Lorentz5Momentum(ZERO,ZERO,GeV,GeV));
  photon->synchronizeLifeLength();
  bool failed=false;try {photon->properLifeLength(mm);}catch(const Exception &e){e.handle();failed=true;}
  check(failed,"nonzero massless flight accepted");
  photon->properLifeLength(ZERO);
  auto policy=new_ptr(DecayHandler());policy->cylinder(true,5*mm,7*mm);
  for (double sign : {-1.,1.}) {
    check(policy->insideCylinder(LorentzPoint(3*mm,4*mm,sign*7*mm,ZERO)),"boundary rejected");
    check(!policy->insideCylinder(LorentzPoint(3*mm,4.00001*mm,ZERO,ZERO)),"radial exterior accepted");
    check(!policy->insideCylinder(LorentzPoint(ZERO,ZERO,sign*7.00001*mm,ZERO)),"axial exterior accepted");
  }
  vector<LorentzRotation> transforms;
  for(double y : {0.,.5*log(46.1/45.1),.5*log(60./35.)})
    transforms.push_back(LorentzRotation(0,0,tanh(y)));
  for(double angle : {0.,1e-12,Constants::pi-1e-12,Constants::pi,
                       Constants::pi+1e-12,2*Constants::pi,2*Constants::pi+1e-12}) {
    LorentzRotation transform;
    transform.setRotate(angle,Axis(1.,2.,-3.));
    transforms.push_back(transform); // includes the spin-1/2 minus identity
    transform.boost(.2,-.1,.3);
    transform.rotateY(.23);
    transform.boost(-.1,.3,.2); // non-collinear composition
    transforms.push_back(transform);
  }
  unsigned int context=0;
  for(const auto & boost : transforms) {
    const LorentzPoint ip(2*mm,-3*mm,4*mm,10*mm);
    const PPair beams(particle(0,40,0),particle(0,-40,0));
    auto event=new_ptr(Event(beams,tcEventBasePtr(),"geometry",9,-.75));
    event->optionalWeight("signed",-.125);
    auto step=event->newStep();
    const PPair incoming(particle(0,40,0),particle(0,-40,0));
    auto group=new_ptr(SubProcessGroup(incoming));
    auto q=particle(2,3,1);group->addOutgoing(q,false);
    auto dep=new_ptr(SubProcess(incoming)); // share incoming pointers (projector)
    auto detached=particle(-2,-3,2);dep->addOutgoing(detached,false);
    group->add(dep);step->addSubProcess(group);
    event->geometry(boost,ip);
    q->setLabVertex(LorentzPoint(mm,ZERO,ZERO,ZERO));
    auto expected=q->labDecayVertex();expected.transform(boost);expected+=ip;
    check(distance(event->laboratoryPoint(q->labDecayVertex()),expected)<1e-10,"candidate frame");
    check(policy->acceptsDecay(q,event.operator->())==policy->insideCylinder(expected),"rotated cylinder decision");
    auto cloned=event->clone();
    auto clonedGroup=dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(cloned->primaryCollision()->subProcesses()[0]);
    check(clonedGroup && clonedGroup->dependent()[0]!=dep,"dependent clone aliases source");
    check(clonedGroup->incoming()==clonedGroup->dependent()[0]->incoming(),"projected sharing lost");
    check(cloned->optionalWeight("signed")==-.125,"clone lost signed weight");
    std::stringstream buffer;{PersistentOStream s(buffer);s<<event;}
    EventPtr loaded;{PersistentIStream s(buffer);s>>loaded;}
    for(int i=0;i<4;++i) {
      LorentzMomentum a(i==0?GeV:ZERO,i==1?GeV:ZERO,i==2?GeV:ZERO,i==3?GeV:ZERO),b=a;
      a.transform(boost);b.transform(loaded->laboratoryTransform());
      check((abs(a.x()-b.x())+abs(a.y()-b.y())+abs(a.z()-b.z())+abs(a.t()-b.t()))/GeV<1e-10,"saved Lorentz matrix");
      for(int j=0;j<4;++j)
        check(abs(boost.half()(i,j)-loaded->laboratoryTransform().half()(i,j))<1e-10,"saved spin-1/2 matrix or sign");
    }
    auto handler=new_ptr(ResumeHandler());auto generator=new_ptr(EventGenerator());
    const std::string handlerName="/GeometryResumeHandler"+std::to_string(context++);
    Repository::Register(handler,handlerName);
    const auto response=Repository::exec("set "+handlerName+":ConsistencyLevel Never",std::cerr);
    check(response.find("Error")==std::string::npos,"configure continuation fixture");
    handler->prepare(generator,event);
    bool stopped=false;try {handler->pause();}catch(const Stop &){stopped=true;}
    check(stopped && !event->geometryComplete(),"Hint stop did not defer geometry");
    const auto randomState=saved(rng);
    handler->continueEvent();
    check(event->geometryComplete(),"continueEvent did not finish geometry");
    check(saved(rng)==randomState,"resume resampled geometry");
    check(distance(q->labDecayVertex(),expected)<1e-10,"completed decay position");
    check(distance(detached->labVertex(),ip)<1e-10,"detached translation");
    check(distance(beams.first->labVertex(),ip)<1e-10,"beam translation");
    const auto after=q->labDecayVertex();event->finishGeometry();
    check(distance(q->labDecayVertex(),after)==0,"double geometry transform");
    loaded->finishGeometry();cloned->finishGeometry();
    check(!clonedGroup->outgoing().empty(),"clone missing outgoing particles");
    auto cq=clonedGroup->outgoing()[0];check(bool(cq),"null clone outgoing");
    check(distance(cq->labDecayVertex(),expected)<1e-10,"clone frame or collision origin");
    check(distance(q->labDecayVertex(),after)==0,"clone changed source geometry");
    check(bool(loaded->primaryCollision()),"loaded missing collision");
    check(!loaded->primaryCollision()->subProcesses().empty(),"loaded missing subprocesses");
    auto ls=loaded->primaryCollision()->subProcesses()[0];check(bool(ls),"null loaded subprocess");
    check(!ls->outgoing().empty(),"loaded empty outgoing");
    auto lq=ls->outgoing()[0];check(bool(lq),"null loaded outgoing");
    check(distance(lq->labDecayVertex(),expected)<1e-10,"save/load geometry");

  }
  {
    auto rotated=new_ptr(Event(PPair()));LorentzRotation rotation;
    rotation.rotateY(.5*Constants::pi);rotated->geometry(rotation,LorentzPoint());
    auto candidate=particle(0,6,1);
    check(policy->insideCylinder(candidate->labDecayVertex()),"rotation edge fixture");
    check(!policy->acceptsDecay(candidate,rotated.operator->()),"rotation did not move decay across radial edge");
  }
  auto e=new_ptr(Event(PPair(),tcEventBasePtr(),"decay-policy",1));
  e->geometry(LorentzRotation(),LorentzPoint(4*mm,ZERO,ZERO,ZERO));
  auto inside=particle(0,0,100); // at rest: time is unrestricted
  check(policy->acceptsDecay(inside,e.operator->()),"rest decay rejected");
  auto outside=particle(2,0,1);
  const auto state=saved(rng);
  check(!policy->acceptsDecay(outside,e.operator->()),"displaced decay accepted");
  check(outside->children().empty() && !data->stable(),"geometry changed particle species/history");
  check(!policy->acceptsDecay(outside,e.operator->()),"repeated rejection changed decision");
  check(saved(rng)==state,"geometry rejection drew random numbers");
  auto vertex=new_ptr(GaussianVertexDistribution());
  check(distance(vertex->sample(),LorentzPoint())==0 && saved(rng)==state,"zero widths drew random numbers");
  vertex->configure({{mm,2*mm,3*mm,4*mm}},{{2*mm,3*mm,4*mm,5*mm}},2.,3.);
  double sums[4]={},squares[4]={};const int n=100000;
  for(int i=0;i<n;++i) {
    auto v=vertex->sample();double x[4]={(v.x()/mm-1)/2,(v.y()/mm-2)/3,(v.z()/mm-3)/4,(v.t()/mm-4)/5};
    check(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<=4 && abs(x[3])<=3,"sampler support");
    for(int j=0;j<4;++j){sums[j]+=x[j];squares[j]+=x[j]*x[j];}
  }
  // E[r^2]/3 for a three-dimensional normal truncated at radius 2.
  const double a=2.,f=erf(a/sqrt(2.))-sqrt(2/M_PI)*a*exp(-a*a/2);
  const double variance=1-sqrt(2/M_PI)*pow(a,3)*exp(-a*a/2)/(3*f);
  for(int j=0;j<4;++j) {
    const double v=j<3?variance:1-6*exp(-4.5)/(sqrt(2*M_PI)*erf(3/sqrt(2.)));
    check(abs(sums[j]/n)<6*sqrt(v/n),"sampler mean");
    check(abs(squares[j]/n-v)<6*sqrt(2./n),"truncated sampler variance");
  }
  std::stringstream vb;{PersistentOStream s(vb);s<<vertex<<rng;}
  Ptr<GaussianVertexDistribution>::pointer v2;Ptr<StandardRandom>::pointer r2;
  {PersistentIStream s(vb);s>>v2>>r2;}
  auto sample=vertex->sample();{UseRandom other(r2);check(distance(sample,v2->sample())==0,"vertex persistence/RNG reproducibility");}
  std::cout<<"PASS: cylinder boundaries, flights, full lab rotations/boosts and spin sign, displaced/rest particles, projected clone, save/load, Hint stop/continueEvent, vertex support and moments\n";
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
