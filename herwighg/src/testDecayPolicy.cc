// Native mixing limits and complete-tree pruning, independent of EvtGen.
#include "Herwig/Decay/HwDecayHandler.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/MixedParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/Repository/StandardRandom.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Command.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using namespace ThePEG;
namespace {
void check(bool x,const char *m) {if(!x)throw std::runtime_error(m);}
class PairDecayer:public Decayer {
public:
  bool accept(const DecayMode &) const override {return true;}
  ParticleVector decay(const DecayMode &mode,const Particle &parent) const override {
    auto result=getChildren(mode,parent);
    check(result.size()==2,"two body fixture");
    const double c=2*UseRandom::rnd()-1,s=sqrt(1-c*c);
    const Energy e=parent.mass()/2;
    Lorentz5Momentum a(e*s,ZERO,e*c,e),b(-e*s,ZERO,-e*c,e);
    a.boost(parent.momentum().boostVector());b.boost(parent.momentum().boostVector());
    result[0]->set5Momentum(a);result[1]->set5Momentum(b);
    return result;
  }
  IBPtr clone()const override{return new_ptr(*this);}
  IBPtr fullclone()const override{return new_ptr(*this);}
  static void Init(){}
};
DescribeNoPIOClass<PairDecayer,Decayer> description("GeometryTest::PairDecayer","");
std::string randomState(const RanGenPtr &r) {
  std::stringstream b;
  {PersistentOStream s(b);r->RandomGenerator::persistentOutput(s);r->persistentOutput(s);}
  return b.str();
}
void command(const std::string &s) {
  auto result=Repository::exec(s,std::cerr);
  check(result.find("Error")==std::string::npos,"repository command");
}
}
int testDecayPolicy() {
try {
  auto handler=new_ptr(Herwig::HwDecayHandler());Repository::Register(handler,"/DecayPolicyTest");
  auto mixed=MixedParticleData::Create(511,"MixFixture","AntiMixFixture");
  auto photon=ParticleData::Create(22,"PhotonFixture");photon->stable(true);
  Repository::Register(photon,"/PhotonFixture");
  Repository::Register(mixed.first,"/MixFixture");Repository::Register(mixed.second,"/AntiMixFixture");
  mixed.first->synchronized(true);
  auto decayer=new_ptr(PairDecayer());Repository::Register(decayer,"/PairFixture");
  command("set /PhotonFixture:NominalMass 0");
  for(const std::string name:{"/MixFixture","/AntiMixFixture"}) {
    command("set "+name+":NominalMass 5");
    std::ostringstream w;w<<std::setprecision(17)<<hbarc/(mm*GeV);
    command("set "+name+":Width "+w.str());
  }
  mixed.first->cTau(mm);mixed.second->cTau(mm);
  {std::ostringstream delta;delta<<std::setprecision(17)<<.8*hbarc/(mm*GeV);
   command("set /MixFixture:DeltaM "+delta.str());}
  command("decaymode /MixFixture->/PhotonFixture,/PhotonFixture; 1.0 1 /PairFixture");
  mixed.first->stable(false);mixed.second->stable(false);
  mixed.first->reset();mixed.second->reset();
  mixed.first->update();mixed.second->update();
  check(!mixed.first->decaySelector().empty(),"empty fixture decay selector");
  mixed.first->reset();mixed.second->reset();
  mixed.first->init();mixed.second->init();
  auto makeParent=[&](tcPDPtr pd,Energy px,Energy mass,Length ct) {
    auto p=pd->produceParticle(LorentzMomentum(px,ZERO,ZERO,sqrt(px*px+mass*mass)));
    p->properLifeLength(ct);return p;
  };
  // Exactly the same random stream with a disabled or all-accepting cylinder.
  for(const std::string option:{"Average","Real"})
    for(double limit:{-1.,0.,.5,2.}) {
      command("set /DecayPolicyTest:LifeTimeOption "+option);
      command("set /DecayPolicyTest:MaxLifeTime "+std::to_string(limit)+"*mm");
      std::string history[2],states[2];
      for(int cylinder=0;cylinder<2;++cylinder) {
        auto rng=new_ptr(StandardRandom());rng->setSeed(18273);UseRandom random(rng);
        handler->cylinder(cylinder,1e10*mm,1e10*mm);
        auto event=new_ptr(Event(PPair(),tcEventBasePtr(),"mix",1));auto step=event->newStep();
        event->geometry(LorentzRotation(0,0,.25),LorentzPoint());
        auto p=makeParent(mixed.first,5*GeV,5*GeV,ZERO);step->addParticle(p);
        handler->performDecay(p,*step);
        if(limit==0 || (option=="Average" && limit==.5)) {
          check(p->children().empty(),"lifetime rejection inserted flavour bookkeeping");
          check(p->decayDisposition()==Particle::AboveMaxLifeTime,"wrong rejection reason");
          const auto state=randomState(rng);handler->performDecay(p,*step);
          check(state==randomState(rng),"rejected mixed particle resampled lifetime");
        }
        if(limit<0 || (option=="Average" && limit==2)) check(!p->children().empty(),"allowed mixed decay missing");
        std::stringstream b;
        ParticleVector queue(1,p);
        for(size_t i=0;i<queue.size();++i) {
          const auto x=queue[i];b<<x->id()<<' '<<x->momentum()/GeV<<' '<<x->lifeLength().tau()/mm<<' '<<x->children().size()<<'\n';
          queue.insert(queue.end(),x->children().begin(),x->children().end());
        }
        history[cylinder]=b.str();states[cylinder]=randomState(rng);
      }
      check(history[0]==history[1] && states[0]==states[1],"all-accepting cylinder changed mixed history or RNG");
    }
  command("set /DecayPolicyTest:MaxLifeTime -1*mm");
  auto rng=new_ptr(StandardRandom());rng->setSeed(36152);UseRandom random(rng);
  int flipped=0,unflipped=0;
  handler->cylinder(false,mm,mm);
  for(int i=0;i<64;++i) {
    auto event=new_ptr(Event(PPair(),tcEventBasePtr(),"flavour",i));auto step=event->newStep();
    auto p=makeParent(mixed.first,5*GeV,5*GeV,ZERO);step->addParticle(p);handler->performDecay(p,*step);
    check(p->children().size()==1 && p->lifeLength().tau()==ZERO && p->lifeLength().t()==ZERO,
          "flavour bookkeeping has a nonzero flight");
    if(p->children()[0]->id()==p->id())++unflipped;else ++flipped;
    check(p->children()[0]->lifeLength().tau()>ZERO,"mixed lifetime was lost");
  }
  check(flipped>0 && unflipped>0,"mixing fixture did not cover both flavour outcomes");
  handler->cylinder(true,1e-12*mm,1e-12*mm);
  auto mixEvent=new_ptr(Event(PPair(),tcEventBasePtr(),"mix-reject",1));auto mixStep=mixEvent->newStep();
  mixEvent->geometry(LorentzRotation(),LorentzPoint());
  auto rejected=makeParent(mixed.first,5*GeV,5*GeV,ZERO);mixStep->addParticle(rejected);
  handler->performDecay(rejected,*mixStep);
  check(rejected->children().empty() && rejected->decayDisposition()==Particle::OutsideCylinder,
        "mixed cylinder rejection inserted a flavour transition");
  // Complete pre-generated tree: change a child's momentum after creating its
  // lifetime. The decision must use the new direction and prune only its branch.
  auto parentData=ParticleData::Create(991002,"ChainFixture");parentData->stable(false);
  auto leafData=ParticleData::Create(991003,"LeafFixture");leafData->stable(true);
  for(bool limited:{false,true}) {
    handler->cylinder(limited,mm,5*mm);
    auto event=new_ptr(Event(PPair(),tcEventBasePtr(),"chain",1));auto step=event->newStep();
    event->geometry(LorentzRotation(),LorentzPoint(.25*mm,ZERO,ZERO,ZERO));
    auto root=makeParent(parentData,ZERO,4*GeV,ZERO);
    auto child=makeParent(parentData,ZERO,2*GeV,2*mm);
    auto leaf1=makeParent(leafData,GeV,GeV,ZERO),leaf2=makeParent(leafData,GeV,GeV,ZERO);
    child->addChild(leaf1);child->addChild(leaf2);root->addChild(child);
    step->addParticle(root);
    child->set5Momentum(Lorentz5Momentum(2*GeV,ZERO,ZERO,2*sqrt(2.)*GeV));
    const auto state=randomState(rng);handler->addDecayedParticle(root,*step);
    check(abs(child->lifeLength().x()/mm-2)<1e-10 && child->lifeLength().tau()==2*mm,"pre-decayed flight not repaired");
    check(state==randomState(rng),"tree policy drew random numbers");
    if(limited) {
      check(child->children().empty() && child->decayDisposition()==Particle::OutsideCylinder,"tree branch not pruned");
      check(step->particles().count(child)==1 && !step->all().count(leaf1),"pruned descendants entered event");
      check(leaf1->parents().empty() && leaf2->parents().empty(),"dangling pruned ancestry");
    } else {
      check(child->children().size()==2,"unrestricted tree pruned");
      check(abs(leaf1->labVertex().x()/mm-2)<1e-10,"descendant production position stale");
    }
  }
  std::cout<<"PASS: mixed-particle Average/Real limits, pre-bookkeeping rejection, native huge-cylinder/RNG invariance, complete-tree pruning and changed momentum\n";
  return 0;
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}

class GeometryPolicyCheck:public Interfaced {
public:
  std::string run(std::string) {
    if(testDecayPolicy()!=0)throw Exception()<<"Native geometry test failed."<<Exception::runerror;
    return "";
  }
  IBPtr clone() const override{return new_ptr(*this);}
  IBPtr fullclone() const override{return new_ptr(*this);}
  static void Init(){
    static Command<GeometryPolicyCheck> run("Run","Run the native decay policy regression tests.",&GeometryPolicyCheck::run,false);
  }
};
DescribeNoPIOClass<GeometryPolicyCheck,Interfaced> policyCheckDescription
  ("GeometryTest::PolicyCheck","TestDecayPolicy.so");
