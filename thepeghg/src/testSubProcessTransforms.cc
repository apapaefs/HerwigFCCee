// Detached NLO counterevents must follow frame changes exactly once.
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"
#include "ThePEG/PDT/ParticleData.h"
#include <iostream>
#include <map>
#include <stdexcept>

using namespace ThePEG;

namespace {
void require(bool condition, const char * message) {
  if (!condition) throw std::runtime_error(message);
}

double difference(const LorentzMomentum & a, const LorentzMomentum & b) {
  const Energy scale = max(GeV, max(abs(a.e()), abs(b.e())));
  return max(max(abs(a.x()-b.x()), abs(a.y()-b.y())),
             max(abs(a.z()-b.z()), abs(a.e()-b.e()))) / scale;
}

void checkTransform(bool collisionOnly, bool projected,
                    const LorentzRotation & rotation) {
  const auto data = ParticleData::Create(11, "transform-test");
  std::map<PPtr,LorentzMomentum> before;
  const auto particle = [&](double x, double z) {
    LorentzMomentum p(x*GeV, ZERO, z*GeV, sqrt(x*x+z*z+1)*GeV);
    PPtr result = data->produceParticle(p);
    before[result] = p;
    return result;
  };
  const PPair beams(particle(0,45.6),particle(0,-45.6));
  auto event = new_ptr(Event(beams,tcEventBasePtr(),"transform-test",1,1.));
  auto step = event->newStep();
  auto collision = event->primaryCollision();
  const PPair incoming(particle(0,35),particle(0,-35));
  auto head = new_ptr(SubProcessGroup(incoming));
  auto out1 = particle(12,23), out2 = particle(-12,-23);
  head->addOutgoing(out1,false);
  head->addOutgoing(out2,false);

  // Model a projector selecting a dependent configuration as the head.
  const PPair depIncoming = projected ? incoming :
    PPair(particle(0,35),particle(0,-35));
  auto dependent = new_ptr(SubProcess(depIncoming));
  dependent->addOutgoing(projected ? out1 : particle(15,20),false);
  dependent->addOutgoing(projected ? out2 : particle(-15,-20),false);
  dependent->addIntermediate(particle(2,4),false);
  dependent->groupWeight(-.75);
  head->add(dependent);

  // A nested group and repeated links exercise uniqueness across traversal.
  auto nested = new_ptr(SubProcessGroup(PPair(particle(0,30),particle(0,-30))));
  nested->addOutgoing(particle(9,17),false);
  nested->add(dependent);
  head->add(nested);
  step->addSubProcess(head);
  ParticleSet registered, collisionRegistered;
  event->select(inserter(registered,registered.end()),SelectAll());
  collision->select(inserter(collisionRegistered,collisionRegistered.end()),SelectAll());
  require(registered.size() < before.size(), "fixture has detached particles");

  if (collisionOnly) collision->transform(rotation);
  else event->transform(rotation);
  for (const auto & entry : before) {
    auto expected = entry.second;
    expected.transform(rotation);
    require(difference(entry.first->momentum(),expected)<1e-12,
            "particle missed a transform or was transformed twice");
  }
  ParticleSet eventAfter, collisionAfter;
  event->select(inserter(eventAfter,eventAfter.end()),SelectAll());
  collision->select(inserter(collisionAfter,collisionAfter.end()),SelectAll());
  require(eventAfter==registered && collisionAfter==collisionRegistered,
          "transform must not register counterevent particles");
  require(dependent->groupWeight()==-.75, "transform changed the subevent weight");

  if (collisionOnly) collision->transform(rotation.inverse());
  else event->transform(rotation.inverse());
  for (const auto & entry : before)
    require(difference(entry.first->momentum(),entry.second)<1e-12,
            "inverse transformation failed");
}
}

int main() {
  try {
    for (bool collisionOnly : {false,true})
      for (bool projected : {false,true}) {
        for (double y : {0.,.02,.5*log(30./920.)})
          checkTransform(collisionOnly,projected,LorentzRotation(0.,0.,tanh(y)));
        checkTransform(collisionOnly,projected,LorentzRotation(.1,-.2,.3));
      }
    std::cout << "PASS: event/collision transforms, detached and shared projected "
                 "particles, nested groups, identity, DIS and inverse boosts\n";
  } catch (const std::exception & error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
