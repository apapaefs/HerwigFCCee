// -*- C++ -*-
//
// DecayHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayHandler class.
//

#include "DecayHandler.h"
#include "Hint.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <cmath>

using namespace ThePEG;

DecayHandler::~DecayHandler() {}

IBPtr DecayHandler::clone() const {
  return new_ptr(*this);
}

IBPtr DecayHandler::fullclone() const {
  return new_ptr(*this);
}

void DecayHandler::
handle(EventHandler &, const tPVector & tagged,
       const Hint &) {
  // First go through to see which of the tagged particles should be
  // decayed: Exit if none are found.
  tPVector parents;
  for ( int i = 0, N = tagged.size(); i < N; ++i )
    if ( tagged[i] && !tagged[i]->data().stable() )
      parents.push_back(tagged[i]);

  if ( parents.empty() ) return;

  // Create a new step, decay all particles and add their children in
  // the new step.
  for ( int i = 0, N = parents.size(); i < N; ++i )
    performDecay(newStep()->find(parents[i]->final()), *newStep());
}

void DecayHandler::
performDecay(tPPtr parent, Step & s) const {
  parent->synchronizeLifeLength();
  if (!acceptsDecay(parent,s.collision() ? s.collision()->event().operator->() : nullptr)) return;
  ParticleVector children = Decayer::DecayParticle(parent, s, maxLoop());
  for (const auto & child : children) {
    child->setLabVertex(parent->labDecayVertex());
    child->synchronizeLifeLength();
  }
  for ( int i = 0, N = children.size(); i < N; ++i )
    if ( !children[i]->data().stable() ) performDecay(children[i], s);
}

void DecayHandler::persistentOutput(PersistentOStream & os) const {
  os << theMaxLoop << ounit(theMaxLifeTime, mm) << theLifeTimeOption;
  os << theLimitCylinder << ounit(theCylinderRadius,mm) << ounit(theCylinderHalfLength,mm);
}

void DecayHandler::persistentInput(PersistentIStream & is, int version) {
  is >> theMaxLoop >> iunit(theMaxLifeTime, mm) >> theLifeTimeOption;
  if (version >= 1)
    is >> theLimitCylinder >> iunit(theCylinderRadius,mm) >> iunit(theCylinderHalfLength,mm);
}

void DecayHandler::cylinder(bool enabled, Length radius, Length halfLength) {
  theLimitCylinder=enabled; theCylinderRadius=radius; theCylinderHalfLength=halfLength;
  validateGeometry();
}

void DecayHandler::validateGeometry() const {
  if (theLimitCylinder && (!(theCylinderRadius>ZERO) || !(theCylinderHalfLength>ZERO) ||
      !std::isfinite(theCylinderRadius/mm) || !std::isfinite(theCylinderHalfLength/mm)))
    throw Exception() << "LimitCylinder requires finite positive CylinderRadius and "
                      << "CylinderHalfLength." << Exception::runerror;
  if (theLimitCylinder && generator() &&
      !dynamic_ptr_cast<Ptr<StandardEventHandler>::tcptr>(generator()->eventHandler()))
    throw Exception() << "LimitCylinder requires StandardEventHandler. "
                      << "LesHouchesEventHandler/LHE input is unsupported."
                      << Exception::runerror;
}

void DecayHandler::doinit() { StepHandler::doinit(); validateGeometry(); }
void DecayHandler::doinitrun() { StepHandler::doinitrun(); validateGeometry(); }

bool DecayHandler::insideCylinder(const LorentzPoint & point) const {
  if (!std::isfinite(point.x()/mm) || !std::isfinite(point.y()/mm) ||
      !std::isfinite(point.z()/mm) || !std::isfinite(point.t()/mm))
    throw Exception() << "Nonfinite candidate decay vertex." << Exception::eventerror;
  return std::hypot(point.x()/mm,point.y()/mm)<=theCylinderRadius/mm &&
         abs(point.z())<=theCylinderHalfLength;
}

bool DecayHandler::acceptsDecay(tPPtr parent, const Event * event) const {
  if (parent->decayDisposition()==Particle::OutsideCylinder ||
      parent->decayDisposition()==Particle::AboveMaxLifeTime) return false;
  if (maxLifeTime()>=ZERO &&
      (lifeTimeOption() ? parent->lifeLength().tau() : parent->data().cTau())>maxLifeTime()) {
    parent->decayDisposition(Particle::AboveMaxLifeTime);
    return false;
  }
  if (theLimitCylinder) {
    if (!event || !event->hasGeometry())
      throw Exception() << "LimitCylinder has no laboratory geometry context; "
                        << "StandardEventHandler is required and LHE is unsupported."
                        << Exception::runerror;
    if (!insideCylinder(event->laboratoryPoint(parent->labDecayVertex()))) {
      parent->decayDisposition(Particle::OutsideCylinder);
      return false;
    }
  }
  parent->decayDisposition(Particle::DecayAllowed);
  return true;
}

ClassDescription<DecayHandler> DecayHandler::initDecayHandler;

void DecayHandler::Init() {

  static Switch<DecayHandler,bool> limitCylinder
    ("LimitCylinder","Only decay particles whose sampled laboratory decay point is inside the cylinder.",
     &DecayHandler::theLimitCylinder,false,false,false);
  static SwitchOption cylinderOn(limitCylinder,"On","Enable the laboratory cylinder.",true);
  static SwitchOption cylinderOff(limitCylinder,"Off","Disable the laboratory cylinder.",false);
  static Parameter<DecayHandler,Length> radius
    ("CylinderRadius","Laboratory cylinder radius in mm; the surface is included.",
     &DecayHandler::theCylinderRadius,mm,2250*mm,ZERO,ZERO,false,false,Interface::lowerlim);
  static Parameter<DecayHandler,Length> halfLength
    ("CylinderHalfLength","Laboratory cylinder half-length in mm; both end surfaces are included.",
     &DecayHandler::theCylinderHalfLength,mm,2500*mm,ZERO,ZERO,false,false,Interface::lowerlim);

  static ClassDocumentation<DecayHandler> documentation
    ("This is the main class handling the decay of unstable particles. Note "
     "that the actual decays will be performed by objects of the "
     "ThePEG::Decayer class.");

  static Parameter<DecayHandler,long> interfaceMaxLoop
    ("MaxLoop",
     "The maximum number of attempts per event when selecting a decay channel.",
     &DecayHandler::theMaxLoop, 100000, 100, 100000000, true, false, true);

  static Parameter<DecayHandler,Length> interfaceMaxLifeTime
    ("MaxLifeTime",
     "The maximum lifetime (c*tau) in units of mm. Particles with longer "
     "lifetimes than this will not be decayed.",
     &DecayHandler::theMaxLifeTime, mm, -1.0*mm, -1.0*mm, ZERO,
     true, false, Interface::lowerlim);


  static Switch<DecayHandler,bool> interfaceLifeTimeOption
    ("LifeTimeOption",
     "Option for how the maximum life time is interpreted",
     &DecayHandler::theLifeTimeOption, false, false, false);
  static SwitchOption interfaceLifeTimeOptionAverage
    (interfaceLifeTimeOption,
     "Average",
     "Cut on the average lifetime of the particle type",
     false);
  static SwitchOption interfaceLifeTimeOptionReal
    (interfaceLifeTimeOption,
     "Real",
     "Cut on the lifetime generated for the given instance",
     true);

}

