// -*- C++ -*-
//
// Event.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined member functions of
// the Event class.
//
#include "Event.h"
#include "TransformSubProcesses.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/EventRecord/ParticleTraits.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <iostream>
#include <cmath>
#include <array>

using namespace ThePEG;

namespace {
// A signed unit quaternion also records the lift to the spin-1/2
// representation. Reading it from the spin-1 matrix would lose that sign.
std::array<double,4> laboratoryRotation(const LorentzRotation & transformation) {
  const Boost beta(transformation.xt()/transformation.tt(),
                   transformation.yt()/transformation.tt(),
                   transformation.zt()/transformation.tt());
  const auto rotation=LorentzRotation(-beta)*transformation;
  const auto & spin=rotation.half();
  std::array<double,4> q={{real(spin.s1s1()),-imag(spin.s1s2()),
                         -real(spin.s1s2()),-imag(spin.s1s1())}};
  const double norm=std::hypot(std::hypot(q[0],q[1]),std::hypot(q[2],q[3]));
  if (!(norm>0.) || !std::isfinite(norm))
    throw Exception() << "Invalid laboratory rotation." << Exception::runerror;
  for (auto & value:q) value/=norm;
  return q;
}

LorentzRotation laboratoryTransformation(double bx,double by,double bz,
                                         double qw,double qx,double qy,double qz) {
  const double norm=std::hypot(std::hypot(qw,qx),std::hypot(qy,qz));
  if (!(norm>0.) || !std::isfinite(norm))
    throw Exception() << "Invalid saved laboratory rotation." << Exception::runerror;
  qw/=norm; qx/=norm; qy/=norm; qz/=norm;
  const double sine=std::hypot(qx,std::hypot(qy,qz));
  // Identity and a full 2*pi spin rotation both have zero vector part.
  const Axis axis=sine>0. ? Axis(qx/sine,qy/sine,qz/sine) : Axis(1.,0.,0.);
  LorentzRotation rotation;
  rotation.setRotate(2.*atan2(sine,qw),axis);
  // The constructor avoids cancellation in boost()'s spin-1/2 matrix at
  // small beta. Preserve the order Lambda = B(beta) R.
  return LorentzRotation(bx,by,bz)*rotation;
}
}

Event::Event(const PPair & newIncoming, tcEventBasePtr newHandler,
	     string newName, long newNumber, double newWeight)
  : Named(newName), theIncoming(newIncoming), theHandler(newHandler),
    theNumber(newNumber), theWeight(newWeight), theParticleNumber(0) {
  addParticle(incoming().first);
  addParticle(incoming().second);
}

Event::Event(const Event & e)
  : Base(e), Named(e), 
    theIncoming(e.theIncoming), theCollisions(e.theCollisions),
    allSteps(e.allSteps), allSubProcesses(e.allSubProcesses),
    allParticles(e.allParticles), theHandler(e.theHandler),
    theNumber(e.theNumber), theWeight(e.theWeight),
    theParticleNumber(e.theParticleNumber),
    theHasGeometry(e.theHasGeometry), theGeometryComplete(e.theGeometryComplete),
    theLaboratoryBoost(e.theLaboratoryBoost), theInteractionPoint(e.theInteractionPoint) {
  theOptionalWeights = e.theOptionalWeights;
}

Event::~Event() {
  for ( int i = 0, N = theCollisions.size(); i < N; ++i )
    if ( theCollisions[i]->event() == this )
      theCollisions[i]->theEvent = tEventPtr();
  theIncoming = PPair();
  theCollisions.clear();
  allSteps.clear();
  allSubProcesses.clear();
  allParticles.clear();
  theHandler = tcEventBasePtr();
  theColourLines.clear();
  theNumber = -1;
  theWeight = 0.0;
}

void Event::setInfo(tcEventBasePtr newHandler, string newName,
		    long newNumber, double weight) {
  theHandler = newHandler;
  name(newName);
  theNumber = newNumber;
  theWeight = weight;
}

double Event::optionalWeight(const string& name) const {
  map<string,double>::const_iterator w =
    theOptionalWeights.find(name);
  if ( w == theOptionalWeights.end() )
    return 0.;
  return w->second;
}

void Event::optionalWeight(const string& name, double value) {
  theOptionalWeights[name] = value;
}

tCollPtr Event::newCollision() {
  theCollisions.push_back(new_ptr(Collision(incoming(), this)));
  return theCollisions.back();
}

tStepPtr Event::newStep() {
  if ( theCollisions.empty() ) newCollision();
  return theCollisions.back()->newStep();
}


void Event::addCollision(tCollPtr c) {
  if ( !c ) return;
  theCollisions.push_back(c);
  addParticles(c->all().begin(), c->all().end());
  allSubProcesses.insert(c->subProcesses().begin(), c->subProcesses().end()); 
}

void Event::addParticle(tPPtr p) {
  if ( !p ) return;
  if ( member(allParticles, p) ) return;
  allParticles.insert(p);
  p->number(++theParticleNumber);
}

void Event::transform(const LorentzRotation & r) {
  EventRecordDetail::transformParticles(allParticles, allSubProcesses.begin(),
                                       allSubProcesses.end(), r);
  for (const auto & collision : collisions()) {
    LorentzPoint position=collision->vertex();
    position.transform(r);
    collision->vertex(position);
  }
  if (theGeometryComplete) theInteractionPoint.transform(r);
}

void Event::translate(const LorentzPoint & offset) {
  if (offset.x()==ZERO && offset.y()==ZERO && offset.z()==ZERO && offset.t()==ZERO)
    return;
  // Only particles with a birth collision use collision-relative vertices.
  // Shift the collision origin for those, and the particle itself otherwise.
  const auto particles=EventRecordDetail::allParticles(allParticles,
                                    allSubProcesses.begin(),allSubProcesses.end());
  for (const auto & p : particles)
    if (!p->birthStep() || !p->birthStep()->collision())
      p->setVertex(p->vertex()+offset);
  for (const auto & collision : collisions())
    collision->vertex(collision->vertex()+offset);
  if (theGeometryComplete) theInteractionPoint+=offset;
}

void Event::geometry(const LorentzRotation & boost, const LorentzPoint & origin) {
  if (theHasGeometry)
    throw Exception() << "Event geometry was initialized twice." << Exception::runerror;
  // Keep the exact transformation in memory, including a custom luminosity
  // function's rotation. Persistence stores both its boost and spin rotation.
  const double bx=boost.one().xt()/boost.one().tt();
  const double by=boost.one().yt()/boost.one().tt();
  const double bz=boost.one().zt()/boost.one().tt();
  if (!std::isfinite(bx+by+bz) || bx*bx+by*by+bz*bz>=1.)
    throw Exception() << "Invalid laboratory boost." << Exception::runerror;
  for (int i=0;i<4;++i) {
    LorentzMomentum a(i==0?GeV:ZERO,i==1?GeV:ZERO,i==2?GeV:ZERO,i==3?GeV:ZERO);
    a.transform(boost);
    for (double value : {a.x()/GeV,a.y()/GeV,a.z()/GeV,a.t()/GeV})
      if (!std::isfinite(value))
        throw Exception() << "Nonfinite laboratory transformation." << Exception::runerror;
    for (int j=0;j<4;++j)
      if (!std::isfinite(real(boost.half()(i,j))) || !std::isfinite(imag(boost.half()(i,j))))
        throw Exception() << "Nonfinite laboratory spin transformation." << Exception::runerror;
  }
  for (Length value : {origin.x(),origin.y(),origin.z(),origin.t()})
    if (!std::isfinite(value/mm))
      throw Exception() << "Nonfinite interaction point." << Exception::runerror;
  theHasGeometry=true;
  theGeometryComplete=false;
  theLaboratoryBoost=boost;
  theInteractionPoint=origin;
}

LorentzPoint Event::laboratoryPoint(const LorentzPoint & point) const {
  if (!theHasGeometry)
    throw Exception() << "No laboratory geometry context: this feature requires "
                      << "StandardEventHandler; LHE input is unsupported."
                      << Exception::runerror;
  if (theGeometryComplete) return point;
  LorentzPoint result=point;
  result.transform(theLaboratoryBoost);
  return result+theInteractionPoint;
}

void Event::finishGeometry() {
  if (!theHasGeometry || theGeometryComplete) return;
  transform(theLaboratoryBoost);
  translate(theInteractionPoint);
  theGeometryComplete=true;
}

int Event::colourLineIndex(tcColinePtr line) const {
  ColourLineMap::const_iterator found = theColourLines.find(line);
  if ( found != theColourLines.end() ) return found->second;
  int index = theColourLines.size() + 1;
  theColourLines[line] = index;
  return index;
}

void Event::primaryCollision(tCollPtr c) {
  if ( !c ) return;
  if ( theCollisions.empty() ) theCollisions.push_back(c);
  else {
    if ( theCollisions[0] )
      for ( ParticleSet::const_iterator it = theCollisions[0]->all().begin();
	    it != theCollisions[0]->all().end(); ++it ) 
	allParticles.erase(*it);
    theCollisions[0] = c;
  }
  addParticles(c->all().begin(), c->all().end());
}

void Event::removeDecay(tPPtr p) {
  while ( !p->children().empty() ) removeParticle(p->children().back());
}

void Event::removeEntry(tPPtr p) {
  ParticleSet::iterator it = allParticles.find(p);
  if ( it == allParticles.end() ) return;
  for ( auto & cit : theCollisions ) cit->removeEntry(p);
  allParticles.erase(it);
}

void Event::removeParticle(tPPtr p) {
  if ( p->next() ) removeParticle(p->next());
  while ( !p->children().empty() ) removeParticle(p->children().back());
  removeEntry(p);
}

void Event::cleanSteps() {
  for ( auto & cit : theCollisions ) cit->cleanSteps();
}

EventPtr Event::clone() const {
  EventPtr newEvent = ptr_new<EventPtr>(*this);
  EventTranslationMap trans;
  trans[this] = newEvent;
  for ( CollisionVector::const_iterator cit = theCollisions.begin();
	cit != theCollisions.end(); ++cit ) trans[*cit] = (**cit).clone();
  ParticleSet detached;
  SubProcessSet subprocesses;
  for (const auto & sub : allSubProcesses)
    EventRecordDetail::collectDetachedParticles(sub,allParticles,detached,subprocesses);
  for (const auto & sub : subprocesses) trans[sub]=sub->clone();
  for ( StepSet::const_iterator sit = allSteps.begin();
        sit != allSteps.end(); ++sit ) trans[*sit] = (**sit).clone();
  for (const auto & collision : collisions())
    for (const auto & step : collision->steps())
      if (!trans.translate(step)) trans[step]=step->clone();
  for ( ParticleSet::const_iterator pit = allParticles.begin();
	pit != allParticles.end(); ++pit )
    trans[*pit] = (**pit).clone();
  for (const auto & particle : detached) trans[particle]=particle->clone();
  // Steps and projected groups can share objects. Rebinding a shared clone a
  // second time must be idempotent, not replace its links with null pointers.
  const auto translations=trans.map();
  for (const auto & entry : translations) trans[entry.second]=entry.second;
  newEvent->rebind(trans);
  for (const auto & particle : allParticles) trans.translate(particle)->rebind(trans);
  for (const auto & particle : detached) trans.translate(particle)->rebind(trans);
  return newEvent;
}

void Event::rebind(const EventTranslationMap & trans) {
  theIncoming.first = trans.translate(theIncoming.first);
  theIncoming.second = trans.translate(theIncoming.second);
  for ( CollisionVector::iterator cit = theCollisions.begin();
	cit != theCollisions.end(); ++cit )
    (*cit = trans.translate(*cit))->rebind(trans);
  SubProcessSet newSubProcesses;
  for ( SubProcessSet::const_iterator spit = allSubProcesses.begin();
	spit != allSubProcesses.end(); ++spit )
    newSubProcesses.insert(trans.translate(*spit));
  allSubProcesses.swap(newSubProcesses);
  StepSet newSteps;
  for ( StepSet::const_iterator sit = allSteps.begin();
	sit != allSteps.end(); ++sit )
    newSteps.insert(trans.translate(*sit));
  allSteps.swap(newSteps);
  ParticleSet newParticles;
  for ( ParticleSet::const_iterator pit = allParticles.begin();
	pit != allParticles.end(); ++pit )
    newParticles.insert(trans.translate(*pit));
  allParticles.swap(newParticles);
}
  
ostream & ThePEG::operator<<(ostream & os, const Event & e) {
  os << string(78, '*') << endl
     << "Event number " << e.number() << " (id: " << e.name() << ") ";
  if ( e.handler() ) os << "performed by "
			<<  EventConfig::nameHandler(e.handler());
  os << endl;
  os << "Weight: " << e.weight();
  if ( ! e.optionalWeights().empty() ) {
    os << " [ ";
    for ( auto w : e.optionalWeights() )
      os << w.first << ':' << w.second << ' ';
    os << " ]";
  }
  os << endl;
  for ( unsigned int i = 0; i < e.collisions().size(); ++i ) {
    os << string(78, '=') << endl;
    if ( e.collisions().size() != 1 ) {
      if ( i ) {
	os << "Secondary Collision " << i;
	if ( e.collisions()[i]->handler() )
	  os << " performed by "
	     << EventConfig::nameHandler(e.collisions()[i]->handler());
      } else
	os << "Primary Collision";
    }
    os << endl << *e.collisions()[i];
  }
  return os;
}

// Helpers for the Graphviz output
namespace {
  static const string header = "digraph test {\nrankdir=LR;\nranksep=1.5;\n";

  inline long startnode(tcPPtr p) {
    return p->parents().empty() ? -p->uniqueId : p->parents()[0]->uniqueId;
  }

  inline long endnode(tcPPtr p) {
    return p->children().empty() ? p->uniqueId : startnode( p->children()[0] );
  }

  static const char * colours[] = {
    "red",
    "green",
    "blue",
    "orange",
    "aquamarine",
    "deeppink",
    "darkviolet",
    "darkolivegreen",
    "cyan"
  };

  struct Vertex {
    tcParticleSet in;
    tcParticleSet out;
  };

  typedef map<long, Vertex> VertexMap;

  template<typename Iter>
  Lorentz5Momentum sumP(Iter a, Iter b) {
    Lorentz5Momentum sum;
    for ( Iter it = a; it != b; ++it ) {
      sum += (*it)->momentum();
    }
    return sum;
  }
}

void ThePEG::Event::printGraphviz() const {
  ThePEG::printGraphviz(cout, this);
}


void ThePEG::printGraphviz(ostream & os, tcEventPtr ev) {
  os << header
     << "node [width=0.03,height=0.03,shape=point,label=\"\"];\n";

  tcParticleSet all;
  ev->select(inserter(all), SelectAll());

  VertexMap vertices;
  for (tcParticleSet::const_iterator it = all.begin();
       it != all.end(); ++it) {
    tcPPtr p = (*it);

    long start = startnode(p);
    long end = endnode(p);

    vertices[start].out.insert(p);
    vertices[end  ].in .insert(p);

    os << start << " -> " << end
       << " [label=\"" << p->number() << ' '
       << p->PDGName() << "\\n"
       << p->momentum().e()/GeV << "\\n"
       << p->momentum().mass()/GeV
       << "\"";

    if ( p->hasColourInfo() &&
	 ( p->colourLine() || p->antiColourLine() )) {
      os << ",penwidth=2,color=\"";

      const vector<tcColinePtr> & clines = p->colourInfo()->colourLines();
      for ( int i = 0, N = clines.size(); i < N; ++i ) {
	int colindex = ev->colourLineIndex(clines[i]) % 9;
	if ( i > 0 )  os << ':';
	os << colours[colindex];
      }
      const vector<tcColinePtr> & aclines = p->colourInfo()->antiColourLines();
      for ( int i = 0, N = aclines.size(); i < N; ++i ) {
	int colindex = ev->colourLineIndex(aclines[i]) % 9;
	if ( i > 0 || !clines.empty() )  os << ':';
	os << colours[colindex];
      }
      os << '"';
    }
    os << "];\n";
  }

  int label = 0;
  for ( VertexMap::const_iterator v = vertices.begin();
	v != vertices.end(); ++v ) {

    const long vertexId = v->first;
    const tcParticleSet & in  = v->second.in;
    const tcParticleSet & out = v->second.out;

    if ( in.empty() || out.empty() ) 
      continue;

    Lorentz5Momentum diff 
      = sumP(out.begin(), out.end()) - sumP(in.begin(), in.end());

    if ( abs(diff.e()) > 1.0*GeV ) {
      ++label;
      std::stringstream tail;
      tail << " [label=\"" << std::setprecision(4) 
	   << abs(diff.e()/GeV)
	   << " GeV\",arrowsize=3,penwidth=5,"
	   << "color=\"#ff000010\",fontcolor=\"#ff000010\"];\n";
      if ( diff.e() > ZERO )
	os << "mom" << label << " -> " << vertexId << tail.str();
      else
	os << vertexId << " -> " << "mom" << label << tail.str();
    }
  }

  os << '}' << endl;
}

void Event::debugme() const {
  cerr << *this;
  EventRecordBase::debugme();
}

void Event::persistentOutput(PersistentOStream & os) const {
  os << theIncoming << theCollisions << allSteps << allSubProcesses
     << allParticles << theNumber << theWeight << theOptionalWeights 
     << theParticleNumber;
  EventConfig::putHandler(os, theHandler);
  const auto q=laboratoryRotation(theLaboratoryBoost);
  // Lambda = B(beta) R. The signed quaternion for R preserves the spin-1/2
  // sign as well as the ordinary spatial rotation.
  os << theHasGeometry << theGeometryComplete
     << theLaboratoryBoost.one().xt()/theLaboratoryBoost.one().tt()
     << theLaboratoryBoost.one().yt()/theLaboratoryBoost.one().tt()
     << theLaboratoryBoost.one().zt()/theLaboratoryBoost.one().tt()
     << q[0] << q[1] << q[2] << q[3]
     << ounit(theInteractionPoint,mm);
}

void Event::persistentInput(PersistentIStream & is, int version) {
  is >> theIncoming >> theCollisions >> allSteps >> allSubProcesses
     >> allParticles >> theNumber >> theWeight >> theOptionalWeights
     >> theParticleNumber;
  EventConfig::getHandler(is, theHandler);
  if (version >= 1) {
    double bx,by,bz;
    is >> theHasGeometry >> theGeometryComplete >> bx >> by >> bz;
    if (version>=2) {
      double qw,qx,qy,qz;
      is >> qw >> qx >> qy >> qz;
      theLaboratoryBoost=laboratoryTransformation(bx,by,bz,qw,qx,qy,qz);
    } else theLaboratoryBoost=LorentzRotation(bx,by,bz);
    is >> iunit(theInteractionPoint,mm);
  }
}

ClassDescription<Event> Event::initEvent;

void Event::Init() {}

