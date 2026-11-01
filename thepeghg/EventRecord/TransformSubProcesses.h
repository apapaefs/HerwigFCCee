// Internal helpers for Event and Collision transformations.
#ifndef ThePEG_TransformSubProcesses_H
#define ThePEG_TransformSubProcesses_H

#include "ThePEG/EventRecord/SubProcessGroup.h"
#include "ThePEG/EventRecord/ParticleTraits.h"

namespace ThePEG {
namespace EventRecordDetail {

inline void collectDetachedParticles(tSubProPtr sub,
                                     const ParticleSet & registered,
                                     ParticleSet & detached,
                                     SubProcessSet & visited) {
  if (!sub || !visited.insert(sub).second) return;
  const auto add = [&](tPPtr p) {
    if (p && registered.find(p) == registered.end()) detached.insert(p);
  };
  add(sub->incoming().first);
  add(sub->incoming().second);
  for (const auto & p : sub->intermediates()) add(p);
  for (const auto & p : sub->outgoing()) add(p);
  const auto group = dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(sub);
  if (group)
    for (const auto & dependent : group->dependent())
      collectDetachedParticles(dependent, registered, detached, visited);
}

template <typename Iterator>
ParticleSet allParticles(const ParticleSet & registered, Iterator first,
                         Iterator last) {
  // Projectors can share particles with the head and with other dependents.
  // Do not register counterevents in the physical event history, and do not
  // invoke SubProcessGroup::transform after transforming registered particles.
  ParticleSet detached;
  SubProcessSet visited;
  for (; first != last; ++first)
    collectDetachedParticles(*first, registered, detached, visited);
  detached.insert(registered.begin(),registered.end());
  return detached;
}

template <typename Iterator>
void transformParticles(const ParticleSet & registered, Iterator first,
                        Iterator last, const LorentzRotation & rotation) {
  for (const auto & p : allParticles(registered,first,last)) p->transform(rotation);
}

}
}
#endif
