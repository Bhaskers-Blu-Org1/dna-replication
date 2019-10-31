#ifndef DNA_REPLICATION_PARTICLEBINDINGBEHAVIOR_H
#define DNA_REPLICATION_PARTICLEBINDINGBEHAVIOR_H

#include "Origin.h"
#include "Particle.h"

namespace DNAReplication {

    /**
     * Abstract class modeling the binding behavior of particles
     */
    class ParticleBindingBehavior {

    public:

        /// Randomly shuffle origins
        virtual void shuffleOrigins(std::vector<Origin *> &origins) = 0;

        /// Randomly shuffle particles
        virtual void shuffleParticles(std::vector<Particle *> &particles) = 0;

        /// Checks whether a particle is in the binding proximity of an origin
        virtual bool inProximity(Particle const &p, Origin const &o) const = 0;

        /// Checks whether a particle should bind to an origin
        virtual bool checkBinding(Particle const &p, Origin const &o) = 0;

    };

}

#endif //DNA_REPLICATION_PARTICLEBINDINGBEHAVIOR_H
