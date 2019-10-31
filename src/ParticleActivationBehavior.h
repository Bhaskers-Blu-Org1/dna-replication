#ifndef DNA_REPLICATION_PARTICLEACTIVATIONBEHAVIOR_H
#define DNA_REPLICATION_PARTICLEACTIVATIONBEHAVIOR_H

#include "Particle.h"

namespace DNAReplication {

    /**
     * Abstract class modeling the activation behavior of particles
     */
    class ParticleActivationBehavior {

    public:

        /// Determines whether the particle is active at the start of a simulation
        virtual bool isActiveInitially(Particle const &particle) = 0;

        /// Determines whether the particle should be activated by the SPB in the current iteration
        virtual bool checkSPBActivation(Particle const &particle) = 0;

        /// Determines whether the particle should be inactivated by the periphery in the current iteration
        virtual bool checkPeripheryInactivation(Particle const &particle) = 0;

    };

}

#endif //DNA_REPLICATION_PARTICLEACTIVATIONBEHAVIOR_H
