#ifndef DNA_REPLICATION_PARTICLEDIFFUSIONBEHAVIOR_H
#define DNA_REPLICATION_PARTICLEDIFFUSIONBEHAVIOR_H

#include "misc/Vector3.h"

namespace DNAReplication {

    /**
     * Abstract class modeling the diffusion behavior of particles
     */
    class ParticleDiffusionBehavior {

    public:

        /// Returns the next time step ("delta_t") in seconds
        virtual double timeStep() = 0;

        /// Returns a random 3D position within the nucleus domain
        virtual Vector3<double> getRandomPosition() = 0;

        /// Checks whether the specified position is within the nucleus domain
        virtual bool inDomain(Vector3<double> const &pos) const = 0;

        /// Checks whether the specified position is within the spindle pole body (SPB) region
        virtual bool inSPB(Vector3<double> const &pos) const = 0;

        /// Checks whether the specified position is within the peripheral region
        virtual bool inPeriphery(Vector3<double> const &pos) const = 0;

        /// Returns the next position (for particle diffusion)
        virtual Vector3<double> diffuse(Vector3<double> const &pos) = 0;

        /// Reflects the specified invalid position into the interior of the nucleus
        virtual Vector3<double> reflect(Vector3<double> const &pos) = 0;

    };

}

#endif //DNA_REPLICATION_PARTICLEDIFFUSIONBEHAVIOR_H
