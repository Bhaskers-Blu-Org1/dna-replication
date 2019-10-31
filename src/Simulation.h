#ifndef DNA_REPLICATION_SIMULATION_H
#define DNA_REPLICATION_SIMULATION_H

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <vector>

#include "Chromosome.h"
#include "Origin.h"
#include "Particle.h"
#include "ParticleActivationBehavior.h"
#include "ParticleBindingBehavior.h"
#include "ParticleDiffusionBehavior.h"
#include "SimulationObserver.h"
#include "misc/Vector3.h"

namespace DNAReplication {

    // forward declarations
    class SimulationObserver;

    /**
     * DNA replication simulation
     */
    class Simulation {

    private:

        double const vFork;

        double tCurrent = 0;
        std::vector<Origin> origins;
        std::vector<Particle> particles;
        std::vector<Chromosome> chromosomes;
        ParticleDiffusionBehavior &diffusionBehavior;
        ParticleActivationBehavior &activationBehavior;
        ParticleBindingBehavior &bindingBehavior;
        std::vector<std::shared_ptr<SimulationObserver>> observers;

        /// "Moves" a specified particle within the nucleus domain
        void diffuseParticle(Particle &particle);

        /// Activates or deactivates a particle (after diffusion step)
        void activateParticle(Particle &particle);

        /// Updates a particle's list of origins in proximity
        std::vector<Origin *> const &updateParticle(Particle &particle);

        /// Binds a particle to the specified (pre-replicative) origin
        /// \return True if the binding was successful, false otherwise
        bool bindParticleToOrigin(Particle &particle, Origin &origin);

        /// Activates a (bound) origin
        void fireOrigin(Origin &origin);

        /// Replicates an (activated) origin
        /// \return The number of passively activated origins
        int replicateOrigin(Origin &origin, Vector3<double> const **particleReleasePos);

        /// Releases the particle bound to an origin at the specified position
        void releaseParticleFromOrigin(Origin &origin, Vector3<double> const &particleReleasePos);

        /// Notifies observers about the start of the simulation
        void notifySimulationStarted() const;

        /// Notifies observers about the completion of an iteration
        void notifyIterationCompleted() const;

        /// Notifies observers about the diffusion of a particle
        void notifyParticleDiffused(Particle const &particle) const;

        /// Notifies observers about changes in the activation state of a particle
        void notifyParticleActivationStateChanged(Particle const &particle) const;

        /// Notifies observers about changes in the binding state of a particle
        void notifyParticleBindingStateChanged(Particle const &particle) const;

        /// Notifies observers about the firing of an origin
        void notifyOriginFired(Origin const &origin) const;

        /// Notifies observers about the replication of an origin
        void notifyOriginReplicated(Origin const &origin) const;

    public:

        /// Constructs a new simulation
        /// \param vFork The replication fork velocity (in bp/s)
        Simulation(double vFork, std::vector<Origin::OriginData> const &originData, std::vector<Chromosome::ChromosomeData> const &chromosomeData,
                   ParticleDiffusionBehavior &diffusionBehavior, ParticleActivationBehavior &activationBehavior, ParticleBindingBehavior &bindingBehavior);

        /// Randomly initializes particles within the entire nucleus domain
        void initializeParticles(unsigned short int nParticles);

        /// Returns the current "experimental" (simulated) time in seconds
        double getCurrentTime() const;

        /// Returns a read-only list of replication origins
        std::vector<Origin> const &getOrigins() const;

        /// Returns a read-only list of particles (activation factors)
        std::vector<Particle> const &getParticles() const;

        /// Clears all observers
        void clearObservers();

        /// Registers an observer
        void registerObserver(std::shared_ptr<SimulationObserver> observer);

        /// Executes the simulation
        void run();

    };

}

#endif //DNA_REPLICATION_SIMULATION_H
