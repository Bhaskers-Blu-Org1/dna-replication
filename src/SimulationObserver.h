#ifndef DNA_REPLICATION_SIMULATIONOBSERVER_H
#define DNA_REPLICATION_SIMULATIONOBSERVER_H

#include "Origin.h"
#include "Particle.h"
#include "Simulation.h"

namespace DNAReplication {

    // forward declarations
    class Simulation;

    /**
     * Interface for observing simulation events
     */
    class SimulationObserver {

    public:

        struct SimulationEvent {
            Simulation const &simulation;
        };

        struct SimulationOriginEvent {
            Simulation const &simulation;
            Origin const &origin;
        };

        struct SimulationParticleEvent {
            Simulation const &simulation;
            Particle const &particle;
        };

        virtual void handleSimulationStarted(SimulationEvent const &e) {};

        virtual void handleIterationCompleted(SimulationEvent const &e) {};

        virtual void handleParticleDiffused(SimulationParticleEvent const &e) {};

        virtual void handleParticleActivationStateChanged(SimulationParticleEvent const &e) {};

        virtual void handleParticleBindingStateChanged(SimulationParticleEvent const &e) {};

        virtual void handleOriginFired(SimulationOriginEvent const &e) {};

        virtual void handleOriginReplicated(SimulationOriginEvent const &e) {};

    };

}

#endif //DNA_REPLICATION_SIMULATIONOBSERVER_H
