#include "Simulation.h"

namespace DNAReplication {

    Simulation::Simulation(double vFork, std::vector<Origin::OriginData> const &originData, std::vector<Chromosome::ChromosomeData> const &chromosomeData,
                           ParticleDiffusionBehavior &diffusionBehavior, ParticleActivationBehavior &activationBehavior, ParticleBindingBehavior &bindingBehavior)
            : vFork(vFork), diffusionBehavior(diffusionBehavior), activationBehavior(activationBehavior), bindingBehavior(bindingBehavior) {
        // create chromosomes
        std::size_t nextChromosomeIndex = 0;
        chromosomes.reserve(chromosomeData.size());
        std::map<std::string, std::size_t> chromosomeIndices;
        for (Chromosome::ChromosomeData const &data : chromosomeData) {
            chromosomes.emplace_back(data);
            chromosomeIndices[data.id] = nextChromosomeIndex++;
        }
        // create origins
        origins.reserve(originData.size());
        for (Origin::OriginData const &data : originData) {
            origins.emplace_back(data);
            std::size_t chromosomeIndex = chromosomeIndices.at(data.chromosomeID);
            origins.back().initializeChromosome(&chromosomes[chromosomeIndex]);
        }
        // link origins
        for (Origin &origin : origins) {
            origin.linkNeighborOrigins(origins);
        }
    }

    void Simulation::initializeParticles(unsigned short int nParticles) {
        particles.reserve(nParticles);
        for (unsigned short int i = 0; i < nParticles; ++i) {
            particles.emplace_back(true, diffusionBehavior.getRandomPosition());
            particles.back().setActive(activationBehavior.isActiveInitially(particles.back()));
        }
    }

    double Simulation::getCurrentTime() const {
        return tCurrent;
    }

    std::vector<Origin> const &Simulation::getOrigins() const {
        return origins;
    }

    std::vector<Particle> const &Simulation::getParticles() const {
        return particles;
    }

    void Simulation::clearObservers() {
        observers.clear();
    }

    void Simulation::registerObserver(std::shared_ptr<SimulationObserver> observer) {
        observers.push_back(std::move(observer));
    }

    void Simulation::run() {
        // initialize caches
        std::vector<Origin *> replicatingOrigins;
        replicatingOrigins.reserve(origins.size());
        std::vector<Particle *> mobileParticles;
        std::vector<Particle *> newMobileParticles;
        mobileParticles.reserve(particles.size());
        newMobileParticles.reserve(particles.size());
        for (Particle &particle : particles) {
            mobileParticles.push_back(&particle);
        }
        notifySimulationStarted();
        // simulate diffusion as long as pre-replicative origins are available
        std::size_t numPreOrigins = origins.size();
        while (numPreOrigins > 0) {
            // time step
            tCurrent += diffusionBehavior.timeStep();
            // optimization for the common case that there is only one mobile particle
            if (mobileParticles.size() == 1) {
                // diffuse and activate particle
                Particle &mobileParticle = *mobileParticles[0];
                diffuseParticle(mobileParticle);
                activateParticle(mobileParticle);
                // bind and fire origin
                if (mobileParticle.isActive()) {
                    std::vector<Origin *> preOriginsInProximity = updateParticle(mobileParticle);
                    bindingBehavior.shuffleOrigins(preOriginsInProximity);
                    for (Origin *preOriginInProximity : preOriginsInProximity) {
                        bool const boundToOrigin = bindParticleToOrigin(mobileParticle, *preOriginInProximity);
                        if (boundToOrigin) {
                            mobileParticles.clear();
                            fireOrigin(*preOriginInProximity);
                            replicatingOrigins.push_back(preOriginInProximity);
                            numPreOrigins--;
                            break;
                        }
                    }
                }
            } else {
                // diffuse and activate particles
                for (Particle *mobileParticle : mobileParticles) {
                    diffuseParticle(*mobileParticle);
                    activateParticle(*mobileParticle);
                }
                // bind and fire origins
                newMobileParticles.clear();
                bindingBehavior.shuffleParticles(mobileParticles);
                for (Particle *mobileParticle : mobileParticles) {
                    if (mobileParticle->isActive()) {
                        std::vector<Origin *> preOriginsInProximity = updateParticle(*mobileParticle);
                        bindingBehavior.shuffleOrigins(preOriginsInProximity);
                        for (Origin *preOriginInProximity : preOriginsInProximity) {
                            // check if origin is still unbound (other particles might have bound to it in the meantime)
                            if (preOriginInProximity->getBoundParticle() == nullptr) {
                                bool const boundToOrigin = bindParticleToOrigin(*mobileParticle, *preOriginInProximity);
                                if (boundToOrigin) {
                                    fireOrigin(*preOriginInProximity);
                                    replicatingOrigins.push_back(preOriginInProximity);
                                    numPreOrigins--;
                                    break;
                                }
                            }
                        }
                    }
                    if (mobileParticle->getBoundOrigin() == nullptr) {
                        newMobileParticles.push_back(mobileParticle);
                    }
                }
                std::swap(mobileParticles, newMobileParticles);
            }
            // optimization: skip diffusion dynamics until a particle is released
            if (mobileParticles.empty()) {
                double nextReleaseTime = std::numeric_limits<double>::max();
                for (Origin const *replicatingOrigin : replicatingOrigins) {
                    double const minCollisionTime = replicatingOrigin->getMinCollisionTime(vFork);
                    if (minCollisionTime < nextReleaseTime) {
                        nextReleaseTime = minCollisionTime;
                    }
                }
                tCurrent = nextReleaseTime;
            }
            // replicate origins
            for (auto replicatingOriginIter = replicatingOrigins.begin(); replicatingOriginIter != replicatingOrigins.end();) {
                Origin &replicatingOrigin = **replicatingOriginIter;
                Vector3<double> const *particleReleasePos = nullptr;
                numPreOrigins -= replicateOrigin(replicatingOrigin, &particleReleasePos);
                // release bound particle
                if (particleReleasePos != nullptr) {
                    mobileParticles.push_back(replicatingOrigin.getBoundParticle());
                    releaseParticleFromOrigin(replicatingOrigin, *particleReleasePos);
                }
                // update replicating origins cache
                if (replicatingOrigin.getState() == Origin::State::Post) {
                    replicatingOriginIter = replicatingOrigins.erase(replicatingOriginIter);
                } else {
                    ++replicatingOriginIter; 
                }
            }
            notifyIterationCompleted();
        }
        // finish "replication-only" phase (without simulating diffusion dynamics)
        if (!replicatingOrigins.empty()) {
            // sort remaining replicating origins by their completion times
            std::sort(replicatingOrigins.begin(), replicatingOrigins.end(), [this](Origin *o1, Origin *o2) {
                return o1->getMaxCollisionTime(vFork) < o2->getMaxCollisionTime(vFork);
            });
            // simulate completion of the remaining replicating origins
            for (Origin *replicatingOrigin : replicatingOrigins) {
                tCurrent = replicatingOrigin->getMaxCollisionTime(vFork);
                Vector3<double> const *particleReleasePos = nullptr;
                replicateOrigin(*replicatingOrigin, &particleReleasePos);
                notifyIterationCompleted();
            }
        }
    }

    void Simulation::diffuseParticle(Particle &particle) {
        // diffuse particle
        Vector3<double> newPos = diffusionBehavior.diffuse(particle.getPos());
        if (!diffusionBehavior.inDomain(newPos)) {
            newPos = diffusionBehavior.reflect(newPos);
        }
        particle.setPos(newPos);
        notifyParticleDiffused(particle);
    }

    void Simulation::activateParticle(Particle &particle) {
        if (!particle.isActive() && activationBehavior.checkSPBActivation(particle) && diffusionBehavior.inSPB(particle.getPos())) {
            particle.setActive(true);
            notifyParticleActivationStateChanged(particle);
        } else if (particle.isActive() && activationBehavior.checkPeripheryInactivation(particle) && diffusionBehavior.inPeriphery(particle.getPos())) {
            particle.setActive(false);
            particle.clearPreOriginsInProximity();
            notifyParticleActivationStateChanged(particle);
        }
    }

    std::vector<Origin *> const &Simulation::updateParticle(Particle &particle) {
        particle.clearPreOriginsInProximity();
        for (Origin &origin : origins) {
            if (origin.getState() == Origin::State::Pre && bindingBehavior.inProximity(particle, origin)) {
                particle.addPreOriginInProximity(&origin);
            }
        }
        return particle.getPreOriginsInProximity();
    }

    bool Simulation::bindParticleToOrigin(Particle &particle, Origin &origin) {
        std::vector<Origin *> const &previous = particle.getPreviousPreOriginsInProximity();
        if (std::find(previous.begin(), previous.end(), &origin) == previous.end() && bindingBehavior.checkBinding(particle, origin)) {
            particle.setBoundOrigin(&origin);
            origin.setBoundParticle(&particle);
            notifyParticleBindingStateChanged(particle);
            return true;
        }
        return false;
    }

    void Simulation::fireOrigin(Origin &origin) {
        origin.fire(tCurrent);
        notifyOriginFired(origin);
    }

    int Simulation::replicateOrigin(Origin &origin, Vector3<double> const **particleReleasePos) {
        int nPassivelyActivated = 0;
        if (origin.getState() == Origin::State::ReplLR || origin.getState() == Origin::State::ReplL) {
            nPassivelyActivated += origin.replicateLeft(tCurrent, vFork, particleReleasePos);
        }
        if (origin.getState() == Origin::State::ReplLR || origin.getState() == Origin::State::ReplR) {
            nPassivelyActivated += origin.replicateRight(tCurrent, vFork, particleReleasePos);
        }
        notifyOriginReplicated(origin);
        return nPassivelyActivated;
    }

    void Simulation::releaseParticleFromOrigin(Origin &origin, Vector3<double> const &particleReleasePos) {
        Particle &particle = *origin.getBoundParticle();
        origin.setBoundParticle(nullptr);
        particle.setBoundOrigin(nullptr);
        particle.setPos(particleReleasePos);
        notifyParticleBindingStateChanged(particle);
    }

    void Simulation::notifySimulationStarted() const {
        for (std::shared_ptr<SimulationObserver> const &o : observers) {
            o->handleSimulationStarted({*this});
        }
    }

    void Simulation::notifyIterationCompleted() const {
        for (std::shared_ptr<SimulationObserver> const &o : observers) {
            o->handleIterationCompleted({*this});
        }
    }

    void Simulation::notifyParticleDiffused(Particle const &particle) const {
        for (std::shared_ptr<SimulationObserver> const &o : observers) {
            o->handleParticleDiffused({*this, particle});
        }
    }

    void Simulation::notifyParticleActivationStateChanged(Particle const &particle) const {
        for (std::shared_ptr<SimulationObserver> const &o : observers) {
            o->handleParticleActivationStateChanged({*this, particle});
        }
    }

    void Simulation::notifyParticleBindingStateChanged(Particle const &particle) const {
        for (std::shared_ptr<SimulationObserver> const &o : observers) {
            o->handleParticleBindingStateChanged({*this, particle});
        }
    }

    void Simulation::notifyOriginFired(Origin const &origin) const {
        for (std::shared_ptr<SimulationObserver> const &o : observers) {
            o->handleOriginFired({*this, origin});
        }
    }

    void Simulation::notifyOriginReplicated(Origin const &origin) const {
        for (std::shared_ptr<SimulationObserver> const &o : observers) {
            o->handleOriginReplicated({*this, origin});
        }
    }

}
