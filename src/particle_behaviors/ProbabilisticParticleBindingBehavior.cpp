#include "ProbabilisticParticleBindingBehavior.h"

namespace DNAReplication {

    ProbabilisticParticleBindingBehavior::ProbabilisticParticleBindingBehavior(double dBind, double pBind)
            : dBind(dBind), pBind(pBind), rng(boost::random::random_device{}()) {
    }

    void ProbabilisticParticleBindingBehavior::shuffleOrigins(std::vector<Origin *> &origins) {
        std::shuffle(origins.begin(), origins.end(), rng);
    }

    void ProbabilisticParticleBindingBehavior::shuffleParticles(std::vector<Particle *> &particles) {
        std::shuffle(particles.begin(), particles.end(), rng);
    }

    bool ProbabilisticParticleBindingBehavior::inProximity(Particle const &p, Origin const &o) const {
        Vector3<double> const diff = p.getPos() - o.getChromosomeGranule()->pos;
        return std::abs(diff.x) <= dBind && std::abs(diff.y) <= dBind && std::abs(diff.z) <= dBind;
    }

    bool ProbabilisticParticleBindingBehavior::checkBinding(Particle const &p, Origin const &o) {
        boost::random::bernoulli_distribution<double> bindDist(pBind);
        return bindDist(rng);
    }

}
