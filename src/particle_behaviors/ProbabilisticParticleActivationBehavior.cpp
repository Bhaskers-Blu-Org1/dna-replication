#include "ProbabilisticParticleActivationBehavior.h"

namespace DNAReplication {

    ProbabilisticParticleActivationBehavior::ProbabilisticParticleActivationBehavior(double pAct, bool spbActivationEnabled, bool peripheryInactivationEnabled)
            : pAct(pAct), spbActivationEnabled(spbActivationEnabled), peripheryInactivationEnabled(peripheryInactivationEnabled), rng(boost::random::random_device{}()) {
    }

    bool ProbabilisticParticleActivationBehavior::isActiveInitially(Particle const &particle) {
        return !spbActivationEnabled;
    }

    bool ProbabilisticParticleActivationBehavior::checkSPBActivation(Particle const &particle) {
        boost::random::bernoulli_distribution<double> activateDist(pAct);
        return activateDist(rng);
    }

    bool ProbabilisticParticleActivationBehavior::checkPeripheryInactivation(Particle const &particle) {
        return peripheryInactivationEnabled;
    }

}