#ifndef DNA_REPLICATION_PROBABILISTICPARTICLEACTIVATIONBEHAVIOR_H
#define DNA_REPLICATION_PROBABILISTICPARTICLEACTIVATIONBEHAVIOR_H

#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

#include "../ParticleActivationBehavior.h"

namespace DNAReplication {

    /**
     * Probabilistic SPB-mediated activation and peripheral inactivation of particles
     */
    class ProbabilisticParticleActivationBehavior : public ParticleActivationBehavior {

    private:

        const double pAct;
        const bool spbActivationEnabled;
        const bool peripheryInactivationEnabled;
        boost::random::mt19937 rng;

    public:

        ProbabilisticParticleActivationBehavior(double pAct, bool spbActivationEnabled, bool peripheryInactivationEnabled);

        bool isActiveInitially(Particle const &particle) override;

        bool checkSPBActivation(Particle const &particle) override;

        bool checkPeripheryInactivation(Particle const &particle) override;

    };

}


#endif //DNA_REPLICATION_PROBABILISTICPARTICLEACTIVATIONBEHAVIOR_H
