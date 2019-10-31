#ifndef DNA_REPLICATION_PROBABILISTICPARTICLEBINDINGBEHAVIOR_H
#define DNA_REPLICATION_PROBABILISTICPARTICLEBINDINGBEHAVIOR_H

#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

#include "../ParticleBindingBehavior.h"
#include "../misc/Vector3.h"

namespace DNAReplication {

    /**
     * Probabilistic binding of particles to replication origins
     */
    class ProbabilisticParticleBindingBehavior : public ParticleBindingBehavior {

    private:

        double const dBind;
        double const pBind;
        boost::random::mt19937 rng;

    public:

        /// Constructs a new instance
        /// \param pBind The constant binding probability
        /// \param dBind The maximal distance between the particle and the replication origin (per vector component, in um)
        ProbabilisticParticleBindingBehavior(double dBind, double pBind);

        void shuffleOrigins(std::vector<Origin *> &origins) override;

        void shuffleParticles(std::vector<Particle *> &particles) override;

        bool inProximity(Particle const &p, Origin const &o) const override;

        bool checkBinding(Particle const &p, Origin const &o) override;

    };

}

#endif //DNA_REPLICATION_PROBABILISTICPARTICLEBINDINGBEHAVIOR_H
