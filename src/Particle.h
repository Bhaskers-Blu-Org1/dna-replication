#ifndef DNA_REPLICATION_PARTICLE_H
#define DNA_REPLICATION_PARTICLE_H

#include <stdexcept>
#include <vector>

#include "misc/Vector3.h"

namespace DNAReplication {

    // forward declarations
    class Origin;

    /**
     * Activation factor ("particle")
     */
    class Particle {

    private:

        bool active = false;
        Vector3<double> pos{};
        Origin *boundOrigin = nullptr;
        std::vector<Origin *> preOriginsInProximity{};
        std::vector<Origin *> previousPreOriginsInProximity{};

    public:

        /// Constructs a new particle
        Particle(bool active, Vector3<double> const &pos);

        /// Determines whether the particle is currently active (i.e., ready to fire origins)
        bool isActive() const;

        /// Activates or inactivates the particle
        void setActive(bool active);

        /// Returns the current 3D position of the particle
        Vector3<double> const &getPos() const;

        /// Sets the current 3D position of the particle
        void setPos(Vector3<double> const &pos);

        // Returns the origin this particle is currently bound to
        Origin *getBoundOrigin() const;

        /// Sets the origin this particle is currently bound to
        void setBoundOrigin(Origin *boundOrigin);

        /// Clears the list of pre-replicative origins in proximity
        void clearPreOriginsInProximity();

        /// Adds an origin to the list of pre-replicative origins in proximity
        void addPreOriginInProximity(Origin *origin);

        /// Returns a read-only list of the pre-replicative origins in proximity
        std::vector<Origin *> const &getPreOriginsInProximity() const;

        /// Returns a read-only list of the old pre-replicative origins in proximity
        std::vector<Origin *> const &getPreviousPreOriginsInProximity() const;

        // Disallow copying

        Particle(Particle const &copy) {
            throw std::runtime_error("Internal error: copying particles is not allowed");
        }

        Particle &operator=(Particle const &copy) {
            throw std::runtime_error("Internal error: copying particles is not allowed");
        }

    };

}

#endif //DNA_REPLICATION_PARTICLE_H
