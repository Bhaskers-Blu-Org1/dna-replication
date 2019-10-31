#ifndef DNA_REPLICATION_ORIGIN_H
#define DNA_REPLICATION_ORIGIN_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "Chromosome.h"
#include "misc/Vector3.h"

namespace DNAReplication {

    // forward declarations
    class Particle;

    /**
     * Replication origin
     */
    class Origin {

    public:

        /// Replication state of an origin
        enum class State : unsigned char {
            Pre, ///< Pre-replicative, origin can be fired
            Pass, ///< Passively replicated, origin will never be fired
            ReplLR, ///< Replicative, in both directions along the chromosome
            ReplL, ///< Replicative, towards the left along the chromosome
            ReplR, ///< Replicative, towards the right along the chromosome
            Post ///< Post-replicative, origin finished replication
        };

        struct OriginData {

            friend class Origin;

        public:

            std::string id;
            std::string chromosomeID;
            unsigned long int pos = 0;

            OriginData(std::string const &id, std::string const &chromosomeID, unsigned long int pos)
                    : id(id), chromosomeID(chromosomeID), pos(pos) {
            }

        private:

            OriginData() = default;

        };

    private:

        OriginData const data;

        State state = State::Pre;
        double firingTime = 0;
        Particle *boundParticle = nullptr;
        Chromosome const *chromosome = nullptr;
        Chromosome::Granule const *chromosomeGranule = nullptr;

        // direct neighbor origins on the chromosome
        Origin *leftOrigin = nullptr;
        Origin *rightOrigin = nullptr;

        // next replicating origins on th chromosome
        Origin *leftReplOrigin = nullptr;
        Origin *rightReplOrigin = nullptr;

        // origins that will be passively replicated next
        Origin *nextLeftPassOrigin = nullptr;
        Origin *nextRightPassOrigin = nullptr;

        /// Finds the closest (in terms of base-pair distance) replicating origin on the left
        Origin *findLeftReplOrigin() const;

        /// Finds the closest (in terms of base-pair distance) replicating origin on the right
        Origin *findRightReplOrigin() const;

        /// Computes the time of collision of the left replication fork with a right replication fork of a replicating origin on the left or with the contig start
        /// \param vFork The replication fork velocity (in bp/s)
        /// \return The time of collision / contig start arrival
        double getLeftCollisionTime(double vFork) const;

        /// Computes the time of collision of the right replication fork with a left replication fork of a replicating origin on the right or with the contig end
        /// \param vFork The replication fork velocity (in bp/s)
        /// \return The time of collision / contig end arrival
        double getRightCollisionTime(double vFork) const;

    public:

        /// Constructs a new origin
        explicit Origin(OriginData data);

        /// Return the ID (i.e., the name) of the origin
        std::string const &getID() const;

        /// Two origins are considered equal if their ids match
        bool operator==(Origin const &other) const;

        /// Returns the state of the origin
        State getState() const;

        /// Returns the firing time of the origin
        double getFiringTime() const;

        /// Returns the particle bound to this origin
        Particle *getBoundParticle() const;

        /// Sets the particle bound to this origin
        void setBoundParticle(Particle *boundParticle);

        /// Returns the chromosome of this origin
        Chromosome const *getChromosome() const;

        /// Returns the granule (3D position) of the origin
        Chromosome::Granule const *getChromosomeGranule() const;

        /// Initializes the chromosome (and granule) of this origin
        void initializeChromosome(Chromosome const *chromosome);

        /// Links the direct neighbor origins
        void linkNeighborOrigins(std::vector<Origin> &origins);

        /// Returns the time of the first collision (left or right)
        double getMinCollisionTime(double vFork) const;

        /// Returns the time of the last collision (left or right)
        double getMaxCollisionTime(double vFork) const;

        /// Fires the origin
        /// \param tFire The time of firing (in seconds)
        void fire(double tFire);

        /// Updates ("moves") the position of the left replication fork
        /// \param tCurrent The current time (in seconds)
        /// \param vFork The replication fork velocity (in bp/s)
        /// \param releasePos Out parameter, not used (imer assumption -> particles are released "to the right")
        /// \return The number of passively replicated origins
        int replicateLeft(double tCurrent, double vFork, Vector3<double> const **releasePos);

        /// Updates ("moves") the position of the right replication fork
        /// \param tCurrent The current time (in seconds)
        /// \param vFork The replication fork velocity (in bp/s)
        /// \param releasePos Out parameter; release position of the bound particle, if replication forks collided
        /// \return The number of passively replicated origins
        int replicateRight(double tCurrent, double vFork, Vector3<double> const **releasePos);

        Origin(Origin const &copy) {
            throw std::runtime_error("Internal error: copying origins is not allowed");
        }

        Origin &operator=(Origin const &copy) {
            throw std::runtime_error("Internal error: copying origins is not allowed");
        }

    };

}

#endif //DNA_REPLICATION_ORIGIN_H
