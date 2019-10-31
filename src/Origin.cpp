#include "Origin.h"

namespace DNAReplication {

    Origin::Origin(Origin::OriginData data) : data(std::move(data)) {
    }

    std::string const &Origin::getID() const {
        return data.id;
    }

    bool Origin::operator==(Origin const &other) const {
        return data.id == other.data.id;
    }

    Origin::State Origin::getState() const {
        return state;
    }

    double Origin::getFiringTime() const {
        return firingTime;
    }

    Particle *Origin::getBoundParticle() const {
        return boundParticle;
    }

    void Origin::setBoundParticle(Particle *boundParticle) {
        this->boundParticle = boundParticle;
    }

    Chromosome const *Origin::getChromosome() const {
        return chromosome;
    }

    Chromosome::Granule const *Origin::getChromosomeGranule() const {
        return chromosomeGranule;
    }

    void Origin::initializeChromosome(Chromosome const *chromosome) {
        this->chromosome = chromosome;
        this->chromosomeGranule = (chromosome != nullptr) ? &chromosome->findGranule(data.pos) : nullptr;
    }

    void Origin::linkNeighborOrigins(std::vector<Origin> &origins) {
        for (Origin &o : origins) {
            bool const closerThanLeft = (leftOrigin == nullptr || o.data.pos > leftOrigin->data.pos);
            if (o.data.chromosomeID == data.chromosomeID && chromosome->inSameContig(o.data.pos, data.pos) && o.data.pos < data.pos && closerThanLeft) {
                leftOrigin = &o;
                nextLeftPassOrigin = &o;
            }
        }
        for (Origin &o : origins) {
            bool const closerThanRight = (rightOrigin == nullptr || o.data.pos < rightOrigin->data.pos);
            if (o.data.chromosomeID == data.chromosomeID && chromosome->inSameContig(o.data.pos, data.pos) &&
                o.data.pos > data.pos && closerThanRight) {
                rightOrigin = &o;
                nextRightPassOrigin = &o;
            }
        }
    }

    double Origin::getMinCollisionTime(double vFork) const {
        if (state == State::ReplL) {
            return getLeftCollisionTime(vFork);
        }
        if (state == State::ReplR) {
            return getRightCollisionTime(vFork);
        }
        if (state == State::ReplLR) {
            return std::min(getLeftCollisionTime(vFork), getRightCollisionTime(vFork));
        }
        throw std::runtime_error("Origin is not replicating");
    }

    double Origin::getMaxCollisionTime(double vFork) const {
        if (state == State::ReplL) {
            return getLeftCollisionTime(vFork);
        }
        if (state == State::ReplR) {
            return getRightCollisionTime(vFork);
        }
        if (state == State::ReplLR) {
            return std::max(getLeftCollisionTime(vFork), getRightCollisionTime(vFork));
        }
        throw std::runtime_error("Origin is not replicating");
    }

    void Origin::fire(double tFire) {
        firingTime = tFire;
        state = State::ReplLR;
        leftReplOrigin = findLeftReplOrigin();
        if (leftReplOrigin != nullptr) {
            leftReplOrigin->rightReplOrigin = this;
        }
        rightReplOrigin = findRightReplOrigin();
        if (rightReplOrigin != nullptr) {
            rightReplOrigin->leftReplOrigin = this;
        }
    }

    int Origin::replicateLeft(double tCurrent, double vFork, Vector3<double> const **releasePos) {
        // replicate origin
        double const leftCollisionTime = getLeftCollisionTime(vFork);
        unsigned long int leftPos = data.pos - (unsigned long int) std::floor((tCurrent - firingTime) * vFork);
        if (leftCollisionTime <= tCurrent) {
            leftPos = data.pos - (unsigned long int) std::floor((leftCollisionTime - firingTime) * vFork);
            state = (state == State::ReplLR) ? State::ReplR : State::Post;
        }
        // passively fire neighboring origins
        int nPassivelyActivated = 0;
        while (nextLeftPassOrigin != nullptr && nextLeftPassOrigin->data.pos >= leftPos && nextLeftPassOrigin->state == State::Pre) {
            nextLeftPassOrigin->state = State::Pass;
            nextLeftPassOrigin->firingTime = firingTime + (data.pos - nextLeftPassOrigin->data.pos) / vFork;
            nextLeftPassOrigin = nextLeftPassOrigin->leftOrigin;
            nPassivelyActivated++;
        }
        return nPassivelyActivated;
    }

    int Origin::replicateRight(double tCurrent, double vFork, Vector3<double> const **releasePos) {
        // replicate origin
        double const rightCollisionTime = getRightCollisionTime(vFork);
        unsigned long int rightPos = data.pos + (unsigned long int) std::floor((tCurrent - firingTime) * vFork);
        if (rightCollisionTime <= tCurrent) {
            rightPos = data.pos + (unsigned long int) std::floor((rightCollisionTime - firingTime) * vFork);
            state = (state == State::ReplLR) ? State::ReplL : State::Post;
            if (rightReplOrigin != nullptr) {
                *releasePos = &chromosome->findGranule(rightPos).pos;
            }
        }
        // passively fire neighboring origins
        int nPassivelyActivated = 0;
        while (nextRightPassOrigin != nullptr && nextRightPassOrigin->data.pos <= rightPos && nextRightPassOrigin->state == State::Pre) {
            nextRightPassOrigin->state = State::Pass;
            nextRightPassOrigin->firingTime = firingTime + (nextRightPassOrigin->data.pos - data.pos) / vFork;
            nextRightPassOrigin = nextRightPassOrigin->rightOrigin;
            nPassivelyActivated++;
        }
        return nPassivelyActivated;
    }

    Origin *Origin::findLeftReplOrigin() const {
        Origin *o = leftOrigin;
        while (o != nullptr && o->state != State::ReplLR && o->state != State::ReplR) {
            o = o->leftOrigin;
        }
        return o;
    }

    Origin *Origin::findRightReplOrigin() const {
        Origin *o = rightOrigin;
        while (o != nullptr && o->state != State::ReplLR && o->state != State::ReplL) {
            o = o->rightOrigin;
        }
        return o;
    }

    double Origin::getLeftCollisionTime(double vFork) const {
        if (leftReplOrigin != nullptr) {
            return (firingTime + leftReplOrigin->firingTime + (data.pos - leftReplOrigin->data.pos) / vFork) / 2;
        }
        return firingTime + (data.pos - chromosome->findContig(data.pos).start) / vFork;
    }

    double Origin::getRightCollisionTime(double vFork) const {
        if (rightReplOrigin != nullptr) {
            return (firingTime + rightReplOrigin->firingTime + (rightReplOrigin->data.pos - data.pos) / vFork) / 2;
        }
        return firingTime + (chromosome->findContig(data.pos).end - data.pos) / vFork;
    }

}