#include "Particle.h"

namespace DNAReplication {

    Particle::Particle(bool active, Vector3<double> const &pos) : active(active), pos(pos) {
    }

    bool Particle::isActive() const {
        return active;
    }

    void Particle::setActive(bool active) {
        this->active = active;
    }

    Vector3<double> const &Particle::getPos() const {
        return pos;
    }

    void Particle::setPos(Vector3<double> const &pos) {
        this->pos = pos;
    }

    Origin *Particle::getBoundOrigin() const {
        return boundOrigin;
    }

    void Particle::setBoundOrigin(Origin *boundOrigin) {
        this->boundOrigin = boundOrigin;
    }

    void Particle::clearPreOriginsInProximity() {
        std::swap(preOriginsInProximity, previousPreOriginsInProximity);
        preOriginsInProximity.clear();
    }

    void Particle::addPreOriginInProximity(Origin *origin) {
        preOriginsInProximity.push_back(origin);
    }

    std::vector<Origin *> const &Particle::getPreOriginsInProximity() const {
        return preOriginsInProximity;
    }

    std::vector<Origin *> const &Particle::getPreviousPreOriginsInProximity() const {
        return previousPreOriginsInProximity;
    }

}