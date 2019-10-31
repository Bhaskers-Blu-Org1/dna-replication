#include "MultiSimulationObserver.h"

namespace DNAReplication {

    void MultiSimulationObserver::handleOriginFired(SimulationOriginEvent const &e) {
        if (e.origin.getState() == Origin::State::ReplLR) {
            originFiringCounts[e.origin.getID()]++;
            originFiringTimeSums[e.origin.getID()] += e.simulation.getCurrentTime();
        }
    }

    std::vector<int> MultiSimulationObserver::getOriginFiringCounts(std::vector<Origin> const &origins) {
        std::vector<int> result;
        result.reserve(origins.size());
        for (Origin const &o : origins) {
            result.push_back(originFiringCounts[o.getID()]);
        }
        return result;
    }

    std::vector<double> MultiSimulationObserver::getOriginFiringTimeSums(const std::vector<Origin> &origins) {
        std::vector<double> result;
        result.reserve(origins.size());
        for (Origin const &o : origins) {
            result.push_back(originFiringTimeSums[o.getID()]);
        }
        return result;
    }

}