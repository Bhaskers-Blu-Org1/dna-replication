#ifndef DNA_REPLICATION_MULTISIMULATIONOBSERVER_H
#define DNA_REPLICATION_MULTISIMULATIONOBSERVER_H

#include <map>
#include <string>

#include "../SimulationObserver.h"

namespace DNAReplication {

    /**
     * Aggregates the origin firing events of multiple simulations
     */
    class MultiSimulationObserver : public SimulationObserver {

    private:

        std::map<std::string, int> originFiringCounts;
        std::map<std::string, double> originFiringTimeSums;

    public:

        void handleOriginFired(SimulationOriginEvent const &e) override;

        /// Returns the activation counts for the specified origins (in the same order)
        std::vector<int> getOriginFiringCounts(std::vector<Origin> const &origins);

        /// Returns the activation time sums for the specified origins (in the same order)
        std::vector<double> getOriginFiringTimeSums(std::vector<Origin> const &origins);

    };

}

#endif //DNA_REPLICATION_MULTISIMULATIONOBSERVER_H
