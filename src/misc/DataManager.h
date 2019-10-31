#ifndef DNA_REPLICATION_DATAMANAGER_H
#define DNA_REPLICATION_DATAMANAGER_H

#include <boost/filesystem.hpp>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../Chromosome.h"
#include "../Origin.h"
#include "Vector3.h"

#define CSV_DELIM ','

namespace DNAReplication {

    /**
     * Static class handling all data loading operations
     */
    class DataManager {

    public:

        /// Returns the canonical paths to all files in a directory (non-recursive)
        static std::vector<std::string> listFiles(std::string const &dirPath);

        /// Loads the origin information from the specified CSV file
        static std::vector<Origin::OriginData> loadOriginData(std::string const &filePath);

        /// Loads the chromosome information from the specified CSV file
        static std::vector<Chromosome::ChromosomeData> loadChromosomeData(std::string const &filePath);

        /// Initializes the granules (3D structure) of the specified chromosomes
        static void initializeChromosomeGranules(std::vector<Chromosome::ChromosomeData> &chromosomeData, std::string const &filePath);

    };

}

#endif //DNA_REPLICATION_DATAMANAGER_H
