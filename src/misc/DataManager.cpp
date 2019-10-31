#include "DataManager.h"

namespace DNAReplication {

    std::vector<std::string> DataManager::listFiles(std::string const &dirPath) {
        std::vector<std::string> files;
        boost::filesystem::directory_iterator const endIter;
        for (boost::filesystem::directory_iterator iter(dirPath); iter != endIter; ++iter) {
            boost::filesystem::path path = boost::filesystem::canonical(iter->path());
            if (boost::filesystem::is_regular_file(path)) {
                files.push_back(path.string());
            }
        }
        return files;
    }

    std::vector<Origin::OriginData> DataManager::loadOriginData(std::string const &filePath) {
        std::ifstream fileStream(filePath);
        if (!fileStream.is_open()) {
            throw std::runtime_error("Could not open file " + filePath);
        }
        std::string fileLine;
        std::vector<Origin::OriginData> originData;
        while (std::getline(fileStream, fileLine)) {
            std::string originID, originChromosomeID, originPos;
            std::istringstream lineStream(fileLine);
            std::getline(lineStream, originID, CSV_DELIM);
            std::getline(lineStream, originChromosomeID, CSV_DELIM);
            std::getline(lineStream, originPos, CSV_DELIM);
            originData.emplace_back(originID, originChromosomeID, std::stoul(originPos));
        }
        return originData;
    }

    std::vector<Chromosome::ChromosomeData> DataManager::loadChromosomeData(std::string const &filePath) {
        std::ifstream fileStream(filePath);
        if (!fileStream.is_open()) {
            throw std::runtime_error("Could not open file " + filePath);
        }
        std::string fileLine;
        std::vector<Chromosome::ChromosomeData> chromosomeData;
        while (std::getline(fileStream, fileLine)) {
            std::string chromosomeID;
            std::istringstream lineStream(fileLine);
            std::getline(lineStream, chromosomeID, CSV_DELIM);
            std::vector<Chromosome::Contig> chromosomeContigs;
            while (!lineStream.eof()) {
                std::string chromosomeContigStart, chromosomeContigEnd;
                std::getline(lineStream, chromosomeContigStart, CSV_DELIM);
                std::getline(lineStream, chromosomeContigEnd, CSV_DELIM);
                chromosomeContigs.emplace_back(std::stoul(chromosomeContigStart), std::stoul(chromosomeContigEnd));
            }
            chromosomeData.emplace_back(chromosomeID, chromosomeContigs);
        }
        return chromosomeData;
    }

    void DataManager::initializeChromosomeGranules(std::vector<Chromosome::ChromosomeData> &chromosomeData, std::string const &filePath) {
        std::ifstream fileStream(filePath);
        if (!fileStream.is_open()) {
            throw std::runtime_error("Could not open file " + filePath);
        }
        std::string fileLine;
        std::map<std::string, std::size_t> chromosomeDataIndices;
        for (std::size_t i = 0; i < chromosomeData.size(); ++i) {
            chromosomeData[i].granules.clear();
            chromosomeDataIndices[chromosomeData[i].id] = i;
        }
        while (std::getline(fileStream, fileLine)) {
            std::string granuleChromosomeID, granulePosX, granulePosY, granulePosZ;
            std::istringstream lineStream(fileLine);
            std::getline(lineStream, granuleChromosomeID, CSV_DELIM);
            std::getline(lineStream, granulePosX, CSV_DELIM);
            std::getline(lineStream, granulePosY, CSV_DELIM);
            std::getline(lineStream, granulePosZ, CSV_DELIM);
            Vector3<double> const granulePos{std::stod(granulePosX), std::stod(granulePosY), std::stod(granulePosZ)};
            chromosomeData[chromosomeDataIndices.at(granuleChromosomeID)].granules.emplace_back(granulePos);
        }
    }

}