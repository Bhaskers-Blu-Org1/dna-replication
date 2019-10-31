#include "Chromosome.h"

namespace DNAReplication {

    Chromosome::Chromosome(ChromosomeData data) : data(std::move(data)) {
    }

    std::string const &Chromosome::getID() const {
        return data.id;
    }

    bool Chromosome::operator==(Chromosome const &other) const {
        return data.id == other.data.id;
    }

    Chromosome::Contig const &Chromosome::findContig(unsigned long int pos) const {
        for (Contig const &contig : data.contigs) {
            if (contig.start <= pos && contig.end >= pos) {
                return contig;
            }
        }
        throw std::runtime_error("Position does not belong to any contig");
    }

    Chromosome::Granule const &Chromosome::findGranule(unsigned long int pos) const {
        auto const granuleIdx = (std::size_t) std::floor(pos / CHROMOSOME_GRANULE_SIZE);
        if (granuleIdx >= data.granules.size()) {
            throw std::runtime_error("Position does not belong to any granule");
        }
        return data.granules[granuleIdx];
    }

    bool Chromosome::inSameContig(unsigned long int pos1, unsigned long int pos2) const {
        for (Contig const &contig : data.contigs) {
            if (contig.start <= pos1 && contig.end >= pos1 && contig.start <= pos2 && contig.end >= pos2) {
                return true;
            }
        }
        return false;
    }

}