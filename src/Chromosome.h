#ifndef DNA_REPLICATION_CHROMOSOME_H
#define DNA_REPLICATION_CHROMOSOME_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "misc/Vector3.h"

#define CHROMOSOME_GRANULE_SIZE 3500

namespace DNAReplication {

    /**
     * Chromosome represented by a chain of granules
     *
     * Each granule (3D position) is representing a region of 3.5kbp on the genome sequence.
     * Additionally, the start- and end-positions of contigs are stored (in base-pairs).
     * DNA replication is simulated for sequenced regions ("contigs") only.
     */
    class Chromosome {

    public:

        /// The start- and end-position of a sequenced region ("contig") in base-pairs
        struct Contig {
            unsigned long int const start;
            unsigned long int const end;

            Contig(unsigned long int start, unsigned long int end) : start(start), end(end) {
            }
        };

        /// 3D position representing a 3.5kbp-stretch on the genome sequence
        struct Granule {
            Vector3<double> const pos;

            explicit Granule(Vector3<double> const &pos) : pos(pos) {
            }
        };

        struct ChromosomeData {

            friend class Chromosome;

        public:

            std::string id;
            std::vector<Contig> contigs{};
            std::vector<Granule> granules{};

            ChromosomeData(std::string const &id, std::vector<Contig> const &contigs) : id(id), contigs(contigs) {
                // granules are initialized later
            }

        private:

            ChromosomeData() = default;

        };

    private:

        ChromosomeData const data;

    public:

        /// Constructs a new chromosome
        explicit Chromosome(ChromosomeData data);

        /// Return the ID (i.e., the number) of the chromosome
        std::string const &getID() const;

        /// Two chromosomes are considered equal if their ids match
        bool operator==(Chromosome const &other) const;

        /// Returns the contig belonging to a given base-pair position
        Contig const &findContig(unsigned long int pos) const;

        /// Returns the granule belonging to a given base-pair position
        Granule const &findGranule(unsigned long int pos) const;

        /// Efficiently checks whether two base-pair positions belong to the same granule
        bool inSameContig(unsigned long int pos1, unsigned long int pos2) const;

        // Disallow copying

        Chromosome(Chromosome const &copy) {
            throw std::runtime_error("Internal error: copying chromosomes is not allowed");
        }

        Chromosome &operator=(Chromosome const &copy) {
            throw std::runtime_error("Internal error: copying chromosomes is not allowed");
        }

    };

}

#endif //DNA_REPLICATION_CHROMOSOME_H
