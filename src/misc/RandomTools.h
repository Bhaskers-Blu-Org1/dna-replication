#ifndef DNA_REPLICATION_RANDOMTOOLS_H
#define DNA_REPLICATION_RANDOMTOOLS_H

namespace DNAReplication {

    /**
     * Tools for random number generation
     */
    class RandomTools {

    public:

        /// Creates a cumulative "probability" vector (not normalized)
        /// \param probFirst Iterator pointing to the first "probability"
        /// \param probLast Iterator pointing to the element after the last "probability"
        /// \param cdfFirst Iterator pointing to the first cumulative "probability" (out parameter)
        /// \return The total sum of the "probabilities" (not normalized, equal to the last element)
        template<typename ProbIterator, typename CDFIterator>
        static double getCDF(ProbIterator probFirst, ProbIterator probLast, CDFIterator cdfFirst) {
            double cumsum = 0;
            while (probFirst != probLast) {
                cumsum += *probFirst;
                *cdfFirst = cumsum;
                probFirst++;
                cdfFirst++;
            }
            return cumsum;
        }

    };

}

#endif //DNA_REPLICATION_RANDOMTOOLS_H
