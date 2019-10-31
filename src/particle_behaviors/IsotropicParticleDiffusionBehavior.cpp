#include "IsotropicParticleDiffusionBehavior.h"

namespace DNAReplication {

    IsotropicParticleDiffusionBehavior::IsotropicParticleDiffusionBehavior(double h, double d, double r, double xNucl, double rSPB, double rPeriphery)
            : h(h), d(d), r(r), xNucl(xNucl), rSPB(rSPB), rPeriphery(rPeriphery), rng(boost::random::random_device{}()) {
    }

    double IsotropicParticleDiffusionBehavior::timeStep() {
        double const lambda = 6 * d / (h * h);
        boost::random::exponential_distribution<double> timeDist(lambda);
        return timeDist(rng);
    }

    Vector3<double> IsotropicParticleDiffusionBehavior::getRandomPosition() {
        boost::random::uniform_real_distribution<double> coordinateDist(-r, r);
        Vector3<double> pos{};
        do {
            pos = {coordinateDist(rng), coordinateDist(rng), coordinateDist(rng)};
        } while (!inDomain(pos));
        return pos;
    }

    bool IsotropicParticleDiffusionBehavior::inDomain(Vector3<double> const &pos) const {
        return inNucleus(pos) && !inNucleolus(pos);
    }

    bool IsotropicParticleDiffusionBehavior::inSPB(Vector3<double> const &pos) const {
        return (pos + Vector3<double> {r - rSPB, 0, 0}).getLengthSquared() <= rSPB * rSPB;
    }

    bool IsotropicParticleDiffusionBehavior::inPeriphery(Vector3<double> const &pos) const {
        if (inSPB(pos) || !inDomain(pos)) {
            return false;
        }
        return pos.getLengthSquared() >= rPeriphery * rPeriphery;
    }

    Vector3<double> IsotropicParticleDiffusionBehavior::diffuse(Vector3<double> const &pos) {
        boost::random::uniform_int_distribution<std::size_t> diffusionDist(0, diffusionMoves.size() - 1);
        return pos + diffusionMoves[diffusionDist(rng)] * h;
    }

    Vector3<double> IsotropicParticleDiffusionBehavior::reflect(Vector3<double> const &pos) {
        // find the (at least) three closest points within the domain
        std::size_t oldSize;
        std::vector<Vector3<double>> candidates{pos};
        do {
            oldSize = reflectExpand(candidates, 0);
            reflectSort(candidates, pos, oldSize);
        } while (candidates.size() <= 3);
        // construct matrix of the three closest candidates
        Eigen::Matrix3d candidateMatrix;
        candidateMatrix << candidates[1].x, candidates[2].x, candidates[3].x,
                candidates[1].y, candidates[2].y, candidates[3].y,
                candidates[1].z, candidates[2].z, candidates[3].z;
        // replace the third candidate until linear system can be solved
        std::size_t thirdCandidateIdx = 3;
        Eigen::ColPivHouseholderQR<Eigen::Matrix3d> qr(candidateMatrix);
        while (qr.rank() != 3) {
            if (++thirdCandidateIdx == candidates.size()) {
                oldSize = reflectExpand(candidates, oldSize);
                reflectSort(candidates, pos, oldSize);
            }
            candidateMatrix(0, 2) = candidates[thirdCandidateIdx].x;
            candidateMatrix(1, 2) = candidates[thirdCandidateIdx].y;
            candidateMatrix(2, 2) = candidates[thirdCandidateIdx].z;
            qr.compute(candidateMatrix);
        }
        // solve linear system (numeric instability --> take component-wise absolute values)
        Eigen::Vector3d const candidateCoef = qr.solve(Eigen::Vector3d(pos.x, pos.y, pos.z)).cwiseAbs();
        // sample reflected point from candidates according to their coefficients ("probabilities")
        std::array<double, 3> candidateCDF{{candidateCoef.x(), candidateCoef.y(), candidateCoef.z()}};
        double const candidateCoefSum = RandomTools::getCDF(candidateCDF.begin(), candidateCDF.end(), candidateCDF.begin());
        boost::random::uniform_real_distribution<double> candidateDist(0, candidateCoefSum);
        auto const candidateIdx = std::upper_bound(candidateCDF.begin(), candidateCDF.end(), candidateDist(rng)) - candidateCDF.begin();
        return candidates[candidateIdx == 2 ? thirdCandidateIdx : candidateIdx + 1];
    }

    bool IsotropicParticleDiffusionBehavior::inNucleus(Vector3<double> const &pos) const {
        return pos.getLengthSquared() <= r * r;
    }

    bool IsotropicParticleDiffusionBehavior::inNucleolus(Vector3<double> const &pos) const {
        return (pos - Vector3<double> {xNucl, 0, 0}).getLengthSquared() <= r * r;
    }

    std::size_t IsotropicParticleDiffusionBehavior::reflectExpand(std::vector<Vector3<double>> &vec, std::size_t startIdx) const {
        std::size_t const oldSize = vec.size();
        for (std::size_t i = startIdx; i < oldSize; ++i) {
            for (Vector3<double> const &reflectionMove : reflectionMoves) {
                Vector3<double> const candidate = vec[i] + reflectionMove * h;
                if (inDomain(candidate)) {
                    auto const last = vec.begin() + oldSize;
                    if (std::find(vec.begin(), last, candidate) == last) {
                        vec.push_back(candidate);
                    }
                }
            }
        }
        if (vec.size() == oldSize) {
            // If this happens, check if all granule positions are within the diffusion domain
            throw std::runtime_error("Could not expand reflection candidates");
        }
        return oldSize;
    }

    void IsotropicParticleDiffusionBehavior::reflectSort(std::vector<Vector3<double>> &vec, Vector3<double> const &refPos, std::size_t startIdx) const {
        if (inNucleolus(refPos)) {
            std::sort(vec.begin() + startIdx, vec.end(), [&refPos](Vector3<double> const &p1, Vector3<double> const &p2) {
                return (refPos - p1).getLengthSquared() < (refPos - p2).getLengthSquared();
            });
        } else {
            std::sort(vec.begin() + startIdx, vec.end(), [&refPos](Vector3<double> const &p1, Vector3<double> const &p2) {
                return p1.dotProduct(-refPos) < p2.dotProduct(-refPos);
            });
        }
    }

}