#ifndef DNA_REPLICATION_ISOTROPICPARTICLEDIFFUSIONBEHAVIOR_H
#define DNA_REPLICATION_ISOTROPICPARTICLEDIFFUSIONBEHAVIOR_H

#include <array>
#include <stdexcept>
#include <vector>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <Eigen/Dense>

#include "../ParticleDiffusionBehavior.h"
#include "../misc/RandomTools.h"
#include "../misc/Vector3.h"

namespace DNAReplication {

    /**
     * Isotropic diffusion in a spherical domain with Kushner-type boundary reflection
     *
     * Each particle diffuses on its own grid with a constant step size.
     * Gridding of time is based on the constant step size (see Kushner).
     * In one time step, particles can only move in one of six directions.
     *
     * delta_t ~ Exp(lambda), with lambda = 6 * d / h ^ 2
     *
     * The domain excludes the intersection with a spherical volume representing the nucleolus
     * (same radius as the spherical nucleus volume, spatially translated along the x-axis).
     *
     * At the boundary, particles are reflected using the following approach:
     * 1. Let x_new be the new position outside the domain boundary
     * 2. Define a plane at the point x_new with a normal n = -x_new (i.e., the radius)
     * 3. "Move" the plane towards the center until it has "hit" three candidate points forming a matrix C
     * 4. Solve the linear system of equations C * p = -x_new to obtain a candidate "probability" vector p
     * 5. As long as the system doesn't have a unique solution, keep replacing candidate points as before
     * 6. Sample the reflected position x_reflected from C according to the "probabilities" in vector p
     */
    class IsotropicParticleDiffusionBehavior : public ParticleDiffusionBehavior {

    private:

        // the six allowed diffusion directions (moves on the grid)
        const std::array<Vector3<double>, 7> diffusionMoves = {{{0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}};

        // moves to all 26 direct neighbors (diffusion directions + diagonal moves)
        const std::array<Vector3<double>, 26> reflectionMoves = {{
                                                                         {0, 0, 1}, {0, 0, -1}, {0, 1, 0}, {0, 1, 1}, {0, 1, -1}, {0, -1, 0}, {0, -1, 1}, {0, -1, -1},
                                                                         {1, 0, 0}, {1, 0, 1}, {1, 0, -1}, {1, 1, 0}, {1, 1, 1}, {1, 1, -1}, {1, -1, 0}, {1, -1, 1}, {1, -1, -1},
                                                                         {-1, 0, 0}, {-1, 0, 1}, {-1, 0, -1}, {-1, 1, 0}, {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, -1, -1}
                                                                 }};

        const double h;
        const double d;
        const double r;
        const double xNucl;
        const double rSPB;
        const double rPeriphery;
        boost::random::mt19937 rng;

        bool inNucleus(Vector3<double> const &pos) const;

        bool inNucleolus(Vector3<double> const &pos) const;

        std::size_t reflectExpand(std::vector<Vector3<double>> &vec, std::size_t startIdx) const;

        void reflectSort(std::vector<Vector3<double>> &vec, const Vector3<double> &refPos, std::size_t startIdx) const;

    public:

        /// Constructs a new instance
        /// \param h The diffusion step size (in um)
        /// \param d The diffusion coefficient (in um2s-1)
        /// \param r The radius of the spherical nucleus domain (in um)
        /// \param xNucl The x-coordinate of the nucleolus center (in um, relative to nucleus center)
        /// \param rSPB The radius of the spindle pole body (SPB) region (in um)
        /// \param rPeriphery The radius of the peripheral region (in um)
        IsotropicParticleDiffusionBehavior(double h, double d, double r, double xNucl, double rSPB, double rPeriphery);

        double timeStep() override;

        Vector3<double> getRandomPosition() override;

        bool inDomain(Vector3<double> const &pos) const override;

        bool inSPB(Vector3<double> const &pos) const override;

        bool inPeriphery(Vector3<double> const &pos) const override;

        Vector3<double> diffuse(const Vector3<double> &pos) override;

        Vector3<double> reflect(const Vector3<double> &pos) override;

    };

}

#endif //DNA_REPLICATION_ISOTROPICPARTICLEDIFFUSIONBEHAVIOR_H