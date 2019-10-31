#include <fstream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "Simulation.h"
#include "misc/DataManager.h"
#include "particle_behaviors/IsotropicParticleDiffusionBehavior.h"
#include "particle_behaviors/ProbabilisticParticleActivationBehavior.h"
#include "particle_behaviors/ProbabilisticParticleBindingBehavior.h"

#define MY_SUCCESS 0
#define MY_ERROR 1

#define MYOPT_HELP "help"

#define MYOPT_ORIFILE "orifile"
#define MYOPT_ASSYFILE "assyfile"
#define MYOPT_STRUCTFILE "structfile"

#define MYOPT_RNUCL "rnucl"
#define MYOPT_XNUCL "xnucl"
#define MYOPT_RPERY "rpery"
#define MYOPT_RSPB "rspb"
#define MYOPT_NPART "npart"
#define MYOPT_HGRID "hgrid"
#define MYOPT_PACT "pact"
#define MYOPT_DCOEF "dcoef"
#define MYOPT_DBIND "dbind"
#define MYOPT_PBIND "pbind"
#define MYOPT_VFORK "vfork"

using namespace DNAReplication;

int main(int argc, char *argv[]) {
    // client options
    boost::program_options::options_description opts_client("Client options");
    opts_client.add_options()
            (MYOPT_HELP ",h", "Display help message");
    // input data
    boost::program_options::options_description opts_indata("Input data");
    opts_indata.add_options()
            (MYOPT_ORIFILE ",o", boost::program_options::value<std::string>()->required(), "Path to origin positions file (required)")
            (MYOPT_ASSYFILE ",c", boost::program_options::value<std::string>()->required(), "Path to assembly information file (required)")
            (MYOPT_STRUCTFILE ",s", boost::program_options::value<std::string>()->required(), "Path genome structure file (required)");
    // simulation parameters
    boost::program_options::options_description opts_simulation("Simulation parameters");
    opts_simulation.add_options()
            (MYOPT_RNUCL ",r", boost::program_options::value<double>()->required(), "Nucleus radius (required, in um)")
            (MYOPT_XNUCL ",x", boost::program_options::value<double>()->required(), "Nucleolus displacement (required, in um)")
            (MYOPT_RPERY ",q", boost::program_options::value<double>(), "Periphery radius (set for peripheral particle inactivation, in um)")
            (MYOPT_RSPB ",z", boost::program_options::value<double>(), "Spindle pole body radius (enables SPB-mediated particle activation, in um)")
            (MYOPT_NPART ",n", boost::program_options::value<unsigned short int>()->required(), "Number of activation factors (required)")
            (MYOPT_HGRID ",g", boost::program_options::value<double>()->required(), "Step size of the diffusion grid (required, in um)")
            (MYOPT_PACT ",a", boost::program_options::value<double>()->required(), "Activation probability (for SPB-mediated particle activation)")
            (MYOPT_DCOEF ",d", boost::program_options::value<double>()->required(), "Effective diffusion coefficient (required, in um2/s)")
            (MYOPT_DBIND ",b", boost::program_options::value<double>()->required(), "Maximal binding distance (required, in um)")
            (MYOPT_PBIND ",p", boost::program_options::value<double>()->required(), "Binding probability (required)")
            (MYOPT_VFORK ",f", boost::program_options::value<double>()->required(), "Replication fork velocity (required, in b/s)");
    // parse command line
    boost::program_options::variables_map vm;
    try {
        boost::program_options::options_description opts;
        opts.add(opts_client).add(opts_indata).add(opts_simulation);
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, opts), vm);
        if (vm.count(MYOPT_HELP) > 0) {
            std::cout << "Runs a single DNA replication simulation." << std::endl;
            std::cout << opts << std::endl;
            return MY_SUCCESS;
        }
        boost::program_options::notify(vm);
    } catch (boost::program_options::error &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "Specify option '--" MYOPT_HELP "' to display a help message." << std::endl;
        return MY_ERROR;
    }           
    // determine parameters
    bool const spbActivationEnabled = vm.count(MYOPT_RSPB) > 0;
    bool const peripheryInactivationEnabled = vm.count(MYOPT_RPERY) > 0;
    double const rSPB = spbActivationEnabled ? vm[MYOPT_RSPB].as<double>() : 0;
    double const rPeriphery = peripheryInactivationEnabled ? vm[MYOPT_RPERY].as<double>() : 0;
    IsotropicParticleDiffusionBehavior diffusionBehavior(vm[MYOPT_HGRID].as<double>(), vm[MYOPT_DCOEF].as<double>(), vm[MYOPT_RNUCL].as<double>(),
                                                         vm[MYOPT_XNUCL].as<double>(), rSPB, rPeriphery);
    ProbabilisticParticleActivationBehavior activationBehavior(vm[MYOPT_PACT].as<double>(), spbActivationEnabled, peripheryInactivationEnabled);
    ProbabilisticParticleBindingBehavior bindingBehavior(vm[MYOPT_DBIND].as<double>(), vm[MYOPT_PBIND].as<double>());
    // load data
    std::vector<Origin::OriginData> originData = DataManager::loadOriginData(vm[MYOPT_ORIFILE].as<std::string>());
    std::vector<Chromosome::ChromosomeData> chromosomeData = DataManager::loadChromosomeData(vm[MYOPT_ASSYFILE].as<std::string>());
    DataManager::initializeChromosomeGranules(chromosomeData, vm[MYOPT_STRUCTFILE].as<std::string>());
    // run simulation
    Simulation sim(vm[MYOPT_VFORK].as<double>(), originData, chromosomeData, diffusionBehavior, activationBehavior, bindingBehavior);
    sim.initializeParticles(vm[MYOPT_NPART].as<unsigned short int>());
    sim.run();
    return MY_SUCCESS;
}
