#include <fstream>
#include <mpi.h>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "Simulation.h"
#include "misc/DataManager.h"
#include "particle_behaviors/IsotropicParticleDiffusionBehavior.h"
#include "particle_behaviors/ProbabilisticParticleActivationBehavior.h"
#include "particle_behaviors/ProbabilisticParticleBindingBehavior.h"
#include "simulation_observers/MultiSimulationObserver.h"

#define MY_SUCCESS 0
#define MY_ERROR 1

#define MYOPT_HELP "help"
#define MYOPT_NITER "niter"
#define MYOPT_OUTDIR "outdir"
#define MYOPT_OUTKEY "outkey"

#define MYOPT_ORIFILE "orifile"
#define MYOPT_ASSYFILE "assyfile"
#define MYOPT_STRUCTDIR "structdir"

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

#define MYMPI_ROOT 0

#define MYMPI_MSG_WORKER_READY "WORKER_READY"
#define MYMPI_MSG_WORKER_RUN "WORKER_RUN"
#define MYMPI_MSG_WORKER_RESULT "WORKER_RESULT"
#define MYMPI_MSG_WORKER_EXIT "WORKER_EXIT"

using namespace DNAReplication;

void MPI_SendString(std::string value, int dest) {
    MPI_Send(&value[0], (int) value.length(), MPI_CHAR, dest, 0, MPI_COMM_WORLD);
}

int MPI_RecvString(std::string &value, int source) {
    int count;
    MPI_Status status{};
    MPI_Probe(source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &count);
    value.resize((std::size_t) count);
    MPI_Recv(&value[0], count, MPI_CHAR, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return status.MPI_SOURCE;
}

void runMaster(boost::program_options::variables_map &vm, std::vector<Origin::OriginData> const &originData, std::vector<Chromosome::ChromosomeData> &chromosomeData) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    auto const &structDir = vm[MYOPT_STRUCTDIR].as<std::string>();
    std::vector<std::string> structFiles = DataManager::listFiles(structDir);
    // prepare files
    auto const &outDir = vm[MYOPT_OUTDIR].as<std::string>();
    auto const &outKey = vm[MYOPT_OUTKEY].as<std::string>();
    std::ofstream completionTimesFileStream((boost::filesystem::path(outDir) / boost::filesystem::path(outKey + "_completionTimes.csv")).string());
    std::ofstream firingTimesFileStream((boost::filesystem::path(outDir) / boost::filesystem::path(outKey + "_firingTimes.csv")).string());
    completionTimesFileStream << "STRUCTURE,ITERATION,COMPLETION_TIME" << std::endl;
    firingTimesFileStream << "STRUCTURE,ITERATION";
    for (Origin::OriginData const &data : originData) {
        firingTimesFileStream << "," + data.id;
    }
    firingTimesFileStream << std::endl;
    // master loop
    int iterCount = 0;
    int numActiveWorkers = size - 1;
    auto structFileIter = structFiles.begin();
    while (numActiveWorkers > 0) {
        // receive message
        std::string workerMessage;
        int const workerRank = MPI_RecvString(workerMessage, MPI_ANY_SOURCE);
        if (workerMessage == std::string(MYMPI_MSG_WORKER_READY)) {
            // start or end worker
            if (iterCount >= vm[MYOPT_NITER].as<int>()) {
                structFileIter++;
                iterCount = 0;
            }
            if (structFileIter != structFiles.end()) {
                MPI_SendString(MYMPI_MSG_WORKER_RUN, workerRank);
                MPI_SendString(*structFileIter, workerRank);
                MPI_Send(&iterCount, 1, MPI_INT, workerRank, 0, MPI_COMM_WORLD);
                iterCount++;
            } else {
                MPI_SendString(MYMPI_MSG_WORKER_EXIT, workerRank);
                numActiveWorkers--;
            }
        } else if (workerMessage == std::string(MYMPI_MSG_WORKER_RESULT)) {
            // receive info
            int workerIteration;
            std::string workerGranuleFile;
            MPI_RecvString(workerGranuleFile, workerRank);
            workerGranuleFile = boost::filesystem::path(workerGranuleFile).filename().string();
            MPI_Recv(&workerIteration, 1, MPI_INT, workerRank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // receive results
            double workerSimulationTime;
            std::vector<int> workerOriginFiringCounts(originData.size(), 0);
            std::vector<double> workerOriginFiringTimes(originData.size(), 0);
            MPI_Recv(&workerSimulationTime, 1, MPI_DOUBLE, workerRank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&workerOriginFiringCounts[0], (int) workerOriginFiringCounts.size(), MPI_INT, workerRank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&workerOriginFiringTimes[0], (int) workerOriginFiringTimes.size(), MPI_DOUBLE, workerRank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // write results
            completionTimesFileStream << workerGranuleFile << ',' << std::to_string(workerIteration) << ',' << std::to_string(workerSimulationTime) << std::endl;
            firingTimesFileStream << workerGranuleFile << ',' << std::to_string(workerIteration);
            for (std::size_t i = 0; i < originData.size(); ++i) {
                firingTimesFileStream << ',' << workerOriginFiringTimes[i];
            }
            firingTimesFileStream << std::endl;
        } else {
            throw std::runtime_error("Unknown MPI message received in master process");
        }
    }
}

void runWorker(boost::program_options::variables_map const &vm, std::vector<Origin::OriginData> &originData, std::vector<Chromosome::ChromosomeData> &chromosomeData) {
    // worker loop
    bool workerActive = true;
    while (workerActive) {
        MPI_SendString(MYMPI_MSG_WORKER_READY, MYMPI_ROOT);
        // receive message
        std::string masterMessage;
        MPI_RecvString(masterMessage, MYMPI_ROOT);
        if (masterMessage == std::string(MYMPI_MSG_WORKER_RUN)) {
            // receive info
            int iteration;
            std::string structFile;
            MPI_RecvString(structFile, MYMPI_ROOT);
            MPI_Recv(&iteration, 1, MPI_INT, MYMPI_ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // determine parameters
            DataManager::initializeChromosomeGranules(chromosomeData, structFile);
            bool const spbActivationEnabled = vm.count(MYOPT_RSPB) > 0;
            bool const peripheryInactivationEnabled = vm.count(MYOPT_RPERY) > 0;
            double const rSPB = spbActivationEnabled ? vm[MYOPT_RSPB].as<double>() : 0;
            double const rPeriphery = peripheryInactivationEnabled ? vm[MYOPT_RPERY].as<double>() : 0;
            IsotropicParticleDiffusionBehavior diffusionBehavior(vm[MYOPT_HGRID].as<double>(), vm[MYOPT_DCOEF].as<double>(), vm[MYOPT_RNUCL].as<double>(),
                                                                 vm[MYOPT_XNUCL].as<double>(), rSPB, rPeriphery);
            ProbabilisticParticleActivationBehavior activationBehavior(vm[MYOPT_PACT].as<double>(), spbActivationEnabled, peripheryInactivationEnabled);
            ProbabilisticParticleBindingBehavior bindingBehavior(vm[MYOPT_DBIND].as<double>(), vm[MYOPT_PBIND].as<double>());
            // run simulation
            Simulation sim(vm[MYOPT_VFORK].as<double>(), originData, chromosomeData, diffusionBehavior, activationBehavior, bindingBehavior);
            sim.initializeParticles(vm[MYOPT_NPART].as<unsigned short int>());
            std::shared_ptr<MultiSimulationObserver> originObserver = std::make_shared<MultiSimulationObserver>();
            sim.registerObserver(originObserver);
            sim.run();
            sim.clearObservers();
            // collect results
            double t = sim.getCurrentTime();
            std::vector<int> originFiringCounts = originObserver->getOriginFiringCounts(sim.getOrigins());
            std::vector<double> originFiringTimes = originObserver->getOriginFiringTimeSums(sim.getOrigins());
            // send info
            MPI_SendString(MYMPI_MSG_WORKER_RESULT, MYMPI_ROOT);
            MPI_SendString(structFile, MYMPI_ROOT);
            MPI_Send(&iteration, 1, MPI_INT, MYMPI_ROOT, 0, MPI_COMM_WORLD);
            // send results
            MPI_Send(&t, 1, MPI_DOUBLE, MYMPI_ROOT, 0, MPI_COMM_WORLD);
            MPI_Send(&originFiringCounts[0], (int) originFiringCounts.size(), MPI_INT, MYMPI_ROOT, 0, MPI_COMM_WORLD);
            MPI_Send(&originFiringTimes[0], (int) originFiringTimes.size(), MPI_DOUBLE, MYMPI_ROOT, 0, MPI_COMM_WORLD);
        } else if (masterMessage == std::string(MYMPI_MSG_WORKER_EXIT)) {
            workerActive = false;
        } else {
            throw std::runtime_error("Unknown MPI message received in worker process");
        }
    }
}

int main(int argc, char *argv[]) {
    // initialize MPI
    int rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // client options
    boost::program_options::options_description opts_client("Client options");
    opts_client.add_options()
            (MYOPT_HELP ",h", "Display help message")
            (MYOPT_NITER ",i", boost::program_options::value<int>()->required(), "Number of iterations (required)")
            (MYOPT_OUTDIR ",w", boost::program_options::value<std::string>()->required(), "Path to output directory (required)")
            (MYOPT_OUTKEY ",k", boost::program_options::value<std::string>()->required(), "Simulation key (for output files, required)");
    // input data
    boost::program_options::options_description opts_indata("Input data");
    opts_indata.add_options()
            (MYOPT_ORIFILE ",o", boost::program_options::value<std::string>()->required(), "Path to origin positions file (required)")
            (MYOPT_ASSYFILE ",c", boost::program_options::value<std::string>()->required(), "Path to assembly information file (required)")
            (MYOPT_STRUCTDIR ",s", boost::program_options::value<std::string>()->required(), "Path genome structure directory (required)");
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
            std::cout << "Runs multiple DNA replication simulations using MPI." << std::endl;
            std::cout << opts << std::endl;
            MPI_Finalize();
            return MY_SUCCESS;
        }
        boost::program_options::notify(vm);
    } catch (boost::program_options::error &e) {
        if (rank == MYMPI_ROOT) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            std::cerr << "Specify option '--" MYOPT_HELP "' to display a help message." << std::endl;
        }
        MPI_Finalize();
        return MY_ERROR;
    }
    // run simulations
    std::vector<Origin::OriginData> originData = DataManager::loadOriginData(vm[MYOPT_ORIFILE].as<std::string>());
    std::vector<Chromosome::ChromosomeData> chromosomeData = DataManager::loadChromosomeData(vm[MYOPT_ASSYFILE].as<std::string>());
    if (rank == MYMPI_ROOT) {
        runMaster(vm, originData, chromosomeData);
    } else {
        runWorker(vm, originData, chromosomeData);
    }
    // finalize MPI
    MPI_Finalize();
    return MY_SUCCESS;
}
