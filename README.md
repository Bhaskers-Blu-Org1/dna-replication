# dna-replication

A stochastic hybrid model for DNA replication incorporating protein mobility dynamics

If you find this code useful in your research, please consider citing:
```
@article{Windhager2019,
  doi = {10.1101/583187},
  url = {https://doi.org/10.1101/583187},
  year = {2019},
  month = mar,
  publisher = {Cold Spring Harbor Laboratory},
  author = {Jonas Windhager and Amelia Paine and Patroula Nathanailidou and Eve Tasiudi and Maria Rodriguez Martinez and Zoi Lygerou and John Lygeros and Maria Anna Rapsomaniki},
  title = {A stochastic hybrid model for {DNA} replication incorporating protein mobility dynamics}
}
```

## Prerequisites

* C++11 compiler
* CMake 3.6 or newer
* [Eigen 3.3.4](https://eigen.tuxfamily.org)

To clone the Eigen library directly from GitHub using git:
```bash
git clone -b 3.3.4 https://github.com/eigenteam/eigen-git-mirror.git /path/to/eigen3
```

## Installation

Replace `/path/to/eigen3` with the path to your local copy of the Eigen library and execute the following commands:
```bash
git clone https://github.com/jwindhager/dna-replication.git
mkdir dna-replication/build
cd dna-replication/build
cmake -DCMAKE_BUILD_TYPE=Release -DEIGEN3_INCLUDE_DIR=/path/to/eigen3 ..
make
make install
```

## Usage

To perform a single DNA replication simulation, call `dnarepl`:
```
Runs a single DNA replication simulation.

Client options:
  -h [ --help ]            Display help message

Input data:
  -o [ --orifile ] arg     Path to origin positions file (required)
  -c [ --assyfile ] arg    Path to assembly information file (required)
  -s [ --structfile ] arg  Path genome structure file (required)

Simulation parameters:
  -r [ --rnucl ] arg       Nucleus radius (required, in um)
  -x [ --xnucl ] arg       Nucleolus displacement (required, in um)
  -q [ --rpery ] arg       Periphery radius (set for peripheral particle inactivation, in um)
  -z [ --rspb ] arg        Spindle pole body radius (enables SPB-mediated particle activation, in um)
  -n [ --npart ] arg       Number of activation factors (required)
  -g [ --hgrid ] arg       Step size of the diffusion grid (required, in um)
  -a [ --pact ] arg        Activation probability (for SPB-mediated particle activation)
  -d [ --dcoef ] arg       Effective diffusion coefficient (required, in um2/s)
  -b [ --dbind ] arg       Maximal binding distance (required, in um)
  -p [ --pbind ] arg       Binding probability (required)
  -f [ --vfork ] arg       Replication fork velocity (required, in b/s)
```

To perform a batch of DNA replication simulations, call `dnarepl-batch` on any MPI-enabled environment:
```
Runs multiple DNA replication simulations using MPI.

Client options:
  -h [ --help ]           Display help message
  -i [ --niter ] arg      Number of iterations (required)
  -w [ --outdir ] arg     Path to output directory (required)
  -k [ --outkey ] arg     Simulation key (for output files, required)

Input data:
  -o [ --orifile ] arg    Path to origin positions file (required)
  -c [ --assyfile ] arg   Path to assembly information file (required)
  -s [ --structdir ] arg  Path genome structure directory (required)

Simulation parameters:
  -r [ --rnucl ] arg      Nucleus radius (required, in um)
  -x [ --xnucl ] arg      Nucleolus displacement (required, in um)
  -q [ --rpery ] arg      Periphery radius (set for peripheral particle inactivation, in um)
  -z [ --rspb ] arg       Spindle pole body radius (enables SPB-mediated particle activation, in um)
  -n [ --npart ] arg      Number of activation factors (required)
  -g [ --hgrid ] arg      Step size of the diffusion grid (required, in um)
  -a [ --pact ] arg       Activation probability (for SPB-mediated particle activation)
  -d [ --dcoef ] arg      Effective diffusion coefficient (required, in um2/s)
  -b [ --dbind ] arg      Maximal binding distance (required, in um)
  -p [ --pbind ] arg      Binding probability (required)
  -f [ --vfork ] arg      Replication fork velocity (required, in b/s)
```

## License

This project is licensed under the terms of the MIT license.
