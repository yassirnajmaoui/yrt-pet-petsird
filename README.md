# PETSIRD reconstruction using YRT-PET

The purpose of this repo is to provide a reconstruction pipeline that uses
PETSIRD as data input and YRT-PET as a reconstruction engine.

## Background

The [Emission Tomography Standardization Initiative (ETSI)](https://etsinitiative.org/)
is working towards establishing a standard for PET Raw Data, called PETSIRD
("PET ETSI Raw Data").
More information is on https://github.com/ETSInitiative/PETSIRD.

## Build

1. Compile or download the binaries for
   [YRT-PET](https://github.com/YaleBioImaging/yrt-pet).
2. Make sure your environment has the dependencies that
   [PETSIRD](https://github.com/ETSInitiative/PETSIRD) requires
3. Put the `libyrtpet.a` library in a directory visible to your linker and
   compiler and YRT-PET's `include` folder in a directory visible to your
   compiler.
    - Note: The compiler reads the `CPLUS_INCLUDE_PATH` (for the includes) and
       `LIBRARY_PATH` (for the libraries). The linker reads `LD_LIBRARY_PATH` to
       link with other libraries.
    - Note: This step is automated if you run ``cmake --install .`` from the
       from your [YRT-PET](https://github.com/YaleBioImaging/yrt-pet) build
       directory
4. Run:
```
git clone <this repository's URL>
cd yrt-pet-petsird
git clone <the PETSIRD repository URL>
cd PETSIRD/model
yardl generate
cd ../..
mkdir build
cd build
cmake ../src -DUSE_CUDA=[ON/OFF]
make
```
5. The `USE_CUDA` option allows GPU reconstruction. Note that if you compiled
   YRT-PET with `USE_CUDA=OFF`, you might get compilation errors.

## Usage

The project compiles into a program names `petsird_yrtpet_reconstruct`, which
reconstructs from a given PETSIRD list-mode. Run `petsird_yrtpet_reconstruct -h`
for more information.
