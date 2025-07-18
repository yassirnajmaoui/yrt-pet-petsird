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
2. Put the `libyrtpet.a` library in a directory visible to your linker and compiler
and YRT-PET's `include` folder in a directory visible to your compiler.

## Usage

