#!/bin/bash
ELECTRON_SIM_PATH="${HOME}/alice/O2Simulation/e_gun_2GeV/generated_files/trddigits.root"
PION_SIM_PATH="${HOME}/alice/O2Simulation/p_gun_2GeV/generated_files/trddigits.root"

for i in {1..50}; do
    make
    ELEC_CSV_PATH=`printf "${HOME}/alice/O2Simulation/electron_csvs/%04d" \"${i}\"`
    PION_CSV_PATH=`printf "${HOME}/alice/O2Simulation/pion_csvs/%04d" \"${i}\"`
    root -b -q -l "tbsumDigits.C(\"${ELECTRON_SIM_PATH}\" , \"${PION_SIM_PATH}\", \"${ELEC_CSV_PATH}\" ,\"${PION_CSV_PATH}\")"
    make clean
done
