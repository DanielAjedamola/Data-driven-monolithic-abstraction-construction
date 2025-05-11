### This repository contains code for the paper titled: 'Data-driven controller synthesis via finite abstractions with formal guarantees'
### Link: https://ieeexplore.ieee.org/abstract/document/10314434
### Conference: LCSS-ACC 2024
### Authors: Daniel Ajeleye, Abolfazl Lavaei, and Majid Zamani
### Contact: daniel.ajeleye@colorado.edu
### Affiliation: University of Colorado Boulder

##### Before doing anything with the provided code, please do the following:
##### 1. visit here: https://github.com/mkhaled87/scots-ready for the details on installation and usage of SCOTS, and ensure the vehicle example is working on your device.
##### 2. Ensure gurobi is installed on your device, and place the header files gurobi_c.h, gurobi_c++.h, and the license file gurobi.lic in the same folder as your vehicle.cc file. 
##### 3. ensure that the directory in the make file are properly adjusted as applicable on your local device.
##### 4. attach the file init_abstr.cc to the same folder containing vehicle.cc, and first run the file with $make via the makefile for the part of obtaining an optimized grid parameter.
##### 5. replace the files vehicle.cc, vehicle.m with the one provided here and include vehicleModel.cc file as well, and proceed with the same instrcution to compile the files as it was done on the referenced repository in step 1.
