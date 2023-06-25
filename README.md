# Sim-UAV
A fixed-wing UAV flight Simulator able to compute flight paths based on input command with a PID controller.

## Executable
The simulator executable (SimUAV.exe) can be found in SimUAV/cmake-build-debug/SimUAV.exe
## Main [SRC]
This source file is the main file where all the sub-function have been called in the code.
Here you can choose the principle parameters of Velocity, Flight height and ramp angle.
The Trim conditions are screen printed, calculated by a previous interpolation of the UAV database, and now you can choose the Path desired.

### Database
Here all the required files for the simulation are uploaded

### Path
This directory contains the prebuilt paths for the UAV considering standard conditions. To run a different simulation, please add the .txt file to this directory before starting the simulation.

### Google_tests
It contains a library of test functions to verify the proper performance of the functions implemented

### Preprocessing
The preprocessing phase of reading and saving all data.




All the errors and warnings are listed at the link below:
[https://docs.google.com/spreadsheets/d/1dff7mPkpre8ac-bXZYAgpPNy8UBMbChXXmcD9r2qzHQ/edit#gid=0]
