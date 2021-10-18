# EventPlaneROOT

Analysis code for calculating event plane resolution for FIT detectors (FV0 and FT0) from O2 simulation data.

## Usage:

 * ./saveQvecs <outfile.root> <\path\to\input> <cent>
    * calculates Q-vectors from simulation files in given path (check the path if files not found). outfile is the output file name, give centtrality class (for example 20-30)
 * ./calcCorrections <infile.root> <cent>
    * calculates the corrections for the even plane from raw q-vectors. Needs input file generated with saveQvecs without corrections (turn bDoCorrections off in SaveQvecs.h)
    * run this if no correction files exist or they need updating
 * ./calcEventPlaneRes <infile.root> <outfile.root>
    * calculates the event planes and the components that are needed to calculate the resolution. These are saved into histograms for each detector
