# EventPlaneROOT

Analysis code for calculating event plane resolution for FIT detectors (FV0 and FT0) from O2 simulation data.

## Usage:

 * ./saveQvecs <output.root> <\path\to\input> <comment> <docorrections>
    * calculates Q-vectors from simulation files in given path (check the path if files not found). outputi.root the output file name, comment gives label to the correction files (for example cent class is good)
 * ./calcCorrections <output.root> <\path\to\input> <comment>
    * same arguments as above, runs first saveQvecs and then calculates corrections from that output
    * run this if no correction files exist or they need updating
 * ./calcEventPlaneRes <output.root> <\path\to\input> <correction file label> <docorrections>
    * calculates the event planes and the components that are needed to calculate the resolution. These are saved into histograms for each detector
    * runs first saveQvecs and saves Q-vectors to output.root, output for resolutions is labeled as "res_output.root"


## TO-DO

 * Nothing at this point
