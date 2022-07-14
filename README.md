# EventPlaneROOT

Analysis code for calculating event plane resolution for FIT detectors (FV0 and FT0) from O2 simulation data.

## Usage:

 * ./saveQvecs \<output.root\> \<path to input\> \<some label\> \<docorrections\>
 * ./calcCorrections \<output.root\> \<path to input\> \<some label\>
 * ./calcEventPlaneRes \<output.root\> \<path to input\> \<some label\> \<docorrections\>

## Comments
 * Basic workflow from the scratch has been the following:
    * Calculate the Q-vectors without any corrections (for example <code>./saveQvecs qvecs_no-corr.root \path\to\indput\dir 10-20 0</code>)
        * In Puhti the input directory had several job directories (<code>PbPb_StandaloneAMPT_o2ver-21-10-12_cent10-20_5500GeV/dig_5Hz-mips8/run_jobY</code>) and <code>saveQvecs</code> reads it correctly if the root directory for the run (<code>PbPb_StandaloneAMPT_o2ver-21-10-12_cent10-20_5500GeV</code> in this case) is given. So if other path structure is used make necessary changes to <code>saveQvecs</code> script
    * Calculate corrections from the uncorrected Q-vectors saved into <code>qvecs_no-corr.root</code> (<code>./calcCorrections qvecs_no-corr.root \path\to\indput\dir 10-20</code>)
        * Corrections are saved to individual text file for each of the three detectors
    * Then calculate Q-vectors again but now with corrections. This is done automaticly when using script <code>calcEventPlaneRes</code> (<code>./calcEventPlaneRes qvecs_with-corr.root \path\to\indput\dir 10-20 1</code>)
        *  runs first saveQvecs and saves Q-vectors to output.root, output for resolutions is labeled as "res_qvecs_with-corr.root"
 * <code>calc-res-all.sh</code> and <code>calc-corr-all.sh</code> I've been using to run the code over all centralities at once
 * Data for one 100 event run (centrality 10-20) can be found from <code>data</code> directory
