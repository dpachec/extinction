# System requirements

This code has been implemented in Matlab 2024b, using EEGLAB (2021.0) and Fieldtrip (2023-12-20). 

# Instructions

To run the scripts, open Matlab and add the additional_functions folder to the path 

Use behavior.m to reproduce the behavioral results (Fig. 1). To generate _trlinfo.mat files for every subject, containing matlab matrices that store the parsed logs, run log2trialinfo.m. 

Use power_analysis to reproduce the theta power results (Fig. 2) 

Use contrast_based_RSA.m to reproduce the main effects of item stability and item specificity (Figs. 3-4). To generate neural RSMs, use the script create_neural_RSMs.m.

Use trial_level_analysis.m to perform the coordination analyses in all ROIs (Figs 3-5). To extract trial level estimates of item stability and context specificity, use extract_tr_lev_stab_cxt_spec.m







