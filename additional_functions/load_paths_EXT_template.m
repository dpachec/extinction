% This template file indicates placeholders for paths 
% Raw data and preprocessed data folders, and output results folders should be specified. 
% This template script and be renamed to load_paths_EXT.m to run with the other scripts in the repository

function [paths] = load_paths()


paths.results.rsa = '/path/to/project/results/rsa/';
paths.results.rsa_perm = '/path/to/project/results/rsa/perm/';
paths.trlinfo = '/path/to/data/preproc/trialinfo/'; %Store here all trlinfo files 
paths.results.behavior = '/path/to/data/behavior/';
paths.results.power = '/path/to/results/bipolar/power/norm_across_trials/50TR/';
paths.results.traces = '/path/to/results/bipolar/traces/';
paths.results.nRDMs = '/path/to/results/bipolar/neural_RDMs/'; 
paths.results.POWfromRT = '/path/to/results/bipolar/POWfromRT/'; 
paths.results.trial_based = '/path/to/results/bipolar/trial_based/'; 
paths.results.clusters = '/path/to/results/bipolar/clusters/';
