function setup_study_paths()
%SETUP_STUDY_PATHS Add only package MATLAB and configuration folders.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'config'));
addpath(genpath(fullfile(root, 'matlab')));
end
