function setup

% SETUP    Add the directories required by this project to the MATLAB path.
%
%   Run this once per MATLAB session before calling roboticArm or circuit.

root = fileparts(mfilename('fullpath'));
toolbox = fullfile(root, 'external', 'DAEPreprocessingToolbox');

if ~exist(fullfile(toolbox, 'src', 'preprocessDAE.m'), 'file')
    error(['DAEPreprocessingToolbox not found. ' ...
           'Run "git submodule update --init" in the repository root.']);
end

addpath(fullfile(toolbox, 'src'));
addpath(fullfile(toolbox, 'test'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'problems'));
addpath(fullfile(root, 'experiments'));

% added last so that it precedes the toolbox on the search path
addpath(fullfile(root, 'overlay'));

clear loadMuPADPackage
loadMuPADPackage;
