% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script adds paths to code folders. 
% OUTPUTs:
% +) currentDir: It is the set of working paths including the added paths
% +) oldPath: It is the set of working paths before adding those of the
% project. 
% *************************************************************************

function [currentDir, oldPath] = addPaths()

oldPath = path;
currentDir = pwd;
idcs   = strfind(currentDir,'\');
parentDir = currentDir(1:idcs(end)-1);

path(parentDir, path);
path([parentDir '\Code\Utilities'    ], path); 
path([parentDir '\Code\Statistics'   ], path);
path([parentDir '\Code\Simulation'   ], path);
path([parentDir '\Code\Expressions'   ], path)
cd(currentDir);