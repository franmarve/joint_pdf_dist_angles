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
currentDir = cd;
path([currentDir '\Code\Utilities'      ], path); 
path([currentDir '\Code\Statistics'     ], path);
path([currentDir '\Code\Simulation'     ], path);
path([currentDir '\Code\Expressions'    ], path);
cd(currentDir);