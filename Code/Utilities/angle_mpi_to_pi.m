% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function translates angles defined within the range from 0 to 2pi to
% the range from -pi to pi
% *************************************************************************

function theta = angle_mpi_to_pi(theta)

theta(theta > pi) = theta(theta > pi) - 2*pi;