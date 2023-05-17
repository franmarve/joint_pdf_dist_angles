% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin-Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This function translates angles defined within the range from -pi to pi 
% to the range from 0 to 2pi
% *************************************************************************

function theta = angle_0_to_2pi(theta)

theta(theta < 0) = theta(theta < 0) + 2*pi;