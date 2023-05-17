# MATLAB code: "Joint Distribution of Distance and Angles in Finite Wireless Networks"
**Corresponding author**: Francisco J. Martín Vega (email: fjmv@uma.es, fjmvega@ic.uma.es)
Communications and Signal Processing group (ComSP),  Telecommunication Research Institute (TELMA), E.T.S.I.T., University of Malaga (Spain).

## Introduction
This repository contains the code validates through Monte Carlo simulation and evaluates numerrically the expressions derived in the following paper:

Francisco J. Martin-Vega, Gerardo Gomez, David Morales-Jimenez, F. Javier Lopez-Martinez and Mari Carmen Aguayo-Torres, "Joint  Distribution of Distance and Angles in Finite Wireless Networks", acepted for publication in IEEE Transactions on Vehicular Technology, 2023.

A pre-print of such a paper is available at:
https://arxiv.org/abs/2203.13510

The main results include the joint distribuition of distance and angles (azimuth and zenith) between a reference node, wich is placed at an arbritary location $u=(u_x,u_y,u_z)$ ,  and randomly distributed nodes that falls within a finite area rectangular region. To illustrate the applications of these expressions, the distrubuition of the SNR with a directional antenna pattern is obtained. 

## Abstract
Directional beamforming will play a paramount role in 5G and beyond networks to combat the higher path losses incurred at millimeter wave bands. Appropriate modeling and analysis of the angles and distances between transmitters and receivers in these networks are thus essential to understand performance and limiting factors. Most existing literature considers either infinite and uniform networks, where nodes are drawn according to a Poisson point process, or finite networks with the reference receiver placed at the origin of a disk. Under either of these assumptions, the distance and azimuth angle between transmitter and receiver are independent, and the angle follows a uniform distribution between 0 and 2𝜋. Here, we consider a more realistic case of finite networks where the reference node is placed at any arbitrary location. We obtain the joint distribution between the distance and azimuth angle and demonstrate that these random variables do exhibit certain correlation, which depends on the shape of the region and the location of the reference node. To conduct the analysis, we present a general mathematical framework that is specialized to exemplify the case of a rectangular region. We also derive the statistics for the 3D case where, considering antenna heights, the joint distribution of distance, azimuth, and zenith angles is obtained. Finally, we describe some immediate applications of the present work, including the design of analog codebooks, wireless routing algorithms, and the analysis of directional beamforming, which is illustrated by analyzing the coverage probability of an indoor scenario considering misaligned beams.

## MATLAB version:
This code has been tested with MATLAB R2022b.

## How to start
### Step 1:
To run a simulation that computes the joint distributions of distance and angle and their marginals for the 2D or 3D cases, you can modify the scripts `scr_joint_distribuitions_2D.m` or `scr_joint_distribuitions_3D.m` respectively.  You will see a number of figures that are drawn comparing theoretical versus simulation results. 
### Step 2:
To get the CCDF of the SNR and the p-percentile of the SNR for the case of a directive antenna pattern, you can modify and run the script `scr_runSimulation.m`. Once this script finishes its execution, it creates a results folder, in the same path where the script is executed, that contains the figures that has been created. 

## Hints about the code & folder structure
The folder `Code\Expressions` contains the MATLAB code that evaluates the main theoretical results. Of particlar interest are the following functions: 
- `ccdfSNRfun.m`: it evaluates **Corrolary 9** and **Approximation 2** of the related paper.
- `anonymousFunctions.m`: it evaluates the expressions given with eq. (22) and those given with **Corollary 4**.

The folder `SampleScripts` contains some scripts that reproduces the results included in the related paper. 

The function `getNumResults.m` allows performing a parameter sweep to draw figures where a single metric, e.g., the 10%-percentile, is obtained for different values of a given parameter that is swept, .e.g, the x-coordinate of the reference node location. All the parameters that can be swept are stored in the struct, `vP`, whereas the parameters that cannot be swept are stored in the struct `p`. 