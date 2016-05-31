%  ---------------------------------------------------------------------
%  Copyright (C) 2016 by the LearnEDFM authors
% 
%  This file is part of LearnEDFM.
% 
%  LearnEDFM is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  LearnEDFM is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with LearnEDFM.  If not, see <http://www.gnu.org/licenses/>.
%  ---------------------------------------------------------------------
% 
%  Authors: Gunnar Jansen, University of Neuchatel, 2016 
%           Ivan Lunati, Rouven Kuenze, University of Lausanne, 2012

%% GRID PARAMETERS ---------------------------------------------------------------------------------%
global Nf Nf_f len dx
len     = [9 9];                    % physical length of the domain in x and y direction [m]
Nf      = [99 99];  	                % number of cells in x and y direction
dx      = len./Nf;                      % cell length [m]

%% FRACTURE NETWORK
% Crossed fracture test case
frac_cross
% Random fractures 
%frac_rand
% Single fracture
%frac_single
% Crossed fractures (rotated 45°)
%frac_cross_rot45
% Thirteen 'random' fractures
%frac_complex_n13;

if (dxf < min(dx))
    error('dxf < dx')
end

%% PERMEABILITY ------------------------------------------------------------------------------------%
global K K_f
K       = ones(Nf(1),Nf(2))*1e-9; 	    % permeability field [m2]
K_f     = ones(Nf_f,1)*1e-6; 	        % fracture permeability field [m2]

%% Porosity ----------------------------------------------------------------------------------------%
global phi phi_f
phi     = ones(Nf(1),Nf(2))*0.3;	    % porosity field
phi_f   = ones(Nf_f,1)*0.3;             % fracture porosity field

%% INITIAL CONDITIONS--------------------------------------------------------------------------------%
global s0 s0f smax   
s0   = zeros(Nf(1),Nf(2));              % Initial saturation (normalized concentration) [-]
s0f  = zeros(Nf_f,1);
smax = 1;                               % maximum concentration [kg/m3] for normalization

global T0 T0f
T0   = zeros(Nf(1),Nf(2));              % Initial matrix temperature [°C]
T0f  = zeros(Nf_f,1);                   % Initial fracture temperature [°C]

%% PHASE PROPERTIES---------------------------------------------------------------------------------%
global viscosity density gravity
viscosity = [0.01 0.01];                % viscosity [kg/m/s] [viscosity(s=0) viscosity(s=smax)]
density   = [935 992];                  % density [kg/m3] [density(s=0) density(s=smax)]
gravity   = 0;                          % gravity acceleration in y [m/s2]

%% SIMULATION PARAMETER FOR TRANSPORT --------------------------------------------------------------%
global dt
timeSim  = 100;                        % total simulation time [s]
dt       = 100;                         % time step length [s]
tolpS    = 1.e-4;                       % saturation tolerance on pressure-saturation loop [-]
maxpS    = 100;                         % maximum number of pressure saturation loops to converge

%% BC FLUID ----------------------------------------------------------------------------------------%
global Fix ibcs                                
ibcs = zeros(2*sum(Nf),1);              % type 0:Neumann(N); 1:Dirichlet(D)
Fix  = zeros(2*sum(Nf),1);              % value N [m2/s] (inflow>0); D [Pa]   

ibcs(1:Nf(2)) = 1;
ibcs(Nf(2)+1:2*Nf(2))=1;
Fix(1:Nf(2)) = 1;
Fix(Nf(2)+1:2*Nf(2))=0;

%% BC TRANSPORT ------------------------------------------------------------------------------------%
global FixT
FixT     = zeros(2*sum(Nf),1);           % normalized concentration of boundary flow [-]
FixT(1:Nf(2)) = 1;

%% DIFFUSION ------------------------------------------------------------------------%
global Dif ibcD
Dif     = 0;                            % [m2/s] molecular diffusion 
ibcD    = zeros(2*sum(Nf),1);           % 1 -> Diffusion on boundary cell

%% SOURCE TERMS ------------------------------------------------------------------------------------%
global Q QT
Q       = zeros(Nf);                    % source term [m2/s]; inflow positive
QT      = zeros(Nf);                    % normalized concentration for source term [-] 