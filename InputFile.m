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
len     = [500 500];                    % physical length of the domain in x and y direction [m]
Nf      = [45 45];  	                % number of cells in x and y direction
dx      = len./Nf;                      % cell length [m]

%% SIMULATION PARAMETER FOR TRANSPORT --------------------------------------------------------------%
global dt
timeSim  = 5e5;                         % total simulation time [s]
dt       = 1000;                           % time step length [s]
tol    = 1.e-4;                       % saturation tolerance on pressure-concentration-heat loop [-]
maxit    = 100;                         % maximum number of pressure concentration-heat loops to converge

%% FRACTURE NETWORK
% Crossed fracture test case
%frac_cross
% Random fractures 
%frac_rand
% Single fracture
%frac_single
% Crossed fractures (rotated 45°)
%frac_cross_rot45
% Thirteen 'random' fractures
frac_complex_n13;
% Leeds 4a granite fracture data
%frac_granite_leeds4a;

if (dxf < min(dx))
    error('dxf < dx')
end

%% INITIAL CONDITIONS--------------------------------------------------------------------------------%
global cmax   
c0   = zeros(Nf(1),Nf(2));              % Initial saturation (normalized concentration) [-]
c0f  = zeros(Nf_f,1);
cmax = 1;                               % maximum concentration [kg/m3] for normalization

T0   = zeros(Nf(1),Nf(2));              % Initial matrix temperature [°C]
T0f  = zeros(Nf_f,1);                   % Initial fracture temperature [°C]
tmax = 10;                              % maximum temperature [°C] for plotting

p0   = zeros(Nf(1),Nf(2));              % Initial matrix pressure [Pa]
p0f  = zeros(Nf_f,1);                   % Initial fracture pressure [Pa]

%% BC FLUID ----------------------------------------------------------------------------------------%
global Fix ibcs                                
ibcs = zeros(2*sum(Nf),1);              % type 0:Neumann(N); 1:Dirichlet(D)
Fix  = zeros(2*sum(Nf),1);              % value N [m2/s] (inflow>0); D [Pa]   

% ibcs(1:Nf(2)) = 1;
% ibcs(Nf(2)+1:2*Nf(2))=1;
% Fix(1:Nf(2)) = 1e5;
% Fix(Nf(2)+1:2*Nf(2))=0;
ibcs(1) = 1;
ibcs(2*Nf(1)+1)=1;
ibcs(2*Nf(2))=1;
ibcs(2*Nf(1)+2*Nf(2))=1;
Fix(1) = 1e7;
Fix(2*Nf(2))=0;
Fix(2*Nf(1)+1) = 1e7;
Fix(2*Nf(1)+2*Nf(2))=0;

%% BC TRANSPORT ------------------------------------------------------------------------------------%
flagTracerTransport = 0;
flagHeatTransport = 0;

global FixT FixC
FixT     = zeros(2*sum(Nf),1);           % normalized concentration of boundary flow [-]
FixC     = zeros(2*sum(Nf),1);           % normalized concentration of boundary flow [-]

%FixT(1:Nf(2)) = 1;
FixT(1) = 1;
FixT(2*Nf(1)+1) = 1;

%% SOURCE TERMS ------------------------------------------------------------------------------------%
Q       = zeros(Nf);                    % source term [m2/s]; inflow positive
QC      = zeros(Nf);                    % normalized concentration for source term [-] 
QT      = zeros(Nf);                    % normalized concentration for source term [-] 

%% GRAVITY---------------------------------------------------------------------------------%
global gravity
gravity   = 0;                          % gravity acceleration in y [m/s2]

%% PERMEABILITY ------------------------------------------------------------------------------------%
K       = ones(Nf(1),Nf(2))*1e-9; 	    % permeability field [m2]
K_f     = ones(Nf_f,1)*1e-7; 	        % fracture permeability field [m2]

%% Porosity ----------------------------------------------------------------------------------------%
phi     = ones(Nf(1),Nf(2))*0.3;	    % porosity field
phi_f   = ones(Nf_f,1)*0.3;             % fracture porosity field

%% MECHANIC PROPERTIES
global k_modulus_l k_modulus_s biot_alpha
k_modulus_l = 2.15e9;                        % Bulk Modulus of the fluid [Pa]
k_modulus_s = 1.8e10;                        % Bulk Modulus of the rock  [Pa]
biot_alpha = 1;

global shear_modulus poisson_ratio
shear_modulus = 29e9;                        % Shear modulus of the rock [Pa]
poisson_ratio = 0.25;                        % Poisson ratio of the rock [-]

global therm_exp_coeff
therm_exp_coeff = 7.9e-6;                    % Thermal expansion coefficient of the rock matrix [-]

%% ROCK DENSITY ------------------------------------------------------------------------%
density_s = 2000*ones(Nf);       % density of the rock [kg/m3]
density_sf = 2000*ones(Nf_f,1); % density of the rock [kg/m3]

%% IN SITU STRESS ------------------------------------------------------------------------%
SH_max = 52e6*ones(Nf_f,1);                      % Maximum principal stress [Pa]
SH_min = 20e6*ones(Nf_f,1);                      % Minimum principal stress [Pa]

%% THERMAL DIFFUSION ------------------------------------------------------------------------%
global lambda_l lambda_s cp_l cp_s ibcD
lambda_l = 0;%0.5;                           % Thermal conductivity of the fluid [W/(m*K)]
lambda_s = 0;%2.0;                           % Thermal conductivity of the rock [W/(m*K)]
cp_l = 0;%1100;                               % Specific heat capacity of the fluid [J/(kg*K)]
cp_s = 0;%250;                               % Specific heat capacity of the rock [J/(kg*K)]
ibcD    = zeros(2*sum(Nf),1);           % 1 -> Diffusion on boundary cells

%% MOLECULAR DIFFUSION ------------------------------------------------------------------------%
global DifC ibcDC
DifC     = 0;                  % [m2/s] molecular diffusion 
ibcDC    = zeros(2*sum(Nf),1);           % 1 -> Diffusion on boundary cells

%% MASS DISPERSION ------------------------------------------------------------------------%
global alphal alphat
alphal   = 0.0;                         % longitudinal dispersivity [m]
alphat   = 0.0;                          % transversal dispersivity [m]