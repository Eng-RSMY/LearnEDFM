function EDFM(inputFile,showPlot)
% LearnEDFM Numerical code to solve flow and tracer transport in fractured porous media
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
%
%  Acknowledgement:  thanks to Manav Tyagi, Brad Mallison and Hadi Hajibeygi
%                    for contributing to the early development of the code. 
%
%  LearnEDFM: Numerical code to solve flow and tracer transport in
%             fractured porous media using the embdedded discrete fracture model
%
%  LearnEDFM(inputFile,showPlot)
%
%  input    ->  must be a string (optional, default 'inputFile')
%  showPlot ->  must be a boolean(optional, default '1')
%
%  LearnEDFM is the main program "driving" the simulation
%  It contains the time loop and calls the following functions:
%  
%  > initialize
%  > initializeDFN     >> intersectionsGrid
%                            >>> calc_d_mean
%                      >> intersectionsSegments
%  > pressureSystem
%  > calcVelocity
%  > transportSystem   >> transportAdvection
%                      >> transportDiffusion
%
%  to see the structure of the input file, open InputFile.m
% 
%%-------------------------------------------------------------------------%
%
%  Conventions used in the code:
%
%  1) fluxes:  outgoing fluxes are negative; ingoing fluxes are positive
%  2) gravity: positive g -> gravity is directed downwards
%  3) scalar fields: p(x,y), s(x,y), Kx(x,y), Ky(x,y), Q(x,y), QT(x,y), ...
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%                             --------------------------
%                             |(1,3)  (2,3) (3,3) (4,3)|
%                          Y  |(1,2)  (2,2) (3,2) (4,2)|
%                             |(1,1)  (2,1) (3,1) (4,1)|
%                             --------------------------
%                                          X
%
%
%  4) boundary conditions (b.c.):
%
%                   ibcs(i) and ibcT(i) define the b.c. type:
%                   1 -> Dirichlet; 0 -> Neumann 
%
%                   Fix(i) and FixT(i) contain the assigned value:
%                   if ibcs(i) = 1 -> Fix(i) is the pressure [Pa]
%                   if ibcs(i) = 0 -> Fix(i) is the flux (per unit 
%                                         of transversal length) [m2/s]
%                   if ibcT(i) = 1 -> FixT(i) is the concentration [kg/m3]
%
%                   the boundary conditions are assigned by a vector that
%                   represents all the cells on the perimeter.
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%
%                              i= 11  i= 12 i= 13 i= 14
%                             --------------------------
%                       i=  3 |(1,3)  (2,3) (3,3) (4,3)| i= 6
%                       i=  2 |(1,2)  (2,2) (3,2) (4,2)| i= 5
%                       i=  1 |(1,1)  (2,1) (3,1) (4,1)| i= 4
%                             --------------------------  
%                              i= 7   i= 8  i=  9 i= 10
%
%  5) vector fields: vx(x,y), vy(x,y), Kx(x,y), Ky(x,y), ...
%                                     [Kx,Ky are diagonal matrix, they can 
%                                              be represented as vector...]
%
%                   vectors are defined at cell interfaces, e.g.,
%                   
%                   vx(i,j) is the velocity between (i,j) and (i+1,j)
%                   vy(i,j) is the velocity between (i,j) and (i,j+1)
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%                   vx:
%                             --------------------------
%                           (1,3)  (2,3) (3,3) (4,3) (5,4)
%                        Y  (1,2)  (2,2) (3,2) (4,2) (5,2)
%                           (1,1)  (2,1) (3,1) (4,1) (5,1)
%                             --------------------------
%                                          X
%                   vy:
%                             |(1,4)  (2,4) (3,4) (4,4)|
%                             |(1,3)  (2,3) (3,3) (4,3)|
%                          Y  |(1,2)  (2,2) (3,2) (4,2)|
%                             |(1,1)  (2,1) (3,1) (4,1)|
%                                          X
%
%  6) linear algebra:
%                   for the solution of the linear system, the scalar
%                   fields are transformed into vectors following the
%                   convention:
%
%
%                   Examples for a 4 x 3 grid (Nf = [4 3])
%
%                             --------------------------
%                             |  9     10    11    12  |
%                          Y  |  5      6     7     8  |
%                             |  1      2     3     4  |
%                             --------------------------
%                                          X
%
%                   thus the vector is
%                                       i= 1 |(1,1)|
%                                       i= 2 |(2,1)|
%                                       i= 3 |(3,1)|
%                                       i= 4 |(4,1)|
%                                       i= 5 |(1,2)|
%                                       i= 6 |(2,2)|
%                                       i= 7 |(3,2)|
%                                       i= 8 |(4,2)|
%                                       i= 9 |(1,3)|
%                                       i=10 |(2,3)|
%                                       i=11 |(3,3)|
%                                       i=12 |(4,3)|
%
%
%-------------------------------------------------------------------------%

% Add path (at beginning of script)
addpath(genpath(pwd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULTS AND GLOBAL PARAMETERS                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Nf len Fix Q Txk Tyk Tfk dt XY1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         INITIALIZATION                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<1)
    warning('No input file given. Using default: InputFile.')
    inputFile  = 'InputFile'; 
end
if (nargin<2), showPlot = 1; end

eval(inputFile)                                                            % Load Input Data 

x = linspace(0,len(1),2*Nf(1)+1); x = x(2:2:end);                          % x grid vector (used in plots)
y = linspace(0,len(2),2*Nf(2)+1); y = y(2:2:end);                          % y grid vector (used in plots)

TNew  = T0  .* ones(Nf);                                                   % initialize matrix   temperature solution vector
TNewf = T0f .* ones(Nf_f,1);                                               % initialize fracture temperature solution vector

sNew  = s0  .* ones(Nf);                                                   % initialize matrix   tracer saturation solution vector
sNewf = s0f .* ones(Nf_f,1);                                               % initialize fracture tracer saturation solution vector
sN    = [sNew(:); sNewf(:)];

% Calculate matrix interface permeability & porosity (harmonic average)
Txk             = sparse(Nf(1)+1,Nf(2));                                                                                                                
Tyk             = sparse(Nf(1),Nf(2)+1);
Txk(2:Nf(1),:)  = 2./(1./K(1:Nf(1)-1,:) + 1./K(2:Nf(1),:));             
Tyk(:,2:Nf(2))  = 2./(1./K(:,1:Nf(2)-1) + 1./K(:,2:Nf(2)));
Txk(1,:)        = K(1,:);                                    
Txk(Nf(1)+1,:)  = K(Nf(1),:);
Tyk(:,1)        = K(:,1);
Tyk(:,Nf(2)+1)  = K(:,Nf(2));

phix            = sparse(Nf(1)+1,Nf(2));
phiy            = sparse(Nf(1),Nf(2)+1);
phix(2:Nf(1),:) = 2./(1./phi(1:Nf(1)-1,:) + 1./phi(2:Nf(1),:));           
phiy(:,2:Nf(2)) = 2./(1./phi(:,1:Nf(2)-1) + 1./phi(:,2:Nf(2)));
phix(1,:)       = phi(1,:);                                    
phix(Nf(1)+1,:) = phi(Nf(1),:);
phiy(:,1)       = phi(:,1);
phiy(:,Nf(2)+1) = phi(:,Nf(2));

% Calculate fracture interface permeability % porosity (harmonic average)
% Each fracture has one more interface than its fractures elements which
% raises an issue with the harmonic averaging. It can be resolved by an 
% index shift by 1 has to be done in each iteration leading to an
% overall dimension for Tfk of (Nf_f + N_fractures,1). This is of course
% true also for the velocity, saturation and porosity fields.
% Note that any changes here have also to be transfered to the other loops
% over the fractures.
Tfk = zeros(Nf_f+N_fractures,1);
phif = zeros(Nf_f+N_fractures,1);
ios = 0;                                                                   % offset counter for positioning in global fracture vector
lios = 0;                                                                  % offset counter for the left side
for i= 1:N_fractures 
    Tfk(lios+2:lios+Nf_i(i))     = 2./(1./K_f(ios+1:ios+Nf_i(i)-1) + 1./K_f(ios+2:ios+Nf_i(i))); 
    Tfk(lios+1)          = K_f(ios+1);
    Tfk(lios+Nf_i(i)+1)     = K_f(ios+Nf_i(i));
    
    phif(lios+2:lios+Nf_i(i))     = 2./(1./phif(ios+1:ios+Nf_i(i)-1) + 1./phif(ios+2:ios+Nf_i(i))); 
    phif(lios+1)          = phif(ios+1);
    phif(lios+Nf_i(i)+1)     = phif(ios+Nf_i(i));
    
    ios = ios + Nf_i(i);
    lios = lios + Nf_i(i) + 1;                                             % the additional offset of +1 is due to the additional interface at the end of each fracture
end

% Initialize the DFN connectivity and grid intersections
% This provides also the off diagonal coupling matrices for the pressure
% solver. 
[~,~,~,~,~,~,~,l,l_f]  = initialize(sNew,sNewf,TNew,TNewf);                % Initialize the matrix transmissivity        
[Tfm,Tmf,Tff,Dfm,Dmf,Z] = initializeDFN(x,y,XY1,K_f,l,l_f);

if showPlot
    s = linspace(0,1,256);
    rgb1=[0.230, 0.299, 0.754];
    rgb2=[0.706, 0.016, 0.150];
    cmap = diverging_map(s,rgb1,rgb2);

    figure(1)
    h1= pcolor(x,y,s0');
    hold on
    colormap(cmap)
    line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r');
    shading interp
    axis square
    
    figure(2)
    h2 = pcolor(x,y,sNew');
    hold on
    colormap(cmap)
    shading interp
    axis square
    caxis([0 1]);
    line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r');
    
    drawnow

    nframe=ceil(timeSim/dt);
    mov(1:nframe)= struct('cdata',[],'colormap',[]);
    set(gca,'nextplot','replacechildren')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         SIMULATION                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = dt;                                                                 % Set the time for the first timestep                                             
i = 0;                                                                     % Count all timesteps

while (time <= timeSim)                                                    % Time loop
    i      = i+1;
    if showPlot, disp(sprintf('Simulation time = %d s',time));end
    sOld      = sNew;                                                      % Update saturation of previous iteration (time loop)
    sOldf     = sNewf;
    eps       = inf;
    innerIter = 1;

    while (eps >= tolpS && innerIter <= maxpS)                             % Inner loop due to saturation dependence of gravity and viscosity

        sIt = sN;                                                          % Update saturation of previous iteration (inner loop)

        [Tx,Ty,Tf,G,g,Gf,gf,~,~]  = initialize(sNew,sNewf,TNew,TNewf);     % Initialize (update) the transmissivity based on new temperature
                                                                           % If the viscosity and density are not temperature dependent, this
                                                                           % function could be moved outside of the timeloop
 
        [A,rhs]= pressureSystem(Fix,ibcs,Tx,Ty,Tf,G,Gf,Q, Tfm, Tmf, Tff, Nf_i); % Construct pressure matrix and rhs vectors

        u = A\rhs;                                                         % Solve the pressure equation

        p = reshape(u(1:prod(Nf)),Nf(1),Nf(2));                            % Assign solution to matrix solution array 
        pf = u(prod(Nf)+1:length(u));                                      % Assign solution to fracture solution vector
        
        
        [vx,vy,vf,Vfm,Vmf,Vff] = calcVelocity(p,pf,Tx,Ty,Tf,Tfm,Tmf,Tff,g,gf); % Calculate Darcy velocities     

        [sN]  = transportSystem(sOld,sOldf,vx,vy,vf,Vfm,Vmf,Vff,Dfm,Dmf,phix,phiy); % Solve transport equation for phase alpha

        sNew = sN(1:prod(Nf));                                             % Assign transport solution to matrix solution array 
        sNew= reshape(sNew,Nf(1),Nf(2));                         
        sNewf = sN(prod(Nf)+1:length(sN));                                 % Assign transport solution to fracture solution array 
        
        eps = norm((abs(sN(:) - sIt(:))),inf); 
        if showPlot
            disp(sprintf('\t Residual at %d. loop: %d', innerIter,eps));       
        end
        innerIter = innerIter+1;


        if (innerIter == 100), error('outer finescale loop did not converge'), end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             Output             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if showPlot
        set(h1,'cdata',p');
        set(h2,'cdata',sNew');
        drawnow
        mov(i)=getframe(gcf);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%     Time Step Control                                             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (time+dt > timeSim) 
    dt = timeSim - time;
    if dt == 0; dt = nan; end
    end

    time=time+dt;
end

if showPlot
    try
        movie2avi(mov, 'output.avi', 'compression', 'None');
    catch
        warning('Problem generating video output. Do not touch the figure windows during simulations for video output!');
    end
end

% Since EDFM is a function, the variables used are not written into the
% workspace after the simulation has ended. This can be achieved by the
% following if relevant.
ListOfVariables = who;
for k = 1:length(ListOfVariables)
   assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
end 