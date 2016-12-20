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

global Nf len dt innerIter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         INITIALIZATION                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<1)
    warning('No input file given. Using default: InputFile.')
    inputFile  = 'InputFile'; 
end
if (nargin<2), showPlot = 1; end

run(inputFile)                                                            % Load Input Data 

x = linspace(0,len(1),2*Nf(1)+1); x = x(2:2:end);                          % x grid vector (used in plots)
y = linspace(0,len(2),2*Nf(2)+1); y = y(2:2:end);                          % y grid vector (used in plots)

global frac_mean_mat frac_grad_mat frac_grad_expand_mat
frac_mean_mat = calc_frac_mean_mat(N_fractures,Nf_f, Nf_i);
frac_grad_mat = calc_frac_grad_mat(N_fractures,Nf_f, Nf_i);
frac_grad_expand_mat = calc_frac_grad_expand_mat(N_fractures,Nf_f, Nf_i);


tNew  = T0  .* ones(Nf);                                                   % initialize matrix   temperature solution vector
tNewf = T0f .* ones(Nf_f,1);                                               % initialize fracture temperature solution vector
tN    = [tNew(:); tNewf(:)];

cNew  = c0  .* ones(Nf);                                                   % initialize matrix   tracer saturation solution vector
cNewf = c0f .* ones(Nf_f,1);                                               % initialize fracture tracer saturation solution vector
cN    = [cNew(:); cNewf(:)];

p  = p0  .* ones(Nf);                                                   % initialize matrix   tracer saturation solution vector
pf = p0f .* ones(Nf_f,1);                                               % initialize fracture tracer saturation solution vector
pN    = [p(:); pf(:)];

%% Find the intersection of each segment with the grid and compute the
%% Connectivity Index based on the segments lengths and the average distance
%% to the fracture
CI = intersectionsGrid(x,y,XY1);
[row,col] = intersectionsSegments(XY1,XY1);

[Tx,Ty,Tf,Tfm,Tff,g,gf,density_l,density_lf]  = initialize(K, K_f,cNew,cNewf,tNew,tNewf,p,pf,CI,row,col);

if showPlot
    s = linspace(0,1,256);
    rgb1=[0.230, 0.299, 0.754];
    rgb2=[0.706, 0.016, 0.150];
    cmap = diverging_map(s,rgb1,rgb2);

    figure(1)
    h1= pcolor(x,y,c0');
    hold on
    colormap(cmap)
    line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r');
    shading interp
    axis square
    
    if(flagTracerTransport)
        figure(2)
        h2 = pcolor(x,y,cNew');
        hold on
        colormap(cmap)
        shading interp
        axis square
        caxis([0 1]);
        line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','k','LineWidth',2);
        h3 = scatter(xe,ye,30,cNewf,'filled');
    end
    
    if(flagHeatTransport)
        figure(3)
        h4 = pcolor(x,y,tNew');
        hold on
        colormap(cmap)
        shading interp
        axis square
        caxis([0 tmax]); 
        line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','k','LineWidth',2);
        h5 = scatter(xe,ye,30,tNewf,'filled');
        drawnow
    end

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
    cOld      = cNew;                                                      % Update saturation of previous iteration (time loop)
    cOldf     = cNewf;
    tOld      = tNew;                                                      % Update saturation of previous iteration (time loop)
    tOldf     = tNewf;
    
    pOld      = p;                                                      % Update saturation of previous iteration (time loop)
    pOldf     = pf;
    
    epsP       = inf;
    epsC       = inf;
    epsT       = inf;
        
    innerIter = 1;
    while ((epsP >= tol || epsC >= tol || epsT >= tol) && innerIter <= maxit)                            % Inner loop due to saturation dependence of gravity and viscosity
        cIt = cN;                                                          % Update saturation of previous iteration (inner loop)
        tIt = tN;                                                          % Update saturation of previous iteration (inner loop)
        pIt = pN;
        
        [Tx,Ty,Tf,Tfm,Tff, g,gf,density_l,density_lf]  = initialize(K, K_f,cNew,cNewf,tNew,tNewf,p,pf,CI,row,col);
                                                                           % Initialize (update) the transmissivity based on new temperature
                                                                           % If the viscosity and density are not temperature dependent, this
                                                                           % function could be moved outside of the timeloop
                                                                           
        [pN]= pressureSystem(pOld,pOldf,Fix,ibcs,Tx,Ty,Tf,g,gf,Q,Tfm,Tff,Nf_i,phi,phi_f,density_l,density_lf,K_f); % Construct pressure matrix and rhs vectors

        
        p = reshape(pN(1:prod(Nf)),Nf(1),Nf(2));                            % Assign solution to matrix solution array 
        pf = pN(prod(Nf)+1:length(pN));                                      % Assign solution to fracture solution vector
        
        [vx,vy,vf,Vfm,Vmf,Vff] = calcVelocity(p,pf,Tx,Ty,Tf,Tfm,Tff,g,gf); % Calculate Darcy velocities     
        
        if(flagTracerTransport)
            [cN]  = transport_mass_System(cOld,cOldf,vx,vy,vf,Vfm,Vmf,Vff,Q,QC,phi,phi_f,K_f); % Solve transport equation for phase alpha

            cNew = cN(1:prod(Nf));                                             % Assign transport solution to matrix solution array 
            cNew= reshape(cNew,Nf(1),Nf(2));                         
            cNewf = cN(prod(Nf)+1:length(cN));                                 % Assign transport solution to fracture solution array 
        end
        
        if(flagHeatTransport)
            [tN]  = transport_heat_System(tOld,tOldf,vx,vy,vf,Vfm,Vmf,Vff,Q,QT,K_f,phi,phi_f,density_s,density_sf,density_l,density_lf); % Solve transport equation for phase alpha

            tNew = tN(1:prod(Nf));                     

            % Assign transport solution to matrix solution array 
            tNew= reshape(tNew,Nf(1),Nf(2));                         
            tNewf = tN(prod(Nf)+1:length(tN));                                 % Assign transport solution to fracture solution array 
        end
        
        epsP = norm((abs(pN(:) - pIt(:))),inf); 
        epsC = norm((abs(cN(:) - cIt(:))),inf); 
        epsT = norm((abs(tN(:) - tIt(:))),inf); 
        
        if showPlot
            disp(sprintf('\t Residual at %d. loop: %d , %d and %d', innerIter,epsP, epsC, epsT));       
        end
        innerIter = innerIter+1;


        if (innerIter == maxit), error('outer finescale loop did not converge'), end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             Output             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if showPlot
        set(h1,'cdata',p');
        if(flagTracerTransport)
            set(h2,'cdata',cNew');
            set(h3,'cdata',cNewf);
        end
        if(flagHeatTransport)
            set(h4,'cdata',tNew');
            set(h5,'cdata',tNewf);
        end
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
