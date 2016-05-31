function [Tx,Ty,Tf,Gr,g,Grf,gf, lambda, lambda_f] = initialize(S,Sf,T,Tf)
%  Initialize updates the transmissivity based on the saturation and
%  temperature
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
%  Initialize(S,Sf,T,Tf)
%
%  Input: 
%        S      (nx,ny)         matrix saturation
%        Sf     (nf,1)          fracture saturation
%        T      (nx,ny)         matrix temperature
%        Tf     (nf,1)          fracture temperature
%
%  Output:
%        Tx     (nx+1,ny)       trasmissibility in the x direction
%        Ty     (nx,ny+1)       trasmissibility in the y direction
%        Tf     (nf,1)          fracture transmissivity
%        Gr     (nx,ny+1)       gravity term corresponding to fine cells
%        g      (nx,ny)         gravity term with respect to the interfaces
%        Grf    (nf+n_fractures,1) gravity term corresponding to fracture cells
%        gf     (nf,1)          gravity term with respect to the fracture interfaces
%        lambda (nx,ny)         total mobility of fine grid cells
%        lambda_f (nf,1)        total mobility of the fracture segments

global Txk Tyk dx viscosity FixT ibcs gravity density smax Nf_f Nf_i dxf Tfk K K_f frac_angle
                
N = size(S);

%-------------------------------------------------------------------------%
%    Initialize                                                           %
%-------------------------------------------------------------------------%

i1 = 1:N(2);             i2  = N(2) + (1:N(2));                            
i3 = 2*N(2) + (1:N(1));  i4  = 2*N(2) + N(1) + (1:N(1));

bcwest(1:N(2))  = FixT(i1);
bceast(1:N(2))  = FixT(i2);
bcsouth(1:N(1)) = FixT(i3);
bcnorth(1:N(1)) = FixT(i4);

Satx           = sparse(N(1)+1,N(2));
Satx(2:N(1),:) = (S(2:N(1),:) + S(1:N(1)-1,:))./2./smax;
Satx(1,:)      = (S(1,:) + bcwest(1:N(2)))./2./smax;
Satx(N(1)+1,:) = (S(N(1),:) + bceast(1:N(2)))./2./smax;

Saty           = sparse(N(1),N(2)+1)./smax;
Saty(:,2:N(2)) = (S(:,2:N(2)) + S(:,1:N(2)-1))./2./smax;
Saty(:,1)      = (S(:,1) + bcsouth(1:N(1))')./2./smax;
Saty(:,N(2)+1) = (S(:,N(2)) + bcnorth(1:N(1))')./2./smax;

N_fractures = length(Nf_i);
Satf           = zeros(Nf_f+N_fractures,1);                                          
ios = 0;                                                                   % offset counter for positioning in global fracture vector
lios = 0;                                                                  % offset counter for the left side
for i= 1:N_fractures
    Satf(lios+2:lios+Nf_i(i)) = (Sf(ios+2:ios+Nf_i(i)) + Sf(ios+1:ios+Nf_i(i)-1))./2./smax; 
    Satf(lios+1)              = Sf(ios+1)./smax;
    Satf(lios+Nf_i(i)+1)      = Sf(ios+Nf_i(i))./smax;
    
    lios = lios + Nf_i(i) + 1;
    ios = ios + Nf_i(i);                                                   % Note the procedure on global indexing used here which is 
                                                                           % explained in the calculation of the fracture interface
                                                                           % permeabilities in the EDFM main file
end  

%-------------------------------------------------------------------------%
%    Calculate upwind Lamda total                                         %
%-------------------------------------------------------------------------%
lambda = K./((1-S).*viscosity(1) + S.*viscosity(2));
lambda_f = K_f./((1-Sf).*viscosity(1) + Sf.*viscosity(2));

mux   = (1-Satx).*viscosity(1) + Satx.*viscosity(2);
muy   = (1-Saty).*viscosity(1) + Saty.*viscosity(2);

muf   = (1-Satf).*viscosity(1) + Satf.*viscosity(2);

%-------------------------------------------------------------------------%
%    Gravity for the matrix                                               %
%-------------------------------------------------------------------------%

Gy    = ((1-Saty).*density(1) + Saty.*density(2)).*gravity;
g            = Tyk./muy.*Gy.*dx(1);                                        % Gravity term with respect to the interfaces
g(:,1)       = g(:,1).*ibcs(2*N(2)+1:2*N(2)+N(1)); 
g(:,N(2)+1)  = g(:,N(2)+1).*ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1)); 
Gr           = g(:,2:N(2)+1) - g(:,1:N(2));                                % Gravity term corresponding to fine cells

%-------------------------------------------------------------------------%
%    Gravity for the fractures                                            %
%-------------------------------------------------------------------------%
% The gravity computations on the fractures are a little bit more
% complicated. A couple of index shifts have to be performed throughout
% this file to achieve all the right dimensions of the requiered fields.
Gyf           = zeros(Nf_f+N_fractures,1);                                          
lios = 0;                                                                  % offset counter for positioning in global fracture vector
for i= 1:N_fractures
    
    Gyf(lios+1:lios+Nf_i(i)+1) = ((1-Satf(lios+1:lios+Nf_i(i)+1)).*density(1) ...
                                + Satf(lios+1:lios+Nf_i(i)+1).*density(2)).*gravity.*sin(frac_angle(i)*pi/180); 
                                                                           
    lios = lios + Nf_i(i) + 1;                                             % Note the procedure on global indexing used here which is 
                                                                           % explained in the calculation of the fracture interface
                                                                           % permeabilities in the EDFM main file
end  
         
gf           = Tfk./muf.*Gyf.*sqrt(12*Tfk);                                % Cubic law is used here instead of fixed aperture
Grf           = zeros(Nf_f,1);                                          
ios = 0;                                                                   % offset counter for positioning in global fracture vector
lios =0;
for i= 1:N_fractures
    Grf(lios+1:lios+Nf_i(i)) = gf(ios+2:ios+Nf_i(i)+1) - gf(ios+1:ios+Nf_i(i));
    ios = ios + Nf_i(i)+1; 
    lios = lios + Nf_i(i);
end  

%-------------------------------------------------------------------------%
%    Transmissivities                                                     %
%-------------------------------------------------------------------------%
Tx           = Txk./mux.*dx(2)./dx(1);                                     % Transmissivities
Ty           = Tyk./muy.*dx(1)./dx(2);
Tf           = Tfk./muf./dxf.*sqrt(12*Tfk);                                % The fracture aperture is calculated by the
                                                                           % cubic law instead of this fixed setting.
                                                                           