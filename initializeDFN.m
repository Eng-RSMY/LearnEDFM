function [Tfm,Tmf,Tff,Dfm,Dmf,Z] = initializeDFN(x,y,XY1,K_f,lambda,lambda_f)
%InitializeDFN updates the connectivity and transmissivity of the DFN
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
%
%  InitializeDFN(x,y,XY1,K,K_f,visc)
%
%  Input: 
%        x      (1, nx)         x-direction grid cell centers
%        y      (1, ny)         y-direction grid cell centers
%        XY1    (nf, 4)         array containing all fracture segments as
%                               [x_begin y_begin x_end y_end]
%        K_f(nf,1)              fracture permeability
%        lambda (nx,ny)         matrix permeability/viscosity
%        lambda_f(nf,1)         fracture permeability/viscosity
%
%  Output:
%        Tfm    (nf,nx*ny)      coupling matrix between fracture and matrix
%        Tmf    (nx*ny,nf)      coupling matrix between matrix and fracture
%        Tff    (nf*nf,nf*nf)   fracture-fracture intersection
%                               transmissivity
%        Dfm    (nf,nx*ny)      coupling matrix between fracture and matrix
%                               for transport diffusion
%        Dmf    (nx*ny,nf)      coupling matrix between matrix and fracture
%                               for transport diffusion
%        Z      (nx*ny,nx*ny)   fracture-grid intersections (boolean array)

global dxf Nf Dif phi phi_f

%% Find the intersection of each segment with the grid and compute the
%% Connectivity Index based on the segments lengths and the average distance
%% to the fracture
[Z,CI] = intersectionsGrid(x,y,XY1);

%% Matrix-Fracture transmissivity
%  This assembles the matrix-fracture transmissivity based on the
%  previously computed connectivity index
l = length(XY1);

I = zeros(Nf(1)*Nf(2)*l,1);
J = zeros(Nf(1)*Nf(2)*l,1);
X = zeros(Nf(1)*Nf(2)*l,1);
X2 = zeros(Nf(1)*Nf(2)*l,1);
nind = 0;

for j=1:Nf(2)
    for i=1:Nf(1)
        for k=1:l
        nind = nind +1;    
            I(nind) = (j-1)*Nf(1)+i;
            J(nind) = k;
            X(nind) = CI(i,j,k)*2*lambda(i,j)*lambda_f(k)/(lambda(i,j)+lambda_f(k));
            X2(nind) = CI(i,j,k)*Dif*2*phi(i,j)*phi_f(k)/(phi(i,j)+phi_f(k));
        end
    end
end
Tfm = sparse(I,J,X,Nf(1)*Nf(2),l);

Tmf = Tfm; %./(dx(1)*dx(2)); % This additional normalization is in the literature from Norbeck and Hajibeygi but 
Tfm = Tfm';%./(dxf);         % does not show up in the discretized equations of Sanders. I think without the scaling it should be correct.

Dfm = sparse(I,J,X2,Nf(1)*Nf(2),l);
Dmf = Dfm;
Dfm = Dfm';

% Find intersection points between the fracture segments
[row,col] = intersectionsSegments(XY1,XY1);
alpha_row = sqrt(12*K_f(row)).*lambda_f(row)./(0.5*dxf);                   % Cubic law is used here instead of fixed aperture
alpha_col = sqrt(12*K_f(col)).*lambda_f(col)./(0.5*dxf);                   % Cubic law is used here instead of fixed aperture
Tflk = alpha_row.*alpha_col./(alpha_row+alpha_col);
Tff = sparse(row,col,Tflk,l,l);
