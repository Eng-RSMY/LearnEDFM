function [Tx,Ty,Tf,Tfm,Tff,g,gf,density_l,density_lf] = initialize(K, K_f,c,cf,T,Tf,p,pf,CI,row,col)
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
%        c      (nx,ny)         matrix tracer concentration
%        cf     (nf,1)          fracture tracer concentration
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

global dx dxf gravity Nf Nf_f Nf_i frac_angle
global frac_mean_mat frac_grad_mat 
global FixT ibcs  

N = size(c);

%-------------------------------------------------------------------------%
%    Initialize                                                           %
%-------------------------------------------------------------------------%
viscosity = 2.414e-5*10.^(247.8./(T+273.15-140));
viscosity_f = 2.414e-5*10.^(247.8./(Tf+273.15-140));

[density_l, density_lf] = calc_density(T,Tf,p,pf);

%-------------------------------------------------------------------------%
% Calculate matrix interface permeability & porosity (harmonic average)
%-------------------------------------------------------------------------%
[Kx, Ky] = calc_interface_values(K);
[Kf] = calc_interface_values_fracture(K_f);

[viscosityx, viscosityy] = calc_interface_values(viscosity);
[viscosityf] = calc_interface_values_fracture(viscosity_f);

[~, densityy] = calc_interface_values(density_l);
[densityf] = calc_interface_values_fracture(density_lf);


%-------------------------------------------------------------------------%
%    Gravity for the matrix                                               %
%-------------------------------------------------------------------------%

Gy    = densityy.*gravity;                    % TODO: Add the temperature dependence here
g            = Ky./viscosityy.*Gy.*dx(1);                                        % Gravity term with respect to the interfaces
g(:,1)       = g(:,1).*ibcs(2*N(2)+1:2*N(2)+N(1)); 
g(:,N(2)+1)  = g(:,N(2)+1).*ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1)); 

%-------------------------------------------------------------------------%
%    Gravity for the fractures                                            %
%-------------------------------------------------------------------------%
interface_frac_angle = 0.5*frac_mean_mat*frac_angle';
Gyf= densityf.*gravity.*sin(interface_frac_angle*pi/180);%.*dx_f;         

gf           = -Kf./viscosityf.*Gyf.*sqrt(12*Kf);                                % Cubic law is used here instead of fixed aperture
Ni=zeros(size(Nf_i));
N_fractures = length(Nf_i);
for i=1:N_fractures
    for j=i:-1:1
        Ni(i) = Ni(i) +Nf_i(j);
    end
end
Nii = 1:1:length(Ni);

gf(Ni-Nf_i+Nii) = 0; 
gf(Ni+Nii) = 0;

%-------------------------------------------------------------------------%
%    Transmissivities                                                     %
%-------------------------------------------------------------------------%
Tx           = Kx./viscosityx.*dx(2)./dx(1);                                     % Transmissivities
Ty           = Ky./viscosityy.*dx(1)./dx(2);
Tf           = Kf./viscosityf./dxf.*sqrt(12*Kf);                                % The fracture aperture is calculated by the
                                                                           % cubic law instead of this fixed setting.
                                                                                                                                             
%% Matrix-Fracture transmissivity
%  This assembles the matrix-fracture transmissivity based on the
%  previously computed connectivity index
l = Nf_f;
n = Nf;
X = zeros(length(CI),1);
for i = 1:length(CI)
    indm = CI(i,1);
    indf = CI(i,2);
    lambda_ij = K(indm)/viscosity(indm);
    lambda_k = K_f(indf)/viscosity_f(indf);
    X(i)  =  CI(i,3)*2*lambda_ij*lambda_k/(lambda_ij+lambda_k);                  %Harmonic mean for the fracture-matrix transmissivity 
    [ii,jj] = ind2sub(Nf,CI(i,1));          
    I(i) = (jj-1)*n(1)+ii;
end
if (isempty(CI))
    Tfm = sparse(zeros(n(1)*n(2),l));
else
    Tfm = sparse(I,CI(:,2),X,n(1)*n(2),l);
end
Tfm = Tfm';


% Find intersection points between the fracture segments
alpha_row = sqrt(12*K_f(row)).*K_f(row)./viscosityf(row)./(0.5*dxf);                   % Cubic law is used here instead of fixed aperture
alpha_col = sqrt(12*K_f(col)).*K_f(col)./viscosityf(col)./(0.5*dxf);                   % Cubic law is used here instead of fixed aperture
Tflk = alpha_row.*alpha_col./(alpha_row+alpha_col);
Tff = sparse(row,col,Tflk,l,l); 