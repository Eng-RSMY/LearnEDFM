function [tNew] = transport_heat_System(tN,tNf,vx,vy,vf,Vfm,Vmf,Vff,Q,QT,K_f,phi,phi_f,density_s,density_sf,density_l,density_lf)
%  Assembles and solves the transport system. It consists of the
%  subfunctions for advective and diffusive transport
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
%  Acknowledgement:  Thanks are due to Manav Tyagi and Hadi Hajibeygi for
%                    contributing to the very early development of the code.
%
%  transportSystem(sN,sNf,vx,vy,vf,Vfm,Vmf,Vff,Dfm,Dmf,phix,phiy)
%
%  Input: 
%        tN     (nx,ny)         old heat transport solution on the matrix
%        tNf    (nf,1)          old heat transport solution in the fractures
%        vx     (nx+1,ny)       matrix velocity in the x direction
%        vy     (nx,ny+1)       matrix velocity in the y direction
%        vf     (nf,1)          fracture velocity
%        Vfm    (nf,nx*ny)      velocity matrix between fracture and matrix
%        Vmf    (nx*ny,nf)      velocity matrix between matrix and fracture
%        Vff    (nf*nf,nf*nf)   fracture-fracture intersection velocity
%        phix   (nx+1,ny)       matrix porosity at the interface in
%                               x-direction
%        phiy   (nx,ny+1)       matrix porosity at the interface in
%                               y-direction
%        cprho  (2,1)           constants for the heat transport:
%                               cprho(1) -> fluid / cprho(2) -> rock
%
%  Output: 
%        tNew  (nx*ny+nf,1)     new heat transport solution vector


global lambda_l lambda_s cp_l cp_s dx dxf dt Nf Nf_f Nf_i 

%-------------------------------------------------------------------------%
% Calculate matrix interface permeability & porosity (harmonic average)
%-------------------------------------------------------------------------%
[phix, phiy] = calc_interface_values(phi);

%-------------------------------------------------------------------------%
%    Build Transport System                                               %
%-------------------------------------------------------------------------%
[Up,r]  = transport_heat_Advection(vx,vy,vf,Vfm,Vmf,Vff,Q,QT,Nf,Nf_f,Nf_i,tN,tNf,cp_l.*density_l,cp_l.*density_lf);        % Convective Upwind Matrix

if (phi*lambda_l+(1-phi)*lambda_s) == 0;
   Di = sparse(prod(Nf)+Nf_f,prod(Nf)+Nf_f);
   ri = sparse(prod(Nf)+Nf_f,1);
else
  [Di,ri] = transport_heat_Diffusion(Nf,Nf_f,dx,phix,phiy);              % Diffusive Matrix 
end

Ac      = sparse((phi.*cp_l.*density_l+(1-phi).*cp_s.*density_s).*prod(dx)/dt);                                        % Accumulation term matrix
Acf = sparse((phi_f.*cp_l.*density_lf+(1-phi_f).*cp_s.*density_sf).*dxf.*sqrt(12.*K_f)/dt);                                % Accumulation term fracture

Ac = [Ac(:); Acf(:)];                                                      % Combine the accumulation terms
snnf = [tN(:); tNf(:)];                                                    % for matrix and fracture

A       = Up + diag(Ac(:)) + Di;                                           % Stiffness matrix
rhs     = r + ri + sparse(snnf(:).*Ac(:));                                 % Right hand side

%-------------------------------------------------------------------------%
%    Solve Transport System                                               %
%-------------------------------------------------------------------------%
tNew = A\rhs;
