function [vx,vy,vf,Vfm,Vmf,Vff] = calcVelocity(p,pf,Tx,Ty,Tf,Tfm,Tmf,Tff,g,gf)
%  calcVelocity Calculates the velocity field from the pressure
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
%  Acknowledgement:  thanks are due to Manav Tyagi and Hadi Hajibeygi for
%                    contributing to the very early development of the code. 
%
%  calcVelocity(p,pf,Tx,Ty,Tf,Tfm,Tmf,Tff,g,gf)
%
%  Input: 
%        p      (nx,ny)         matrix pressure
%        pf     (nf,1)          fracture pressure
%        Tx     (nx+1,ny)       trasmissibility in the x direction
%        Ty     (nx,ny+1)       trasmissibility in the y direction
%        Tf     (nf,1)          fracture transmissivity
%        Tfm    (nf,nx*ny)      coupling matrix between fracture and matrix
%        Tmf    (nx*ny,nf)      coupling matrix between matrix and fracture
%        Tff    (nf*nf,nf*nf)   fracture-fracture intersection
%                               transmissivity
%        g      (nx,ny)         gravity term with respect to the interfaces
%        gf     (nf,1)          gravity term with respect to the fracture interfaces
%  Output:
%        vx     (nx+1,ny)       matrix velocity in the x direction
%        vy     (nx,ny+1)       matrix velocity in the y direction
%        vf     (nf,1)          fracture velocity
%        Vfm    (nf,nx*ny)      velocity matrix between fracture and matrix
%        Vmf    (nx*ny,nf)      velocity matrix between matrix and fracture
%        Vff    (nf*nf,nf*nf)   fracture-fracture intersection velocity
 
global Fix ibcs Nf_i

N = size(p);
Nf = size(pf);

%-------------------------------------------------------------------------%
%   Calculate velocity direction for alpha                                %
%-------------------------------------------------------------------------%

%-------------- velocity inside matrix cells -----------------------------%
vx = sparse(N(1)+1,N(2));
vy = sparse(N(1),N(2)+1);

vx(2:N(1),:) = -Tx(2:N(1),:).*(p(2:N(1),:)-p(1:N(1)-1,:));                 % velocity direction in x direction
vy(:,2:N(2)) = -Ty(:,2:N(2)).*(p(:,2:N(2))-p(:,1:N(2)-1)) - g(:,2:N(2));   % velocity direction in y direction

%--------------------- velocity on matrix boundaries ---------------------%                                 
  
vx(1,:) = -Tx(1,:).*(p(1,:)-Fix(1:N(2))').*ibcs(1:N(2))';                                                        % West Dirichlet
vx(1,:) = vx(1,:)+(Fix(1:N(2)).*~ibcs(1:N(2)))';                                                                 % West Neumann

vx(N(1)+1,:) = -Tx(N(1)+1,:).*(Fix(N(2)+1:2*N(2))'-p(N(1),:)).*ibcs(N(2)+1:2*N(2))';                             % East (adjust Neumann sign!)
vx(N(1)+1,:) = vx(N(1)+1,:)-(Fix(N(2)+1:2*N(2)).*~ibcs(N(2)+1:2*N(2)))';

vy(:,1) = (-Ty(:,1).*(p(:,1)-Fix(2*N(2)+1:2*N(2)+N(1))) - g(:,1)).*ibcs(2*N(2)+1:2*N(2)+N(1));                   % South
vy(:,1) = vy(:,1)+(Fix(2*N(2)+1:2*N(2)+N(1)).*~ibcs(2*N(2)+1:2*N(2)+N(1)));   

vy(:,N(2)+1) = (-Ty(:,N(2)+1).*(Fix(2*N(2)+N(1)+1:2*N(2)+2*N(1))-p(:,N(2))) - g(:,N(2)+1)).*ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1));   % North (adjust Neumann sign!)
vy(:,N(2)+1) = vy(:,N(2)+1)-(Fix(2*N(2)+N(1)+1:2*N(2)+2*N(1)).*~ibcs(2*N(2)+N(1)+1:2*N(2)+2*N(1)));

%--------------------- velocity in fractures -----------------------------%                                 
vf = zeros(Nf(1)+length(Nf_i),Nf(2));
ios = 0;                                                                   % offset counter for positioning in global fracture vector
lios = 0;
N_fractures = length(Nf_i);
for i= 1:N_fractures
    vf(lios+2:lios+Nf_i(i)) = -Tf(ios+2:ios+Nf_i(i)).*(pf(ios+2:ios+Nf_i(i))-pf(ios+1:ios+Nf_i(i)-1) - gf(ios+2:ios+Nf_i(i))); 
    lios = lios + Nf_i(i) + 1;
    ios = ios + Nf_i(i);                                                   % Note the procedure on global indexing used here which is 
                                                                           % explained in the calculation of the fracture interface
                                                                           % permeabilities in the EDFM main file
end  

%--------------------- velocity fracture-fracture-------------------------%                                 
Vff = bsxfun(@times, pf, Tff) - bsxfun(@times,pf',Tff);

%--------------------- velocity fracture-matrix---------------------------% 
Vfm = bsxfun(@times, pf, Tfm) - bsxfun(@times, p(:)',Tfm);
%--------------------- velocity matrix-fracture---------------------------%
Vfm = -Vfm;
Vmf = Vfm';