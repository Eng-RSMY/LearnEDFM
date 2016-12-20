function [Di,ri] = transport_heat_Diffusion(Nf,Nf_f,dx,phix,phiy)
%  Creates the diffusion matrix  
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
%  transportDiffusion(Nf,Nf_f,dx,phix,phiy)
%
%  Input: 
%        Nf     (2,1)           matrix grid dimensions
%        Nf_f                   total number of fracture segments
%        dx     (2,1)           matrix grid spacings
%        phix   (nx+1,ny)       matrix porosity at the interface in
%                               x-direction
%        phiy   (nx,ny+1)       matrix porosity at the interface in
%                               y-direction
%
%  Output: 
%        Di     (nx*ny+nf,nx*ny+nf) diffusion matrix of the transport system
%        ri     (nx*ny+nf,1)    diffusion rhs of the transport system

global lambda_l lambda_s ibcD FixT

Dx = dx(2).*(phix.*lambda_l+(1-phix).*lambda_s)./dx(1);
Dy = dx(1).*(phiy.*lambda_l+(1-phiy).*lambda_s)./dx(2);

Tdeast  = sparse(Nf(1),Nf(2));
Tdwest  = sparse(Nf(1),Nf(2));
Tdsouth = sparse(Nf(1),Nf(2));
Tdnorth = sparse(Nf(1),Nf(2));

%-------------------------------------------------------------%
%        preparation & setting of matrix Di                   %
%-------------------------------------------------------------%
Tdeast (2:Nf(1),:)   = Dx(2:Nf(1),:); Tdeast(1,:)      = 0;
Tdnorth(:,2:Nf(2))   = Dy(:,2:Nf(2)); Tdnorth(:,1)     = 0;
Tdwest (1:Nf(1)-1,:) = Dx(2:Nf(1),:); Tdwest(Nf(1),:)  = 0;
Tdsouth(:,1:Nf(2)-1) = Dy(:,2:Nf(2)); Tdsouth(:,Nf(2)) = 0;

Ds      = [Tdsouth(:) Tdwest(:) zeros(prod(Nf),1) Tdeast(:) Tdnorth(:)];
Ds(:,3) = -sum(Ds,2);
Di      = spdiags(-Ds,[-Nf(1),-1,0,1,Nf(1)],Nf(1)*Nf(2),Nf(1)*Nf(2));

%-------------------------------------------------------------%
%        boundary conditions (matrix only)                    %
%-------------------------------------------------------------%

ri   = sparse(prod(Nf),1);
rif  = zeros(Nf_f,1);

i1   = 1:Nf(2);                      i2  = Nf(2) + (1:Nf(2)); 
i3   = 2*Nf(2) + (1:Nf(1));          i4  = 2*Nf(2) + Nf(1) + (1:Nf(1));

ic1  = 1:Nf(1):prod(Nf);             ic2 = Nf(1):Nf(1):prod(Nf); 
ic3  = 1:Nf(1);                      ic4 = (Nf(2)-1)*Nf(1)+1:prod(Nf);

t1(1:Nf(2),1)  = Dx(1,1:Nf(2))'.*ibcD(i1);   
t2(1:Nf(2),1)  = Dx(Nf(1)+1,1:Nf(2))'.*ibcD(i2);
t3(1:Nf(1),1)  = Dy(1:Nf(1),1) .*ibcD(i3);   
t4(1:Nf(1),1)  = Dy(1:Nf(1),Nf(2)+1) .*ibcD(i4);   

iD  = [ic1';ic2';ic3';ic4'];
tD  = [t1;t2;t3;t4];
Di  = Di + sparse(iD,iD,tD,prod(Nf),prod(Nf));                     

ri(ic1,1) = ri(ic1,1) + t1.*FixT(i1,1);
ri(ic2,1) = ri(ic2,1) + t2.*FixT(i2,1);
ri(ic3,1) = ri(ic3,1) + t3.*FixT(i3,1);
ri(ic4,1) = ri(ic4,1) + t4.*FixT(i4,1);

%-------------------------------------------------------------%
%        merging matrix and fracture matrix and rhs           %
%-------------------------------------------------------------%
Di_f = spdiags(zeros(Nf_f,1),0,Nf_f,Nf_f);
Di = blkdiag(Di, Di_f);
ri = vertcat(ri,rif);
