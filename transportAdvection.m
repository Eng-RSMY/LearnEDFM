function [Up,rhs] = transportAdvection(vx,vy,vf,Vfm,Vmf,Vff,Nf,Nf_f,Nf_i)
%UPMAT
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
%  Acknowledgement:  Thanks are due to Hadi Hajibeygi for contributing to 
%                    the very early development of the code. 
%
%  transportAdvection(vx,vy,vf,Vfm,Vmf,Vff,Nf,Nf_f,Nf_i)
%
%  Input: 
%        vx     (nx+1,ny)       matrix velocity in the x direction
%        vy     (nx,ny+1)       matrix velocity in the y direction
%        vf     (nf,1)          fracture velocity
%        Vfm    (nf,nx*ny)      velocity matrix between fracture and matrix
%        Vmf    (nx*ny,nf)      velocity matrix between matrix and fracture
%        Vff    (nf*nf,nf*nf)   fracture-fracture intersection velocity
%        Nf     (2,1)           matrix grid dimensions
%        Nf_f                   total number of fracture segments
%        Nf_i                   indivdidual numer of fracture segments
%
%  Output: 
%        Up    (nx*ny+nf,nx*ny+nf) advection matrix of the transport system
%        rhs   (nx*ny+nf,1)      advection rhs of the transport system

global FixT Q QT

n    = Nf;
nf   = Nf_f;
N_fractures = length(Nf_i);                                                % Calculate the number of fractures


%-------------------------------------------------------------------------%
%    Matrix                                                               %
%-------------------------------------------------------------------------%
ix      = (1+sign(vx))/2;                                             % Velocity indicator for upwinding
iy      = (1+sign(vy))/2;                                             % for aritmetic mean set ix = iy = 1/2
     
upw(1:Nf(1),:) = -ix(1:Nf(1),:)   .*vx(1:Nf(1),:);
ups(:,1:Nf(2)) = -iy(:,1:Nf(2))   .*vy(:,1:Nf(2));
upe(1:Nf(1),:) =  (1-ix(2:Nf(1)+1,:)) .*vx(2:Nf(1)+1,:);
upn(:,1:Nf(2)) =  (1-iy(:,2:Nf(2)+1)) .*vy(:,2:Nf(2)+1);

Txeast(2:n(1),:)   = upe(1:n(1)-1,:);    Txeast(1,:)     = 0;
Tynorth(:,2:n(2))  = upn(:,1:n(2)-1);    Tynorth(:,1)    = 0;
Txwest(1:n(1)-1,:) = upw(2:n(1),:);      Txwest(n(1),:)  = 0;
Tysouth(:,1:n(2)-1)= ups(:,2:n(2));      Tysouth(:,n(2)) = 0;

upd(1:Nf(1),1:Nf(2)) = ix(2:Nf(1)+1,1:Nf(2))   .*vx(2:Nf(1)+1,1:Nf(2))...
                      +iy(1:Nf(1),2:Nf(2)+1)   .*vy(1:Nf(1),2:Nf(2)+1)...
                      -(1-ix(1:Nf(1),1:Nf(2))) .*vx(1:Nf(1),1:Nf(2))...
                      -(1-iy(1:Nf(1),1:Nf(2))) .*vy(1:Nf(1),1:Nf(2));
Ds      = [Tysouth(:) Txwest(:) upd(:) Txeast(:) Tynorth(:)];

Up   = spdiags(Ds,[-n(1) -1 0 1 n(1)],n(1)*n(2),n(1)*n(2));

rhs = sparse(zeros(n(1)*n(2),1));

%------------------------------------------------------%
%                   boundary conditions                %
%------------------------------------------------------%

i1   = 1:n(2);                      i2  = n(2) + (1:n(2)); 
i3   = 2*n(2) + (1:n(1));           i4  = 2*n(2) + n(1) + (1:n(1));
ic1  = 1:n(1):prod(n);              ic2 = n(1):n(1):prod(n); 
ic3  = 1:n(1);                      ic4 = (n(2)-1)*n(1)+1:prod(n);

%---------------------%
% Dirichlet Transport %
%---------------------%

iD1  = ic1';              iD2 = ic2'; 
iD3  = ic3';              iD4 = ic4';
t1         = upw(1,:)';   
t2         = upe(n(1),:)';
t3         = ups(:,1);   
t4         = upn(:,n(2));

rhs(iD1,1) = rhs(iD1,1) - t1.*FixT(i1,1);
rhs(iD2,1) = rhs(iD2,1) - t2.*FixT(i2,1);
rhs(iD3,1) = rhs(iD3,1) - t3.*FixT(i3,1);
rhs(iD4,1) = rhs(iD4,1) - t4.*FixT(i4,1);

%------------------------------------------------------%
%    source terms                                      %
%------------------------------------------------------%                                                    
iQ  = sparse((1+sign(Q(:)))/2);                                              % Indicator for in or outflow

out = spdiags((1-iQ).*Q(:),0,prod(n),prod(n));
in  = iQ.*Q(:).*QT(:);

rhs = sparse(rhs + in);
Up  = sparse(Up  - out);

%-------------------------------------------------------------------------%
%    Fracture                                                             %
%-------------------------------------------------------------------------%
Tfeast = zeros(nf,1);
Tfwest = zeros(nf,1);
updf   = zeros(nf,1);

upwf = zeros(nf,1);
upef = zeros(nf,1);

wf = min(-vf,0); 
pf = max(-vf,0);

lios=0;
ios = 0;                                                                   % offset counter for positioning in global fracture vector
for i= 1:N_fractures
    upef(ios+1:ios+Nf_i(i)) = -pf(lios+2:lios+Nf_i(i)+1);
    upwf(ios+1:ios+Nf_i(i)) = +wf(lios+1:lios+Nf_i(i));
    
    updf(ios+1:ios+Nf_i(i)) = pf(ios+2:ios+Nf_i(i)+1)...
                      -wf(ios+1:ios+Nf_i(i));
        
    Tfeast(ios+2:ios+Nf_i(i))     = upef(ios+1:ios+Nf_i(i)-1);
    Tfwest(ios+1:ios+Nf_i(i)-1)   = upwf(ios+2:ios+Nf_i(i));
    Tfeast(ios+1)     = 0;
    Tfwest(ios+Nf_i(i)) =0;                                                % Note the procedure on global indexing used here which is 
                                                                           % explained in the calculation of the fracture interface
                                                                           % permeabilities in the EDFM main file
    ios = ios + Nf_i(i);
    lios = lios + Nf_i(i) +0 ;
end 
       
Dsf  = [Tfwest updf Tfeast];        
Upf  = spdiags(Dsf,[-1 0 1],nf,nf);
rhsf = sparse(zeros(nf,1));

%-------------------------------------------------------------------------%
%    Matrix-Fracture                                                      %
%-------------------------------------------------------------------------%
Upmf = -max(-Vmf,0);
Ds =  sum(max(Vmf,0),2);
[m,n]=size(Up);
B = spdiags(Ds,0,m,n);                                                     % Diagonal contribution of Amf to the main diagonal of A 
Up = Up + B;

%-------------------------------------------------------------------------%
%    Fracture-Matrix                                                      %
%-------------------------------------------------------------------------%
Upfm = -max(Vfm,0);
DsT = sum(max(Vfm,0),2);
[m,n]=size(Upf);

%-------------------------------------------------------------------------%
%    Fracture-Fracture                                                    %
%-------------------------------------------------------------------------%
DsTT = sum(max(Vff,0),2);
Tff = -max(Vff,0);
DsT = DsT + DsTT;
Bf = spdiags(DsT,0,m,n);                                                   % Diagonal contribution of Afm and Aff to the main diagonal of A 
Upf = Upf + Bf + Tff;

%-------------------------------------------------------------%
%        merging matrix and fracture matrix and rhs           %
%-------------------------------------------------------------%
Up = [Up Upmf; Upfm Upf]; 
rhs = vertcat(rhs,rhsf);                                                   % Concenate the RHS vectors of matrix and fracture
