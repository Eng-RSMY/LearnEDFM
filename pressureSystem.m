function [u] = pressureSystem(pOld,pOldf,Fix,ibcs,Tx,Ty,Tf,g,gf,Q, Tfm,Tff, Nf_i,phi, phi_f, density_l,density_lf, K_f)
%  pressureSystem constructs the coefficient matrix and the right hand side
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
%  pressureSystem(Fix,ibcs,Tx,Ty,Tf,G,Q, Afm, Amf)
%
%  Input: 
%        Fix    (2*(nx+ny),1)   Fixed boundary values (Neumann or Dirichlet)
%                               (see note on ordering below)
%                               (here nr is the number of problems to solve)
%        ibcs   (2*(nx+ny),1)   Boundary condition type (1-Dirichlet, 0-Neumann)
%        Tx     (nx+1,ny)       trasmissibility in the x direction
%        Ty     (nx,ny+1)       trasmissibility in the y direction
%        Tf     (nf,1)          fracture transmissivity
%        G      (nx,ny)         gravity term appearing on the right hand side
%        Gf     (nf,1)          fracture gravity term (RHS)
%        Q      (nx,ny)         right hand side of the linear problem
%        Tfm    (nf,nx*ny)      coupling matrix between fracture and matrix
%        Tmf    (nx*ny,nf)      coupling matrix between matrix and fracture
%        Tff    (nf*nf,nf*nf)   fracture-fracture intersection
%                               transmissivity matrix
%        Ni_i   (X,1)           fracture length information array                
%
%                    Ordering for cells
%                        (6 x 7 grid)
%                    -------------------
%                    |37 38 39 40 41 42|
%                    |31 32 33 34 35 36|
%                    |25 26 27 28 29 30|
%                 Y  |19 20 21 22 23 24|
%                    |13 14 15 16 17 18|                
%                    | 7  8  9 10 11 12|
%                    | 1  2  3  4  5  6|
%                    -------------------    
%                             X
%
%                     21 22 23 24 25 25       
%                    -------------------    
%                   7|                 |14     
%                   6|     ordering    |13    
%                   5|        for      |12                    
%                y  4|     boundary    |11    
%                   3|      faces      |10     
%                   2|   (6 x 7 grid)  |9               
%                   1|                 |8         
%                    -------------------    
%                     15 16 17 18 19 20
%                             x
%  Output: A     (nx*ny+nf,nx*ny+nf)   coefficient matrix
%          rhs   (nx*ny+nf)         right hand side of the linear system
%
%
%  Acknowledgement:  thanks are due to Brad Mallison (Chevron ETC) for
%                    providing the initial core of this function
%
%-------------------------------------------------------------------------%
n    = size(Tx) - [1 0];                                                   % Calculate the logical grid dimension (Tx has dimension (nx+1,ny)...)
nf   = length(Tf) -length(Nf_i);                                           % Calculate the logical grid dimension (Tf has dimension (nf+1,1)...)
ibcs = logical(ibcs);                                                      % Convert to logical

N_fractures = length(Nf_i);                                                % Calculate the number of fractures

%-------------------------------------------------------------------------%
%    creating the transmissibility matrix                                 %
%-------------------------------------------------------------------------%

Txeast  = sparse(n(1),n(2));
Txwest  = sparse(n(1),n(2));
Tysouth = sparse(n(1),n(2));
Tynorth = sparse(n(1),n(2));

Txeast (2:n(1),:)   = Tx(2:n(1),:); Txeast(1,:)     = 0;
Tynorth(:,2:n(2))   = Ty(:,2:n(2)); Tynorth(:,1)    = 0;
Txwest (1:n(1)-1,:) = Tx(2:n(1),:); Txwest(n(1),:)  = 0;
Tysouth(:,1:n(2)-1) = Ty(:,2:n(2)); Tysouth(:,n(2)) = 0;

[frac_left_flux_mat, frac_right_flux_mat] = calc_frac_flux_mat(N_fractures,nf, Nf_i);

%-------------------------------------------------------------------------%
%    preparing and setting the diagonals of the matrix A                  %
%-------------------------------------------------------------------------%
Ds      = [Tysouth(:) Txwest(:) zeros(prod(n),1) Txeast(:) Tynorth(:)];
Ds(:,3) = -sum(Ds,2);
A       = spdiags(-Ds,[-n(1),-1,0,1,n(1)],n(1)*n(2),n(1)*n(2));

Dsf      = [frac_right_flux_mat*Tf zeros(nf,1) frac_left_flux_mat*Tf];
Dsf(:,2) = -sum(Dsf,2);
Af       = spdiags(-Dsf,[-1,0,1],nf,nf);

%-------------------------------------------------------------------------%
%    preparing the right hand side of the equation                        %
%-------------------------------------------------------------------------%

rhs = sparse(reshape(Q + (g(:,2:n(2)+1) - g(:,1:n(2))),prod(n),1));

global frac_grad_mat
rhsf = frac_grad_mat*gf;

%-------------------------------------------------------------------------%
%    boundary conditions                                                  %
%-------------------------------------------------------------------------%

i1   = 1:n(2);                      i2  = n(2) + (1:n(2));                 % Ordering of bc-vector [1 : 2nf(2)+2nf(1)]
i3   = 2*n(2) + (1:n(1));           i4  = 2*n(2) + n(1) + (1:n(1));

ic1  = 1:n(1):prod(n);              ic2 = n(1):n(1):prod(n);               % Ordering of Stiffness Matrix 
ic3  = 1:n(1);                      ic4 = (n(2)-1)*n(1)+1:prod(n);

%-----------%
% Dirichlet %
%-----------%

iD1  = ic1(ibcs(i1))';              iD2 = ic2(ibcs(i2))'; 
iD3  = ic3(ibcs(i3))';              iD4 = ic4(ibcs(i4))';

t1   = Tx(     1,ibcs(i1))';   
t2   = Tx(n(1)+1,ibcs(i2))';
t3   = Ty(ibcs(i3),     1);   
t4   = Ty(ibcs(i4),n(2)+1);
 
iD  = [iD1;iD2;iD3;iD4];
tD  = [t1;t2;t3;t4];
A   = A + sparse(iD,iD,tD,prod(n),prod(n));

rhs(iD1,:) = rhs(iD1,:) + t1.*Fix(i1(ibcs(i1)),:);
rhs(iD2,:) = rhs(iD2,:) + t2.*Fix(i2(ibcs(i2)),:);
rhs(iD3,:) = rhs(iD3,:) + t3.*Fix(i3(ibcs(i3)),:);
rhs(iD4,:) = rhs(iD4,:) + t4.*Fix(i4(ibcs(i4)),:);

%-----------%
%  Neumann  %
%-----------%

iN1  = ic1(~ibcs(i1))';   iN2 = ic2(~ibcs(i2))'; 
iN3  = ic3(~ibcs(i3))';   iN4 = ic4(~ibcs(i4))';

rhs(iN1,:) = rhs(iN1,:) + Fix(i1(~ibcs(i1)),:);
rhs(iN2,:) = rhs(iN2,:) + Fix(i2(~ibcs(i2)),:);
rhs(iN3,:) = rhs(iN3,:) + Fix(i3(~ibcs(i3)),:);
rhs(iN4,:) = rhs(iN4,:) + Fix(i4(~ibcs(i4)),:);

%-------------------------------------------------------------------------%
%    merging matrix and fracture matrix and rhs                           %
%-------------------------------------------------------------------------%
Tmf = Tfm';
Ds = sum(Tmf,2);
DsT = sum(Tfm,2);
[m,n]=size(A);
B = spdiags(Ds,0,m,n);                                                     % Diagonal contribution of Amf to the main diagonal of A 
[m,n]=size(Af);
Bf = spdiags(DsT,0,m,n);                                                   % Diagonal contribution of Afm to the main diagonal of A 
Ds = -sum(Tff,2);                                                          % Note: The - sign is needed here!
Bff = spdiags(Ds,0,nf,nf);                                                 % Diagonal contribution of Aff to the main diagonal of Af 
Tff = Tff +Bff;

A = A + B;
Af = Af + Bf - Tff;

A = [A -Tmf; Tfm -Af];                                                     % Build the coupled pressure block matrix
rhs = vertcat(rhs,rhsf);  

u = A\rhs;                                                         % Solve the pressure equation