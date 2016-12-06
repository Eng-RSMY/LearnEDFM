function [vfx,vfy] = calcFractureVelocity(XY1,vf)
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
%  calcFractureVelocity(XI,vf)
%
%  Input: 
%        XY1      (nx,ny)         matrix pressure
%        vf     (nf,1)          fracture pressure
%  Output:
%        vfx     (nf,1)       fracture velocity in the x direction
%        vfy     (nf,1)       fracture velocity in the y direction
 

%-------------------------------------------------------------------------%
%   Calculate velocity direction for alpha                                %
%-------------------------------------------------------------------------%
global Nf_i
xb = XY1(:,1);
yb = XY1(:,2);
xe = XY1(:,3);
ye = XY1(:,4);

alpha = atan2((ye-yb),(xe-xb));

ios = 0;                                                                   % offset counter for positioning in global fracture vector
lios = 0;
N_fractures = length(Nf_i);
vfx = zeros(length(Nf_i),1);
vfy = zeros(length(Nf_i),1);
for i= 1:N_fractures
    vfx(ios+1:ios+Nf_i(i))= vf(lios+1:lios+Nf_i(i)).*cos(alpha(ios+1:ios+Nf_i(i))); 
    vfy(ios+1:ios+Nf_i(i))= vf(lios+1:lios+Nf_i(i)).*sin(alpha(ios+1:ios+Nf_i(i))); 
    lios = lios + Nf_i(i) + 1;
    ios = ios + Nf_i(i);                                                   % Note the procedure on global indexing used here which is 
                                                                           % explained in the calculation of the fracture interface
                                                                           % permeabilities in the EDFM main file
end  