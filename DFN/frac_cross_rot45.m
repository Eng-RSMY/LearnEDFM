%  Generates two perpendicular intersecting fractures
%  Rotated 45° from the coordinate axis
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

global dx dxf Nf_i frac_angle len

N_fractures = 2;
Nf_i = [];
dxf = sqrt(dx(1)*dx(1)+dx(2)*dx(2));

xc = []; yc = [];
xb = []; yb = [];
xe = []; ye = [];

beta = [-45 45]; % Fracture angle [°]
xy = [0.225 0.725; 0.225 0.225];

for i=1:N_fractures 
    pass = 0;
    while pass < 1
        
        x = xy(i,1)*len(1);
        y = xy(i,2)*len(2);
        r = 18;
        
        dyi = sin(beta(i)*pi/180)*dxf;
        dxi = sqrt(dxf*dxf-dyi*dyi);

        xbi = x:dxi:x+r*dxi-dxi; 
        xei = x+dxi:dxi:x+r*dxi;
        ybi = y:dyi:y+r*dyi-dyi; 
        yei = y+dyi:dyi:y+r*dyi;

        if (all(xbi <= len(1)) && all(xei <= len(1)) && ...
            all(ybi <= len(2)) && all(yei <= len(2)) && ...
            all(xbi >= 0) && all(xei >= 0) && ...
            all(ybi >= 0) && all(yei >= 0))
            pass = 1;
        end
    end
    
    Nf_i(i) = length(xei);

    xb =[xb xbi];
    xe =[xe xei];
    yb =[yb ybi];
    ye =[ye yei];
end

Nf_f = length(xe);
frac_angle = beta+90;
XY1 = [xb' yb' xe' ye'];
XE = [xe; ye]';

% Comment in the following code to visualize the fracture network before
% the simulation is started
% -------------------------------------------------------------------------
%figure
%line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r');
%pause()
