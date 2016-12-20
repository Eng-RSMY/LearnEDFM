%  Generates a single fracture (either horizontal or vertical)
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

global dxf Nf_i frac_angle XY1

N_fractures = 1;

vertical = 1;

dx_f = 1.2*sqrt(dx(1)*dx(1)+dx(2)*dx(2));

xc = [];
yc = [];
xb = [];
yb = [];
xe = [];
ye = [];

Nf_i = [];


if(vertical)
    dxl = 0;
    dy = dx_f;

    ybi = 2:dy:7-dy;
    yei = 2+dy:dy:7;
    yci = 2+dy/2:dy:7;
    xci = 4.5*ones(size(yci));
    xc =[xc round(xci,4)];
    yc =[yc round(yci,4)];
    xb =[xb round(xci,4)];
    xe =[xe round(xci,4)];
    yb =[yb round(ybi,4)];
    ye =[ye round(yei,4)];

    Nf_i(1) = length(xci);
    frac_angle = [0 0];
else
    dxl = dx_f;
    dy = 0;
    xbi = 2:dxl:7-dxl; 
    xei = 2+dxl:dxl:7;
    xci = 2+dxl/2:dxl:7;
    yci = 4.5*ones(size(xci));

    Nf_i(1) = length(xei);

    xc =[xc round(xci,4)];
    yc =[yc round(yci,4)];
    xb =[xb round(xbi,4)];
    xe =[xe round(xei,4)];
    yb =[yb round(yci,4)];
    ye =[ye round(yci,4)];
    frac_angle = [90 0];
end

XY1 = [xb' yb' xe' ye'];
XE = [xe; ye]';

Nf_f = length(xe);
frac_angle = atan2((ye-yb),(xe-xb))*180/pi;
dxf = sqrt((xb(1)-xe(1))^2 + (yb(1)-ye(1))^2);