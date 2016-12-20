%  Generates two perpendicular intersecting fractures
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

N_fractures = 2;

dx_f = sqrt(dx(1)*dx(1)+dx(2)*dx(2));

xc = [];
yc = [];
xb = [];
yb = [];
xe = [];
ye = [];

Nf_i = [];

dxi = 0;
dyi = dx_f;

x = 4.5; y=2;
r = floor(5/9.*min(len(1),len(2))./dx_f);
ybi = y:dyi:y+r*dyi-dyi; 
yei = y+dyi:dyi:y+r*dyi;


%ybi = 2:dy:7-dy;
%yei = 2+dy:dy:7;
xci = 4.5*ones(size(yei));

xb =[xb round(xci,4)];
xe =[xe round(xci,4)];
yb =[yb round(ybi,4)];
ye =[ye round(yei,4)];

Nf_i(1) = length(xci);

dxi = dx_f;
dyi = 0;

x = 2; y=4.5;
xbi = x:dxi:x+r*dxi-dxi; 
xei = x+dxi:dxi:x+r*dxi;

%xbi = 2:dxl:7-dxl; 
%xei = 2+dxl:dxl:7;
%xci = 2+dxl/2:dxl:7;
yci = 4.5*ones(size(xci));

Nf_i(2) = length(xei);

xc =[xc round(xci,4)];
yc =[yc round(yci,4)];
xb =[xb round(xbi,4)];
xe =[xe round(xei,4)];
yb =[yb round(yci,4)];
ye =[ye round(yci,4)];

ListOfVariables = who;
for k = 1:length(ListOfVariables)
   assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
end 

XY1 = [xb' yb' xe' ye'];
XY2(:,1,1) = xb';
XY2(:,1,2) = yb';
XY2(:,2,1) = xe';
XY2(:,2,2) = ye';
XE = [xe; ye]';

Nf_f = length(xe);
frac_angle = atan2((ye-yb),(xe-xb))*180/pi;
dxf = sqrt((xb(1)-xe(1))^2 + (yb(1)-ye(1))^2);