%  Loads a discrete fracture network from a fracsim csv file
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

[x1,y1,x2,y2] = read_fracman_dfn_file(filename);

xshift = min(min(x1),min(x2));
yshift = min(min(y1),min(y2));

x1 = x1 + abs(xshift);
x2 = x2 + abs(xshift);
y1 = y1 + abs(yshift);
y2 = y2 + abs(yshift);

maxx = max(max(x1),max(x2));
maxy = max(max(y1),max(y2));

len     = [maxx maxy];                    % physical length of the domain in x and y direction [m]
dx      = len./Nf;                      % cell length [m]


dfn = [x1 y1 x2 y2];

dfn = round(dfn,3);

N_fractures = length(dfn);

for i=1:N_fractures
    if dfn(i,1) > dfn(i,3) 
        if dfn(i,2) ~= dfn(i,4) 
        tmp = dfn(i,3);
        dfn(i,3) = dfn(i,1);
        dfn(i,1)= tmp;
        tmp = dfn(i,2);
        dfn(i,2) = dfn(i,4);
        dfn(i,4)= tmp;
        end
    end
end

% Specify you conditions
TF1 = dfn(:,1) < 0.0; 
TF2 = dfn(:,2) < 0.0;
TF3 = dfn(:,1) > len(1); 
TF4 = dfn(:,2) > len(2);
TF5 = dfn(:,3) < 0.0; 
TF6 = dfn(:,4) < 0.0;
TF7 = dfn(:,3) > len(1); 
TF8 = dfn(:,4) > len(2);

% combine them
TFall = TF1 | TF2 | TF3 | TF4 | TF5 | TF6 | TF7 | TF8;
% remove
dfn(TFall,:) = [];

% Cut the fractures so that they don't intersect with the lower and west
% boundaries (TESTING PURPOSE)
dfn(:,1) = max(dfn(:,1),0.01*len(1));
dfn(:,2) = max(dfn(:,2),0.01*len(2));
dfn(:,3) = max(dfn(:,3),0.01*len(1));
dfn(:,4) = max(dfn(:,4),0.01*len(2));
dfn(:,1) = min(dfn(:,1),0.99*len(1));
dfn(:,2) = min(dfn(:,2),0.99*len(2));
dfn(:,3) = min(dfn(:,3),0.99*len(1));
dfn(:,4) = min(dfn(:,4),0.99*len(2));


N_fractures = length(dfn);
Nf_i = [];
dxf = sqrt(dx(1)*dx(1)+dx(2)*dx(2));

xc = []; yc = [];
xb = []; yb = [];
xe = []; ye = [];
x = [];y = [];

beta =  atan2((dfn(:,4)-dfn(:,2)),(dfn(:,3)-dfn(:,1)))*180/pi;
L = sqrt( (abs(dfn(:,3)-dfn(:,1))).^2 + (abs(dfn(:,4)-dfn(:,2))).^2 ) ;
L = round(L./dxf);

dfn(:,3) = dfn(:,1)+cos(beta*pi/180).*L.*dxf;
dfn(:,4) = dfn(:,+2)+sin(beta*pi/180).*L.*dxf;

j = 1; Nfi =0;
for i=1:N_fractures 
    if L(i) > 8
        xi  = linspace(dfn(i,1),dfn(i,3),L(i)+1);
        yi  = linspace(dfn(i,2),dfn(i,4),L(i)+1);
        dxi = xi(2:L(i)) - xi(1:L(i)-1);
        dyi = yi(2:L(i)) - yi(1:L(i)-1);
        if (any(round(sqrt(dxi.*dxi + dyi.*dyi),4) ~= round(dxf,4)))
            error('Error in dxf')
        end

        Nfi = Nfi + L(i)-1;

        xb = [xb xi(1:L(i)-1)];
        xe = [xe xi(2:L(i))];
        yb = [yb yi(1:L(i)-1)];
        ye = [ye yi(2:L(i))];
        if (i < N_fractures)
            if (sqrt( (abs(dfn(i+1,1)-dfn(i,3))).^2 + (abs(dfn(i+1,2)-dfn(i,4))).^2 ) > 1e-3*max(len(1),len(2))) ;
            Nf_i(j) = Nfi;
            j = j+1;
            Nfi = 0;
            end
        end
    end
end
N_fractures = j-1;

Nf_f = length(xb);

frac_angle = atan2((ye-yb),(xe-xb))*180/pi;

XY1 = [xb' yb' xe' ye'];

% Comment in the following code to visualize the fracture network before
% the simulation is started
% -------------------------------------------------------------------------
% ListOfVariables = who;
% for k = 1:length(ListOfVariables)
%   assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
% end 
% figure(221)
% line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','[0 0.4470 0.7410]','LineWidth',1);
% hold on
% scatter(xe,ye,30,'filled')
% scatter(xb,yb,30,'filled')
% 
% hold off
% xlabel('x [m]');
% ylabel('y [m]');
% axis tight