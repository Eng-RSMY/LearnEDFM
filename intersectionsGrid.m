function [CI] = intersectionsGrid(xd,yd,XY)
%  Finds the intersections of the fracture segments with the grid and
%  calculates the connectivity index CI between fracture and matrix grids
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
%
%  grid_intersections(xd,yd,XY)
%
%  Input: 
%        xd      (1, nx)         x-direction grid cell centers
%        yd      (1, ny)         y-direction grid cell centers
%        XY      (nf, 4)         array containing all fracture segments as
%                               [x_begin y_begin x_end y_end]
%
%  Output:
%        Aijk   (nx*ny,nx*ny,nf)fracture-grid connectivity indix array
%        Z      (nx*ny,nx*ny)   fracture-grid intersections (boolean array)


global dx len Nf
x = linspace(0,len(1),Nf(1)+1); 
y = linspace(0,len(2),Nf(2)+1);


[Xq,Yq] = ndgrid(xd,yd);

%% // example for first line
n = length(XY);
xl = [XY(:,1) XY(:,3)];
yl = [XY(:,2) XY(:,4)];

CI = [];

% intersections2 = [];% This is used only in the
%                       optional plot at the end of the function

for i = 1:n
    xli = xl(i,:);
    yli = yl(i,:);
    
    isvert = 0;
    ishorz = 0;
    % check for vertical lines
    if  (all(xli == xli(1)))
        yi = NaN(size(x));
        isvert = 1;
    else
        % reinterp the Y values over the X-Grid defining the domain
        yi = qinterp1( xli , yli , x ) ;
        yi = yi';
    end
    % check for horizontal lines
    if (all(yli == yli(1)))
        xi = NaN(size(y));
        ishorz = 1;
    else
        % reinterp the X values over the Y-Grid defining the domain
        xi = qinterp1( yli , xli , y ) ;
        xi = xi';
    end
    
    i1 = [xi x xli]';
    d1 = [y yi yli]';
    intersections = [i1 d1];
    inter = intersections(~any(isnan(intersections),2),:);
    
    if isvert
        % sort by y axis values
        [~,ii] = sort(inter(:,2));
    else
        % sort by x axis values
        [~,ii] = sort(inter(:,1));
    end
    inter = inter(ii,:);
    
    %intersections2 = [intersections2; inter]; % This is used only in the
    %optional plot at the end of the function

    lengths = zeros(length(inter)-1,1);
    dm = zeros(length(inter)-1,1);
    inds = zeros(length(inter)-1,1);
    for j = 1:length(inter)-1
        midx = 0.5*(inter(j,1)+inter(j+1,1));
        midy = 0.5*(inter(j,2)+inter(j+1,2));

        % get segments lengths
        dxx  = (inter(j,1)-inter(j+1,1)).^2;
        dyy  = (inter(j,2)-inter(j+1,2)).^2;
        
        % get the index of the cell that this segment belongs to
        d = (midx-Xq).^2+(midy-Yq).^2; %// compute squared distances
        [~, ind] = min(d(:)); %// minimize distance and obtain (linear) index of minimum
        
        lengths(j) = sqrt(dxx + dyy);   
        inds(j) = ind;
        % Compute mean distances <d>
        % -----------------------------------------------------------------
        if (~ishorz && ~isvert)

                dm(j) =  calc_d_mean(sqrt(dx(1).^2),sqrt(dx(2).^2));
            if(dm(j) < 0)
                ListOfVariables = who;
                for k = 1:length(ListOfVariables)
                    assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
                end 
                error('l <0 in d_mean calculation')
            end
            
        else
            if ishorz
                d = (inter(j,2)-y).^2;
                d = d(:);
                [~,indy] = min(d);
                dy1 = abs(inter(j,2) - y(indy));
                dy2 = dx(2) - dy1;
                dm(j) = (dy1.^2 + dy2.^2)./(2.*dx(2));
            elseif isvert
                d = (inter(j,1)-x).^2;
                d = d(:);
                [~,indx] = min(d);
                dx1 = abs(inter(j,1) - x(indx));
                dx2 = dx(1) - dx1;
                dm(j) = (dx1.^2 + dx2.^2)./(2.*dx(1));
            end
        end
        
    end
    CIij = lengths./ dm;
    CIij(isnan(CIij)) = 0 ;
    CI = vertcat(CI, [inds repmat(i,[length(inds) 1]) CIij]);
end

% % Comment in the following code to visualize the intersections between the
% % fracture network and the matrix grid before the simulation is started
% % -------------------------------------------------------------------------
% figure(2001)
% hold on
% test = sum(Aijk,3);
% pcolor(xd,yd,test');
% %scatter(Xq(:),Yq(:),'k*');
% %scatter(X(:),Y(:),'k.');
% %scatter(intersections2(:,1),intersections2(:,2),'go');
% scatter(xl(:),yl(:),'rx');
% line([XY(:,1)';XY(:,3)'],[XY(:,2)';XY(:,4)'],'Color','r');
% myColorMap = parula(256);
% myColorMap(1,:) = 0.8;
% colormap(myColorMap);
% colorbar
% 
% ListOfVariables = who;
% for k = 1:length(ListOfVariables)
%    assignin('base',ListOfVariables{k},eval(ListOfVariables{k}))
% end 
% 
% pause()
