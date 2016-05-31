function [Z,Aijk] = intersectionsGrid(xd,yd,XY)
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

x = linspace(0,len(1),Nf(1)+1);                                            % x grid vector (used in plots)
y = linspace(0,len(2),Nf(2)+1);                                            % y grid vector (used in plots)
[Xq,Yq] = ndgrid(xd,yd);

%% // example for first line
n = length(XY);
xl = [XY(:,1) XY(:,3)];
yl = [XY(:,2) XY(:,4)];

dims = size(Xq);
Z = zeros(size(Xq));
Aijk = zeros(dims(1), dims(2), n);

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

    lengths = zeros(numel(Xq'),1);
    dm = zeros(numel(Xq'),1);
    for j = 1:length(inter)-1
        midx = 0.5*(inter(j,1)+inter(j+1,1));
        midy = 0.5*(inter(j,2)+inter(j+1,2));

        % get segments lengths
        dxx  = (inter(j,1)-inter(j+1,1)).^2;
        dyy  = (inter(j,2)-inter(j+1,2)).^2;
        l = sqrt(dxx + dyy);   
        
        d = (midx-Xq).^2+(midy-Yq).^2; %// compute squared distances
        [~, ind] = min(d(:)); %// minimize distance and obtain (linear) index of minimum
        lengths(ind) = lengths(ind) + l;
        
        % Compute mean distances <d>
        % -----------------------------------------------------------------
        if (~ishorz && ~isvert)
            %Compute the line through the intersections
            x1  = [inter(j,1) inter(j+1,1)];
            y1  = [inter(j,2) inter(j+1,2)];
            
            A = zeros(4,1);
            d_mean = zeros(4,1);
            d = sqrt((x1(1)-x1(2)).^2 +(y1(1)-y1(2)).^2);
            if (d ~= 0)
                p  = polyfit(x1,y1,1);

                dxk = [0 0]; dyk = [0 0];
                orientation = [0 0];
                x2 = [0 0]; y2 = [0 0];
                
                for k = 1:2
                    l = j-1+k;
                    %Min x grid distance
                    d = (inter(l,1)-x).^2;
                    d = d(:);
                    [valx,indx] = min(d);
                    % Min y grid distance
                    d = (inter(l,2)-y).^2;
                    d = d(:);
                    [valy,indy] = min(d);

                    %Check if intersection is exactly on the corner
                    if ( indx == indy)
                        warning('Exactly corner! Dont worry, it still works.');
                    end

                    %Check if x or y grid intersection is closer
                    if (valx < valy) 
                        % The intersection is on the vertical
                        x2(k) = x(indx); y2(k)=p(1)*x2(k)+p(2);
                        orientation(k) = 1;
                    else
                        % The intersection is on the horizontal
                        y2(k)=y(indy); x2(k) = (y2(k)-p(2))/p(1);
                    end
                    
                    dxk(k) = abs(x2(k)-x1(k));
                    dyk(k) = abs(y2(k)-y1(k));
                    % Calculate small A
                    A(k+2) = 0.5*(abs(x1(k)-x2(k)) * abs(y1(k)-y2(k)) );
                    d_mean(k+2) = calc_d_mean(sqrt((x1(k)-x2(k)).^2),sqrt((y1(k)-y2(k)).^2));

                end
                % Compute the dxk and dyk based on the information
                if orientation(1) % First intersection is vertical
                    dy1 = dx(2) - dyk(1);
                    dx2 = dx(1) + dxk(1);
                else              % First intersection is horizontal
                    dy1 = dx(2) + dyk(1);
                    dx2 = dx(1) - dxk(1);
                end
                if orientation(2) % Second intersection is vertical
                    dx1 = dx(1) + dxk(2);
                    dy2 = dx(2) - dyk(2);
                else              % Second intersection is horizontal
                    dx1 = dx(1) - dxk(2);
                    dy2 = dx(2) + dyk(2);
                end
                    
                % Calculate big A
                A(1) = 0.5 * abs(dx1)*abs(dy1);
                A(2) = 0.5 * abs(dy1)*abs(dy2);
                d_mean(1) = calc_d_mean(sqrt(dx1.^2),sqrt(dy1.^2));
                d_mean(2) = calc_d_mean(sqrt(dx2.^2),sqrt(dy2.^2));
                
                dm(ind) =  (A(1)*d_mean(1)+A(2)*d_mean(2)-A(3)*d_mean(3)-A(4)*d_mean(4))/(A(1)+A(2)-A(3)-A(4));
            end
            if(dm(ind) < 0)
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
                dm(ind) = (dy1.^2 + dy2.^2)./(2.*dx(2));
            elseif isvert
                d = (inter(j,1)-x).^2;
                d = d(:);
                [~,indx] = min(d);
                dx1 = abs(inter(j,1) - x(indx));
                dx2 = dx(1) - dx1;
                dm(ind) = (dx1.^2 + dx2.^2)./(2.*dx(1));
            end
        end
        
    end
    CIij = lengths ./ dm;
    CIij(isnan(CIij)) = 0 ;
    Aij = zeros(size(Xq));
    Aij(:)= Aij(:) +  CIij; 
    Aijk(:,:,i) = Aij;    
    Z(:) = logical(Z(:) + Aij(:));
end

% Comment in the following code to visualize the intersections between the
% fracture network and the matrix grid before the simulation is started
% -------------------------------------------------------------------------
% figure(2001)
% hold on
% test = sum(Aijk,3);
% imagesc(xd,yd,test');
% %scatter(Xq(:),Yq(:),'k*');
% scatter(X(:),Y(:),'k.');
% scatter(intersections2(:,1),intersections2(:,2),'go');
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
