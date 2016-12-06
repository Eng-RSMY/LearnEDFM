function [Disp] = transport_mass_Dispersion(Nf,Nf_f,dx,vx,vy)
%
%DISPERSION Creates Dispersive Matrix  
%

%-------------------------------------------------------------------------%
%
%              %-----------------------------------------------%
%              %  (c) Ivan Lunati, Univerity of Lausanne       %
%              %      rebelott@gmail.com; ivan.lunati@unil.ch  %
%              %      Ruven Kuenze, University of Lausanne     %
%              %      ruven.kunze@unil.ch                      %
%              %-----------------------------------------------%
%
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     PARAMETERS      DEFINITION      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global alphal alphat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Preparing velocities         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qx = vx./dx(2);                                                            % velocity in m/s bcause vx is in m2/s
qy = vy./dx(1);

vyx = sparse(Nf(1)+1,Nf(2));
vxy = sparse(Nf(1),Nf(2)+1);

vyx(2:Nf(1),:) = (qy(1:Nf(1)-1,1:Nf(2)) + qy(1:Nf(1)-1,2:Nf(2)+1) + qy(2:Nf(1),1:Nf(2)) + qy(2:Nf(1),2:Nf(2)+1))./4; % projext y-velocity to vertical boundaries       
vyx(1,:)       = (qy(1,1:Nf(2)) + qy(1,2:Nf(2)+1))./2;
vyx(Nf(1)+1,:) = (qy(Nf(1),1:Nf(2)) + qy(Nf(1),2:Nf(2)+1))./2;

vxy(:,2:Nf(2)) = (qx(1:Nf(1),1:Nf(2)-1) + qx(2:Nf(1)+1,1:Nf(2)-1) + qx(1:Nf(1),2:Nf(2)) + qx(2:Nf(1)+1,2:Nf(2)))./4; % projext x-velocity to horizontal boundaries
vxy(:,1)       = (qx(1:Nf(1),1) + qx(2:Nf(1)+1,1))./2;
vxy(:,Nf(2)+1) = (qx(1:Nf(1),Nf(2)) + qx(2:Nf(1)+1,Nf(2)))./2;

ux = (qx.^2 + vyx.^2).^0.5;                                                % Absolute velocity value x-direction
uy = (vxy.^2 + qy.^2).^0.5;                                                % Absolute velocity value y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Preparing dispersion coefficient   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dxx = alphal.*qx.^2./ux + alphat.*vyx.^2./ux;
Dyy = alphal.*qy.^2./uy + alphat.*vxy.^2./uy;

Dxy = (alphal - alphat).*qx.*vyx./ux;
Dyx = (alphal - alphat).*vxy.*qy./uy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Preparing dispersion matrix    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DyyN = zeros(Nf(1),Nf(2));
DyyS = zeros(Nf(1),Nf(2));
DxxW = zeros(Nf(1),Nf(2));
DxxE = zeros(Nf(1),Nf(2));

DyxN = zeros(Nf(1),Nf(2));
DyxS = zeros(Nf(1),Nf(2));
DxyW = zeros(Nf(1),Nf(2));
DxyE = zeros(Nf(1),Nf(2));

DyxNboundW = sparse(Nf(1),Nf(2)); DyxNboundE = sparse(Nf(1),Nf(2));
DyxSboundW = sparse(Nf(1),Nf(2)); DyxSboundE = sparse(Nf(1),Nf(2));
DxyWboundN = sparse(Nf(1),Nf(2)); DxyWboundS = sparse(Nf(1),Nf(2));
DxyEboundN = sparse(Nf(1),Nf(2)); DxyEboundS = sparse(Nf(1),Nf(2));


DyyN(1:Nf(1),2:Nf(2))   = dx(1)*Dyy(:,2:Nf(2))./dx(2);
DyyS(1:Nf(1),1:Nf(2)-1) = dx(1)*Dyy(:,2:Nf(2))./dx(2);
DxxE(2:Nf(1),1:Nf(2))   = dx(2)*Dxx(2:Nf(1),:)./dx(1);
DxxW(1:Nf(1)-1,1:Nf(2)) = dx(2)*Dxx(2:Nf(1),:)./dx(1);

DyxN      (2:Nf(1)-1,2:Nf(2)  ) = (1/4)*Dyx(2:Nf(1)-1,2:Nf(2)  );
DyxNboundW(  Nf(1)  ,2:Nf(2)  ) = (1/2)*Dyx(  Nf(1)  ,2:Nf(2)  );
DyxNboundE(1        ,2:Nf(2)  ) = (1/2)*Dyx(1        ,2:Nf(2)  );

DyxS      (2:Nf(1)-1,1:Nf(2)-1) = (1/4)*Dyx(2:Nf(1)-1,2:Nf(2)  );
DyxSboundW(  Nf(1)  ,1:Nf(2)-1) = (1/2)*Dyx(  Nf(1)  ,2:Nf(2)  );
DyxSboundE(1        ,1:Nf(2)-1) = (1/2)*Dyx(1        ,2:Nf(2)  );

DxyW      (1:Nf(1)-1,2:Nf(2)-1) = (1/4)*Dxy(2:Nf(1)  ,2:Nf(2)-1);
DxyWboundN(1:Nf(1)-1,1        ) = (1/2)*Dxy(2:Nf(1)  ,1        );
DxyWboundS(1:Nf(1)-1,  Nf(2)  ) = (1/2)*Dxy(2:Nf(1)  ,  Nf(2)  );

DxyE      (2:Nf(1)  ,2:Nf(2)-1) = (1/4)*Dxy(2:Nf(1)  ,2:Nf(2)-1);
DxyEboundN(2:Nf(1)  ,1        ) = (1/2)*Dxy(2:Nf(1)  ,1        );
DxyEboundS(2:Nf(1)  ,  Nf(2)  ) = (1/2)*Dxy(2:Nf(1)  ,  Nf(2)  );


Dinner = [ DyxS(:)+DxyW(:),                  DyyS(:)-DxyE(:)+DxyW(:),   -DxyE(:)-DyxS(:), ...
           DxxW(:)-DyxN(:)+DyxS(:), -DyyN(:)-DyyS(:)-DxxE(:)-DxxW(:),    DxxE(:)+DyxN(:)-DyxS(:), ...
          -DyxN(:)-DxyW(:),                  DyyN(:)+DxyE(:)-DxyW(:),    DxyE(:)+DyxN(:)];
      
Dbound = [+DyxSboundE(:)+DxyWboundN(:),... 
          -DyxSboundW(:)-DyxSboundE(:)-DxyEboundN(:)+DxyWboundN(:),...
          -DxyEboundN(:)+DyxSboundW(:),...
          -DyxNboundE(:)+DyxSboundE(:)+DxyWboundS(:)-DxyWboundN(:),...
          +DyxNboundW(:)+DyxSboundW(:)+DyxNboundE(:)-DyxSboundE(:)-DxyEboundS(:)+DxyWboundS(:)+DxyEboundN(:)-DxyWboundN(:),...
          +DxyEboundS(:)+DxyEboundN(:)-DyxNboundW(:)-DyxSboundW(:),...
          -DyxNboundE(:)-DxyWboundS(:),...
          -DyxNboundW(:)+DyxNboundE(:)-DxyWboundS(:)-DxyWboundS(:),...
          +DyxNboundW(:)-DxyWboundS(:)];

Disp = spdiags(-(Dinner + Dbound),[-Nf(1)-1,-Nf(1),-Nf(1)+1,-1,0,1,Nf(1)-1,Nf(1),Nf(1)+1],Nf(1)*Nf(2),Nf(1)*Nf(2));

Disp_f = spdiags(zeros(Nf_f,1),0,Nf_f,Nf_f);
Disp = blkdiag(Disp, Disp_f);


end
