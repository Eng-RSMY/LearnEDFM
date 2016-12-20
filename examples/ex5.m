addpath ../

EDFM('InputEx5',1)

global Nf

figure 
hold on;
pcolor(x,y,p'); shading interp
%axis square
[C,hfigc] = contour(x, y, p','ShowText','off');
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [0 0 0]);
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r','LineWidth',2);
hold off;

c=colorbar;
ylabel(c,'Pressure [Pa]')
xlabel('x [m]')
ylabel('x [m]')

figure 
hold on;
pcolor(x,y,cNew'); shading interp
%axis square
[C,hfigc] = contour(x, y, cNew',5,'ShowText','off');
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [0 0 0]);
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r','LineWidth',2);
scatter(xe,ye,30,cNewf,'filled')
hold off;

c=colorbar;
ylabel(c,'Concentration [kg/m^3]')
xlabel('x [m]')
ylabel('x [m]')

figure 
hold on;
pcolor(x,y,tNew'); shading interp
%axis square
[C,hfigc] = contour(x, y, tNew',5,'ShowText','off');
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [0 0 0]);
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r','LineWidth',2);
scatter(xe,ye)
hold off;

c=colorbar;
ylabel(c,'Temperature [°C]')
xlabel('x [m]')
ylabel('x [m]')

[vfx,vfy] = calcFractureVelocity(XY1,vf);
figure
quiver(x',y',vx(1:Nf(1),:)',vy(:,1:Nf(2))')
hold on; 
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color',[0.8 0.8 0.8],'LineWidth',4);
quiver(xb',yb',vfx,vfy,0.5,'LineWidth',1.5)
hold off
axis equal

rmpath ../
