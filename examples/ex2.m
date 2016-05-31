addpath ../

EDFM('InputEx2',0)


figure 
hold on;
pcolor(x,y,p'); shading interp
axis square
[C,hfigc] = contour(x, y, p','ShowText','off');
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [0 0 0]);
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r','LineWidth',2);
hold off;

c=colorbar
ylabel(c,'Pressure [Pa]')
xlabel('x [m]')
ylabel('x [m]')

figure 
hold on;
pcolor(x,y,sNew'); shading interp
axis square
[C,hfigc] = contour(x, y, sNew',5,'ShowText','off');
set(hfigc, ...
    'LineWidth',1.0, ...
    'Color', [0 0 0]);
line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r','LineWidth',2);
hold off;

c=colorbar
ylabel(c,'Concentration [kg/m^3]')
xlabel('x [m]')
ylabel('x [m]')

rmpath ../
