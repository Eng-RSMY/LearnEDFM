addpath ../

EDFM('InputEx1',0)

load ex1_reference_sanders.mat

xf = 2:dxf:7;
xf = xf(1:32);
figure; plot(xf,pf(33:64),'-',ref_f1e3_x,ref_f1e3_y,'x','LineWidth',1.5)
xlabel('x [m]')
ylabel('Pressure [Pa]')
legend('EDFM','Reference (Sanders,2015)')
figure; plot(x,p(:,50),'-',ref_m1e3_x,ref_m1e3_y,'x','LineWidth',1.5)
xlabel('x [m]')
ylabel('Pressure [Pa]')
legend('EDFM','Reference (Sanders,2015)')


rmpath ../
