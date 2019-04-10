
my_dir = '~/Documents/MATLAB/Courses/Ocean513';
load([my_dir,'/rho']);

rho  = fliplr(rho');
rho0 = 1025; % [kg/m^3]
g    = 9.81;
f    = gsw_f(40); % [radians/s]

%% Coordinates

[lz,ly] = size(rho);
y = linspace(-100,100,ly); % meridional distance, centered at 40N [km]
z = linspace(-500,0,lz); % vertical depth [m]
[Y,Z] = meshgrid(y,z);

%% Density

figure('position',[0 0 501 201])
%cmocean('thermal')

contourf(Y,Z,rho,'LineWidth',1,'LineColor',[.3 .3 .3])
caxis([1022 1028])
colorbar('TickLabelInterpreter','latex','fontsize',11)

title('Potential Density $[kg/m^3]$','Interpreter','latex','fontsize',15)
xlabel('meridional distance from 40 N [km]','Interpreter','latex','fontsize',13)
ylabel('depth [m]','Interpreter','latex','fontsize',13)
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
    
export_fig('pden','-png','-transparent','-painters','-m2')

%% Buoyancy frequency

drho_dz =  center_diff(rho,z,1);
Nsq     = -g/rho0*drho_dz; % buoyancy frequency [1/s^2]
zi      =  (z(1:end-1) + z(2:end))/2;
[Yi,Zi] =  meshgrid(y,zi);

figure('position',[0 0 501 201])
cmocean('amp')
contourf(Yi,Zi,Nsq,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([0 3.4e-4])
ylim([-500 0])

title('Buoyancy Frequency $N^2$ $[1/s^2]$','Interpreter','latex','fontsize',15)
xlabel('meridional distance from 40 N [km]','Interpreter','latex','fontsize',13)
ylabel('depth [m]','Interpreter','latex','fontsize',13)
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
    
export_fig('Nsq','-png','-transparent','-painters','-m2')

%% Geostrophic velocity

drho_dy = center_diff(rho,1000*y,2);
yj      = (y(1:end-1) + y(2:end))/2;
[Yj,Zj] = meshgrid(yj,z);
Ug      = cumsum(g/rho0/f*drho_dy,1); % geostrophic flow in x-direction [m/s]

figure('position',[0 0 501 201])
cmocean('speed')
contourf(Yj,Zj,Ug,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([0 .13])
xlim([-100 100])

title('Geostrophic velocity (eastward) $U_g$ $[m/s]$','Interpreter','latex','fontsize',15)
xlabel('meridional distance from 40 N [km]','Interpreter','latex','fontsize',13)
ylabel('depth [m]','Interpreter','latex','fontsize',13)
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
    
export_fig('Ug','-png','-transparent','-painters','-m2')

%% Relative vorticity from geostrophic flow

zeta2     =  g/rho0/f*drho_dy; % relative vorticity in y-direction [1/s]
zeta3     = -center_diff(Ug,y,2); % relative vorticity in z-direction [1/s]
yjj       =  (yj(1:end-1) + yj(2:end))/2;
[Yjj,Zjj] =  meshgrid(yjj,z);

figure('position',[0 0 501 482])
    subplot(2,1,1)
cmocean('balance')
contourf(Yj,Zj,zeta2,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-0.008 0.008])
xlim([-100 100])

setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
title('Relative vorticity (northward) $\omega_2$ $[1/s]$','Interpreter','latex','fontsize',15)

    subplot(2,1,2)
cmocean('balance')
contourf(Yjj,Zjj,zeta3,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-0.008 0.008])
xlim([-100 100])

setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
title('Relative vorticity (upward) $\omega_3$ $[1/s]$','Interpreter','latex','fontsize',15)

[~,hx] = suplabel('meridional distance from 40 N [km]');
[~,hy] = suplabel('depth [m]','y');
set(hx,'Interpreter','latex','fontsize',15)
set(hy,'Interpreter','latex','fontsize',15)

export_fig('zeta','-png','-transparent','-painters','-m2')

%% Potential vorticity (Ertel)

