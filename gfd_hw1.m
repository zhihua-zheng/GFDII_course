
my_dir = '~/Documents/MATLAB/Courses/Ocean513';
load([my_dir,'/rho']);

rho  = fliplr(rho');
rho0 = 1025; % [kg/m^3]
g    = 9.81;
f0   = gsw_f(40); % [radians/s]


%% Coordinates

[lz,ly] = size(rho);
y = linspace(-100,100,ly); % meridional distance, centered at 40N [km]
z = linspace(-500,0,lz); % vertical depth [m]
[Y,Z] = meshgrid(y,z);
lat   = 40 + y/111.03; 
f     = gsw_f(lat); % [radians/s]

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
Ug      = g/rho0/f0*cumtrapz(z,drho_dy,1); % geostrophic flow in x-direction [m/s]

figure('position',[0 0 501 201])
cmocean('tempo')
contourf(Yj,Zj,Ug,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([0 .65])
xlim([-100 100])

title('Geostrophic velocity (eastward) $U_g$ $[m/s]$','Interpreter','latex','fontsize',15)
xlabel('meridional distance from 40 N [km]','Interpreter','latex','fontsize',13)
ylabel('depth [m]','Interpreter','latex','fontsize',13)
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
    
export_fig('Ug','-png','-transparent','-painters','-m2')

%% Relative vorticity from geostrophic flow

omega2    =  g/rho0/f0*drho_dy; % relative vorticity in y-direction [1/s]
omega3    = -center_diff(Ug,1000*y,2); % relative vorticity in z-direction [1/s]
yjj       =  (yj(1:end-1) + yj(2:end))/2;
[Yjj,Zjj] =  meshgrid(yjj,z);

figure('position',[0 0 501 482])
    subplot(2,1,1)
cmocean('balance')
contourf(Yj,Zj,omega2,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-8e-3 8e-3])
xlim([-100 100])
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
title('Relative vorticity (northward) $\omega_2$ $[1/s]$','Interpreter','latex','fontsize',15)

    subplot(2,1,2)
cmocean('balance')
contourf(Yjj,Zjj,omega3,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-3e-5 3e-5])
xlim([-100 100])
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
title('Relative vorticity (upward) $\omega_3$ $[1/s]$','Interpreter','latex','fontsize',15)

[~,hx] = suplabel('meridional distance from 40 N [km]');
[~,hy] = suplabel('depth [m]','y');
set(hx,'Interpreter','latex','fontsize',15)
set(hy,'Interpreter','latex','fontsize',15)

export_fig('omega','-png','-transparent','-painters','-m2')

%% Potential vorticity (Ertel)

% interpolate rho to zi level to match EPV2
rhoj    = (rho(:,1:end-1) + rho(:,2:end))/2;
EPV2    = omega2 .* drho_dy ./ rhoj; % meridional component of Ertel PV [1/(m*s)]

% interpolate EPV2 to zi level and yjj level to match EPV3
EPV2i   = (EPV2(1:end-1,:) + EPV2(2:end,:))/2;
EPV2ij  = (EPV2i(:,1:end-1) + EPV2i(:,2:end))/2;


% interpolate rho to zi level to match EPV3 (implicitly yjj level)
rhoijj  = (rho(1:end-1,2:end-1) + rho(2:end,2:end-1))/2;

% interpolate omega3 to zi level to match drho_dz
omega3i = (omega3(1:end-1,:) + omega3(2:end,:))/2;
EPV3    = (omega3i + f0) .* drho_dz(:,2:end-1) ./ rhoijj; % vertical component of Ertel PV [1/(m*s)]

% Add up Ertel PV components
EPV = EPV2ij + EPV3; % full Ertel PV [1/(m*s)]
[Yijj,Zijj] = meshgrid(yjj,zi);


figure('position',[0 0 501 703])
    subplot(3,1,1)
cmocean('balance')
contourf(Yj,Zj,EPV2,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-8e-10 8e-10])
xlim([-100 100])
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
title('Ertel PV component $\frac{\omega_2}{\rho}\frac{\partial\rho_{\theta}}{\partial y}$ $[m^{-1}s^{-1}]$','Interpreter','latex','fontsize',16)

    subplot(3,1,2)
cmocean('balance')
contourf(Yijj,Zijj,EPV3,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-8e-9 8e-9])
xlim([-100 100])
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
title('Ertel PV component $\frac{\omega_3 + f}{\rho}\frac{\partial\rho_{\theta}}{\partial z}$ $[m^{-1}s^{-1}]$','Interpreter','latex','fontsize',16)

    subplot(3,1,3)
cmocean('balance')
contourf(Yijj,Zijj,EPV,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-8e-9 8e-9])
xlim([-100 100])
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
title('Full Ertel PV $[m^{-1}s^{-1}]$','Interpreter','latex','fontsize',16)

[~,hx] = suplabel('meridional distance from 40 N [km]');
[~,hy] = suplabel('depth [m]','y');
set(hx,'Interpreter','latex','fontsize',16)
set(hy,'Interpreter','latex','fontsize',16)

export_fig('ErtelPV','-png','-transparent','-painters','-m2')

%% Large scale potential vorticity

lsPV  = f0*Nsq; % large scale PV [1/s^3]
lsPVe = lsPV/g; % scaled large scale PV [1/(m*s)]

figure('position',[0 0 501 482])
    subplot(2,1,1)
cmocean('balance')
contourf(Yi,Zi,lsPV,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-8e-8 8e-8])
ylim([-500 0])
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
title('large scale potential vorticity $fN^2$ $[1/s]$','Interpreter','latex','fontsize',15)

    subplot(2,1,2)
cmocean('balance')
contourf(Yi,Zi,lsPVe,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-8e-9 8e-9])
ylim([-500 0])
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
title('scaled large scale PV $fN^2/g$ $[m^{-1}s^{-1}]$','Interpreter','latex','fontsize',15)

[~,hx] = suplabel('meridional distance from 40 N [km]');
[~,hy] = suplabel('depth [m]','y');
set(hx,'Interpreter','latex','fontsize',15)
set(hy,'Interpreter','latex','fontsize',15)

export_fig('lsPV','-png','-transparent','-painters','-m2')

%% Rossby number

% L  = 200*1000;
% Ro = Ug/f0/L; % Rossby number

% interpolate omega2 to yjj level to match omega3
omega2j = (omega2(:,1:end-1) + omega2(:,2:end))/2;
omega   = sqrt(omega2j.^2 + omega3.^2); % relative vorticity [1/s]
Ro      = abs(omega3/f0); % Rossby number
clist   = [-6 -5 -4 -3 -2 -1 0];

figure('position',[0 0 501 201])
cmocean('dense')
contourf(Yjj,Zjj,log10(Ro),'LineStyle','none')
hold on
contour(Yjj,Zjj,log10(Ro),clist,'LineWidth',1,'LineColor',[.3 .3 .3])
colorbar('TickLabelInterpreter','latex','fontsize',11)
caxis([-6 1])
xlim([-100 100])

title('Rossby number $log_{10}R_o$','Interpreter','latex','fontsize',15)
xlabel('meridional distance from 40 N [km]','Interpreter','latex','fontsize',13)
ylabel('depth [m]','Interpreter','latex','fontsize',13)
setDateAxes(gca,'fontsize',11,'TickLabelInterpreter','latex')
    
export_fig('Ro','-png','-transparent','-painters','-m2')

%% Quasi-geostrophic stream function

