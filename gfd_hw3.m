
%% Loading

clear
my_dir = '~/Documents/MATLAB/course_code/Ocean513/data/';
Fname  = [my_dir,'QGmodel.cdf'];

% Note the dimension in cdf is X-Y-T, but I want Y-X-Z-T in matlab

PSI(:,:,1,:) = permute(ncread(Fname,'P1'),[2 1 3]); % [m^2/s]
PSI(:,:,2,:) = permute(ncread(Fname,'P2'),[2 1 3]); % [m^2/s]

f0 = 1e-4; % [s^-1]
H  = [500 500]; % [m]
g  = 9.81; % gravity [m/s^2]
gp = 0.1; % reduced gravity [m/s^2]
kd = 1/100e3; % 1/k_d [m^-1]

alpha = gsw_alpha(35,20,0); % [1/degree C]
% rho0  = 1025; % [kg/m^3]

%% Grid

% i - x dimension 
% j - y dimension

[ny,nx,nz,nt] = size(PSI);

dt = 0.5; % [day]
dx = 10; % [km]
dy = 10; % [km]

t   = dt*(0:nt-1); % [hr]
xi  = dx*(-31:31); % [km]
yj  = dy*(-31:31); % [km]

xii = (xi(1:end-1) + xi(2:end))/2; % [km]
yjj = (yj(1:end-1) + yj(2:end))/2; % [km]

xiii = (xii(1:end-1) + xii(2:end))/2; % [km]
yjjj = (yjj(1:end-1) + yjj(2:end))/2; % [km]

x   = [-315 xii 315]; % [km]
y   = [-315 yjj 315]; % [km]

Lx   = x(end)-x(1); % [km]
Lxi  = xi(end)-xi(1); % [km]
Lxii = xii(end)-xii(1); % [km]
Ly   = y(end)-y(1); % [km]

%% Initial condition

% PSI_ini = squeeze(PSI(:,:,:,1));
% U_ini   = -center_diff(PSI_ini(:,:,1),1000*y,1);

%% Zonal averge of zonal velocity

% my center_diff can't dela with 4-D matrix
U(:,:,1,:) = -center_diff(PSI(:,:,1,:),1000*y,1);
U(:,:,2,:) = -center_diff(PSI(:,:,2,:),1000*y,1);

% zonal average
U_bar = trapz(x,U,2)/Lx;
u     = U-U_bar;

[Tj,Yj] = meshgrid(t,yj);

%% 

figure('position',[0 0 600 600])

subplot(2,1,1)
contourf(Tj,Yj,squeeze(U_bar(:,1,1,:)),6,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',28)
colorbar
caxis([-1.4 1.4])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged zonal velocity - Layer 1 [m/s]')

subplot(2,1,2)
contourf(Tj,Yj,squeeze(U_bar(:,1,2,:)),6,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',28)
colorbar
caxis([-1.4 1.4])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged zonal velocity - Layer 2 [m/s]')

% saveas(gcf,'./figs/hw3/Ubar','png')

%% Zonal average of interface height

eta     = squeeze(f0/gp*diff(PSI,1,3)); % [m]
eta_bar = trapz(x,eta,2)/Lx;

[T,Y] = meshgrid(t,y);

%%

figure('position',[0 0 600 300])

contourf(T,Y,squeeze(eta_bar),10,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',28)
colorbar
caxis([-140 140])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged interface height [m]')

%% Zonal average of potential vorticity

% lp - Laplacian operator

lpx_PSI(:,:,1,:) = center_diff(center_diff(PSI(:,:,1,:),1000*x,2),1000*x,2);
lpx_PSI(:,:,2,:) = center_diff(center_diff(PSI(:,:,2,:),1000*x,2),1000*x,2);

lpy_PSI(:,:,1,:) = center_diff(center_diff(PSI(:,:,1,:),1000*y,1),1000*y,1);
lpy_PSI(:,:,2,:) = center_diff(center_diff(PSI(:,:,2,:),1000*y,1),1000*y,1);

q = lpx_PSI(2:end-1,:,:,:) + lpy_PSI(:,2:end-1,:,:) - kd^2*(diff(PSI(2:end-1,2:end-1,:,:),1,3))/2;

q_bar = trapz(xii,q,2)/Lxii;
q_p   = q - q_bar;

[Tjj,Yjj] = meshgrid(t,yjj);

%%

figure('position',[0 0 600 600])

subplot(2,1,1)
contourf(Tjj,Yjj,squeeze(q_bar(:,1,1,:)),10,'LineStyle','none')
cmocean('curl',20)
colorbar
caxis([-3.1e-5 3.1e-5])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged potential vorticity - Layer 1 [1/s]')

subplot(2,1,2)
contourf(Tjj,Yjj,squeeze(q_bar(:,1,2,:)),10,'LineStyle','none')
cmocean('curl',20)
colorbar
caxis([-3.1e-5 3.1e-5])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged potential vorticity - Layer 2 [1/s]')

%% Zonal average of potential vorticity gradient

dq_dy(:,:,1,:) = center_diff(q(:,:,1,:),1000*y,1);
dq_dy(:,:,2,:) = center_diff(q(:,:,2,:),1000*y,1);
dq_dy_bar = squeeze(trapz(xii,dq_dy,2))/Lxii;

% dq_dx(:,:,1,:) = center_diff(q(:,:,1,:),1000*x,2);
% dq_dx(:,:,2,:) = center_diff(q(:,:,2,:),1000*x,2);
% dq_dx_bar = squeeze(trapz(xiii,dq_dx,2))/(xiii(end)-xiii(1));

[Tjjj,Yjjj] = meshgrid(t,yjjj);

%%

figure('position',[0 0 600 600])

subplot(2,1,1)
contourf(Tjjj,Yjjj,squeeze(dq_dy_bar(:,1,:)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',28)
colorbar
caxis([-1.4e-9 1.4e-9])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged potential vorticity gradient - Layer 1 [m^{-1}s^{-1}]')

subplot(2,1,2)
contourf(Tjjj,Yjjj,squeeze(dq_dy_bar(:,2,:)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',28)
colorbar
caxis([-1.4e-9 1.4e-9])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged potential vorticity gradient - Layer 2 [m^{-1}s^{-1}]')

%% Perturbation streamfunction

% zonal mean streamfunction
PSI_bar = trapz(x,PSI,2)/Lx;

% perturbation streamfunction
PSI_p = PSI - PSI_bar;

t_ind(1) = find(t==1);
t_ind(2) = find(t==6);

[X,Y] = meshgrid(x,y);

%%

figure('position',[0 0 750 600])

subplot(2,2,1)
contourf(X,Y,PSI_p(:,:,1,t_ind(1)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',22)
colorbar
caxis([-1.1e4 1.1e4])
xlabel('along-jet distance [km]')
ylabel('cross-jet distance [km]')
legend('day 1')
title('Perturbation streamfunction - Layer 1 [m^{2}s^{-1}]')

subplot(2,2,2)
contourf(X,Y,PSI_p(:,:,2,t_ind(1)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',22)
colorbar
caxis([-1.1e4 1.1e4])
xlabel('along-jet distance [km]')
ylabel('cross-jet distance [km]')
legend('day 1')
title('Perturbation streamfunction - Layer 2 [m^{2}s^{-1}]')

subplot(2,2,3)
contourf(X,Y,PSI_p(:,:,1,t_ind(2)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',22)
colorbar
caxis([-1.1e4 1.1e4])
xlabel('along-jet distance [km]')
ylabel('cross-jet distance [km]')
legend('day 6')
title('Perturbation streamfunction - Layer 1 [m^{2}s^{-1}]')

subplot(2,2,4)
contourf(X,Y,PSI_p(:,:,2,t_ind(2)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',22)
colorbar
caxis([-1.1e4 1.1e4])
xlabel('along-jet distance [km]')
ylabel('cross-jet distance [km]')
legend('day 6')
title('Perturbation streamfunction - Layer 2 [m^{2}s^{-1}]')

%% Perturbation temperature and meridional velocity

eta_p   = eta - eta_bar;

theta_p = -gp/H(1)/g/alpha*eta_p; % psudo-temperature [degree C]

% V_bar = 0
% meridional velocity perturbation

v(:,:,1,:) = center_diff(PSI_p(:,:,1,:),1000*x,2);
v(:,:,2,:) = center_diff(PSI_p(:,:,2,:),1000*x,2);

v_eta   = squeeze(mean(v,3));

[Xi,Yi] = meshgrid(xi,y);

%% 

% figure('position',[0 0 460 340])
% contourf(Xi,Yi,v_eta(:,:,t_ind(1)),10,'LineWidth',.8,'LineColor',[.2 .2 .2])
% cmocean('balance',26)
% colorbar
% caxis([-1.3e-2 1.3e-2])
% xlabel('along-jet distance [km]')
% ylabel('cross-jet distance [km]')
% legend('day 1')
% title('Meridional velocity at interface [m/s]')

figure('position',[0 0 500 240])
xlabel('along-jet distance [km]')
yyaxis left
plot(xi,v_eta(ny/2,:,t_ind(1)),'linewidth',1.5)
ylabel('meridional velocity at interface [m/s]')
yyaxis right
plot(x,theta_p(ny/2,:,t_ind(1)),'linewidth',1.5)
ylabel('temperature perturbation ')
xlim([x(1) x(end)])

%% movie

[X,Y] = meshgrid(x,y);
[Xj,Yj] = meshgrid(x,yj);
[Xiijj,Yiijj] = meshgrid(xii,yjj);

figure('position',[0 0 460 340])

for i=1:nt
    
    pcolor(X,Y,eta(:,:,1,i));
    shading interp
    text(120,300,['day ',num2str(t(i))]);
    pause(0.1)
end
    
%% Zonal average of meridional PV flux

vi   = (v(:,1:end-1,:,:) + v(:,2:end,:,:))/2; % project v to xii grid
vq_p = vi(2:end-1,:,:,:) .* q_p; % v'q' [m/s^2]

% zonal average
vq_p_bar = trapz(xii,vq_p,2)/Lxii; % <v'q'> [m/s^2]

%%

figure('position',[0 0 600 600])

subplot(2,1,1)
contourf(Tjj,Yjj,squeeze(vq_p_bar(:,1,1,:)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',26)
colorbar
caxis([-1.3e-5 1.3e-5])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged meridional PV flux - Layer 1 [m*s^{-2}]')

subplot(2,1,2)
contourf(Tjj,Yjj,squeeze(vq_p_bar(:,1,2,:)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',26)
colorbar
caxis([-1.3e-5 1.3e-5])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged meridional PV flux - Layer 2 [m*s^{-2}]')

%% Zonal average of meridional momentum flux

vj   = (v(1:end-1,:,:,:) + v(2:end,:,:,:))/2; % project v to yj grid
ui   = (u(:,1:end-1,:,:) + u(:,2:end,:,:))/2; % project u to xi grid

vu_p = vj .* ui; % v'u' [m^2/s^2]

% Zonal average
vu_p_bar = trapz(xi,vu_p,2)/Lxi;

[Tj,Yj] = meshgrid(t,yj);

%%

figure('position',[0 0 600 600])

subplot(2,1,1)
contourf(Tj,Yj,squeeze(vu_p_bar(:,1,1,:)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',34)
colorbar
caxis([-1.7e-1 1.7e-1])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged meridional momentum flux - Layer 1 [m^{2}*s^{-2}]')

subplot(2,1,2)
contourf(Tj,Yj,squeeze(vu_p_bar(:,1,2,:)),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',34)
colorbar
caxis([-1.7e-1 1.7e-1])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged meridional momentum flux - Layer 2 [m^{2}*s^{-2}]')

%% Zonal average of meridional momentum flux

eta_pi = (eta_p(:,1:end-1,:,:) + eta_p(:,2:end,:,:))/2; % project eta_p to xi grid
veta_p = v_eta .* eta_pi; % v_eta'\eta' [m^2/s]

% Zonal average
veta_p_bar = trapz(xi,veta_p,2)/Lxi;

[T,Y] = meshgrid(t,y);

%%

figure('position',[0 0 600 300])

contourf(T,Y,squeeze(veta_p_bar),8,'LineWidth',.8,'LineColor',[.2 .2 .2])
cmocean('balance',24)
colorbar
caxis([-6e1 6e1])
xlabel('time [day]')
ylabel('cross-jet distance [km]')
title('Zonally averaged meridional interface height flux [m^{2}/s]')

%%