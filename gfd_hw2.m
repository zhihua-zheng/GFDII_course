
%% General setup

clear
my_dir = '~/Documents/Courses/OCN513/data/';
N      = flipud(load([my_dir,'N.txt'])); % [1/s]
Nsq    = N.^2; % [1/s^2]
M      = 1 ./ Nsq; % [s^2]
[I,J]  = size(M); % # of vertical levels, # of profiles

H      = 4000; % [m]
z      = linspace(-H,0,I)';
dz     = z(2) - z(1);

%% Nsq profiles

plot(Nsq,z)
% plot(M(:,1),z)

%% Constant Nsq case

n      = (1:4);
lambda = (n .* pi/N(1)/H).^2;

c_grav = 1 ./ sqrt(lambda);

lat = (5:60); % latitude
f   = gsw_f(lat); % Coriolis parameter [1/s]
Ld  = c_grav(1) ./ f;

R      = 6400*1e3; % Earth's radius [m]
Omega  = 7.292115e-5; % (Groten, 2004) Earth's rotation rate [1/s]
beta_f = 2*Omega*cosd(lat) ./ R; % beta parameter for Coriolis effect [1/(m*s)]
c_Ro   = -beta_f ./ (f.^2) .* c_grav(1).^2; % phase speed for long Rossby wave

dis = 100*111*cosd(lat)*1e3; % distance across Pacific ocean (100 degree) [m]
t   = dis ./ abs(c_Ro); % time [s]

figure('position',[0 0 650 1080])
subplot(4,1,1)
stem(n,c_grav,'LineWidth',2)
set(gca,'fontsize',12,'TickLabelInterpreter','latex')
xlabel('vertical mode number','Interpreter','latex','fontsize',14)
ylabel('gravity wave $c_g$ [m/s]','Interpreter','latex','fontsize',13)
xlim([0 5])
ylim([0 15])
xticks(1:4)
% saveas(gcf,'./figs/c_grav','png')

subplot(4,1,2)
plot(lat,Ld/1000,'LineWidth',2)
set(gca,'fontsize',12,'TickLabelInterpreter','latex')
xlabel('latitude [degree]','Interpreter','latex','fontsize',14)
ylabel('deformation radius $L_d$ [km]','Interpreter','latex','fontsize',13)
legend('first baroclinic mode','Interpreter','latex','fontsize',14,'location','east')
xlim([5 60])
ylim([0 1200])
% saveas(gcf,'./figs/Ld','png')

subplot(4,1,3)
plot(lat,c_Ro,'LineWidth',2)
set(gca,'fontsize',12,'TickLabelInterpreter','latex')
xlabel('latitude [degree]','Interpreter','latex','fontsize',14)
ylabel('long Rossby wave $c_{Ro}$ [m/s]','Interpreter','latex','fontsize',13)
legend('first baroclinic mode','Interpreter','latex','fontsize',14,'location','east')
xlim([5 60])
% saveas(gcf,'./figs/c_Ro','png')

subplot(4,1,4)
plot(lat,t/3600/24,'LineWidth',2)
set(gca,'fontsize',12,'TickLabelInterpreter','latex')
title('time for long Rossby wave to across the Pacific','Interpreter','latex','fontsize',15)
xlabel('latitude [degree]','Interpreter','latex','fontsize',14)
ylabel('time [day]','Interpreter','latex','fontsize',13)
legend('first baroclinic mode','Interpreter','latex','fontsize',14,'location','southeast')
xlim([5 60])

saveas(gcf,'./figs/hw2/const_N','png')

%% Horizontal velocities / pressure - multiple N profiles

A = ( zeros(I) + diag([-ones(I-2,1);0],-1) + diag([0;ones(I-2,1)],1)) ./ (2*dz);
B = (-2*eye(I) + diag([ ones(I-2,1);2],-1) + diag([2;ones(I-2,1)],1)) ./ (dz^2);

% eigen-pairs for horizontal velocities, pressure...
Phi    = zeros(I,I,J); % eigenfuncions
Lambda = zeros(I,I,J); % eigenvalues

for j = 1:J
        
    Lj = diag(A*M(:,j));
    Mj = diag(M(:,j));
    C  = -(Lj*A + Mj*B);
    
    [phi,lambda] = eig(C);
    
    % sort the eigen-pairs
    [~,ind]       = sort(diag(lambda));    
    Lambda(:,:,j) = lambda(ind,ind);
    Phi(:,:,j)    = phi(:,ind);
end

%% Visualize phi modes - multiple N profiles

figure('position',[0 0 600 740])

% Number of rows and columns of axes
ncols = 5;
nrows = J;

row = zeros(ncols,nrows);
col = zeros(ncols,nrows);

% width and height of each subplot axis in normalized units
% axisw = (1 / ncols) * 0.95;
% axish = (1 / nrows) * 0.95;

for k = 1:ncols*nrows % count for sub-figures

    % calculate the row and column of the subplot
    row(k) = ceil( k/ncols );
    col(k) = mod( k-1, ncols ) + 1;
    
    % calculate the left, bottom coordinate of this subplot
%     axisl = (axisw+0.02) * (col-1);
%     axisb = (axish+0.02) * (row-1);

    subplot(J,5,k)

    if col(k)==1
        plot(Nsq(:,row(k)),z/1000,'linewidth',2,'color','k')
        grid on
        xlim([0 4e-4])
        set(gca,'fontsize',10,'TickLabelInterpreter','latex','GridLineStyle','--')
        
        if row(k) ==1
            title('$N^2$ [$1/s^2$]','Interpreter','latex','fontsize',14)
        end
    else
        phi = Phi(:,col(k),row(k));
        phi_l = max(abs(phi)); % max magnitude
        plot(phi/phi_l,z/1000,'linewidth',1.5,'color','r')
        hold on
        plot([0 0],[-4 0],'linewidth',1.8,'color',[.3 .3 .3],'linestyle',':')
        hold off
        grid on
        xlim([-1 1])
        set(gca,'fontsize',10,'TickLabelInterpreter','latex','GridLineStyle','--')
        
        if row(k) ==1
            title(['$\phi_{',num2str(col(k)-1),'}$'],'Interpreter','latex','fontsize',14)
        end
    end
end

[~,hy] = suplabel('z [km]','y'); 
set(hy,'Interpreter','latex','fontsize',15)

saveas(gcf,'./figs/hw2/phi_modes','png')

%% Fisrt baroclinic mode for phi - multiple N profiles

k      = find(lat==40);
c_grav = 1 ./ sqrt(squeeze(Lambda(2,2,:)));
Ld     = c_grav ./ f(k);
c_Ro   = -beta_f(k)/(f(k)^2) .* (c_grav.^2);
t      = 100*111*cosd(lat(k))*1e3 ./ abs(c_Ro);

%% Vertical velocity - multiple N profiles

% eigen-pairs for vertical velocity
S     = zeros(I,I,J); % eigenfuncions
Gamma = zeros(I,I,J); % eigenvalues

B = (diag([0;-2*ones(I-2,1);0]) + diag([ones(I-2,1);0],-1) + diag([0;ones(I-2,1)],1)) ./ (dz^2);

for j = 1:J
        
    Mj = diag(M(:,j));
    C  = -Mj*B;
    
    [s,gamma] = eig(C);
    
    % sort the eigen-pairs
    [~,ind]      = sort(diag(gamma));    
    Gamma(:,:,j) = gamma(ind,ind);
    S(:,:,j)     = s(:,ind);
end

%% Visualize s modes - multiple N profiles

figure('position',[0 0 300 840])

% Number of rows and columns of axes
ncols = 2;
nrows = J;

row = zeros(ncols,nrows);
col = zeros(ncols,nrows);

for k = 1:ncols*nrows % count for sub-figures

    % calculate the row and column of the subplot
    row(k) = ceil( k/ncols );
    col(k) = mod( k-1, ncols ) + 1;
    
    subplot(J,2,k)

    if col(k)==1
        plot(Nsq(:,row(k)),z/1000,'linewidth',2,'color','k')
        grid on
        xlim([0 4e-4])
        set(gca,'fontsize',10,'TickLabelInterpreter','latex','GridLineStyle','--')
        ylabel('z [km]','Interpreter','latex','fontsize',14)
        if row(k) ==1
            title('$N^2$ [$1/s^2$]','Interpreter','latex','fontsize',14)
        end
    else
        s = S(:,col(k)+1,row(k));
        s_l = max(abs(s)); % max magnitude
        plot(s/s_l,z/1000,'linewidth',1.5,'color','b')
        hold on
        plot([0 0],[-4 0],'linewidth',1.8,'color',[.3 .3 .3],'linestyle',':')
        hold off
        grid on
        xlim([-1 1])
        set(gca,'fontsize',10,'TickLabelInterpreter','latex','GridLineStyle','--')
        
        if row(k) ==1
            title(['$S_{',num2str(col(k)-1),'}$'],'Interpreter','latex','fontsize',14)
        end
    end
end

saveas(gcf,'./figs/hw2/s_modes','png')
