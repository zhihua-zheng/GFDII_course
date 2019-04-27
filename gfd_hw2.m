
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
% plot(M,z)

%% Constant Nsq case

n      = (1:4);
lambda = (n*pi/N(1)/H).^2;

c_grav = 1/sqrt(lambda);

figure
stem(c_grav)

%% Eigenvalues & eigenvectors

A = ( zeros(I) + diag([-ones(I-2,1);0],-1) + diag([0;ones(I-2,1)],1)) ./ (2*dz);
B = (-2*eye(I) + diag([ ones(I-2,1);2],-1) + diag([2;ones(I-2,1)],1)) ./ (dz^2);

Phi    = zeros(I,I,J);
Lambda = zeros(I,I,J);

for j = 1:J
        
    Lj = diag(B*M(:,j),0);
    Mj = diag(M(:,j),0);
    C  = -(A*Lj + B*Mj);
    
    [Phi(:,:,j),Lambda(:,:,j)] = eig(C);
end
