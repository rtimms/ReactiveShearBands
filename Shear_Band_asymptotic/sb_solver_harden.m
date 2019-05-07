% Exponential elastic-plastic solution. Comparison between "cohesive"
% numerical model and exponential asymptotic solution

close all
clear all


%% Model set-up

% Non-dimensional parameters
eps = 1e-3;

rho_hat = 1e-4;                            % Density 
omega = 1e-2;                              % Rescaled nominal strain-rate
lambda = 1e-2;                             % Plastic work

% Plastic strain-rate parameters
T_p = 1.01;
s_p = 1.01;
beta1 =  1;
beta2 = 1;
beta4 = 0;
gammadot0 = 1;

% Filename
filename = sprintf('SB_harden_beta4_0_lambda1em2');
%% Centreline heat flux

% Heat flux scaling
delta = 0.1;

% Heat flux inhomogeneity (non-dimensional)
h_flux = @(t) 4*t.*exp(-4*t.^2);
dhdt = @(t)  4*(1 - 8*t.^2).*exp(-4*t.^2);

q_flux = @(t) delta*h_flux(t);


%% Elastic-platic parameters

% Critical plastic timescale
t_p = fzero(@(t) ...
    beta1*(T_p - T_elastic_centreline_fun(delta,t,h_flux)) + beta2*(s_p - s_elastic_fun(omega,t)),1e-6);

% EP parameters
a =  delta*integral(@(t) (pi*(t_p-t)).^(-1/2).*dhdt(t),1e-12,t_p);
b = delta*h_flux(t_p);
beta3 = beta1*a + beta2*omega;
rho_hat_0 = rho_hat/eps;

% Combined parameters
Lambda_p = (beta2*(rho_hat_0*beta3)^(1/2))/(beta1*lambda*(1+omega*t_p));
Lambda_k = (b*beta4)/(lambda*gammadot0*(1+omega*t_p));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Numerical solution "cohesive scheme"

% Non-dimensional elastic wave speed
S = sqrt(1/rho_hat);

% Computational grid
N = 201;                       % Number of gridpoints
y = linspace(0,5,N)';          % Spatial domain (half space) 
dy = (y(end)-y(1))/(N-1);      % Grid size
dt = dy/S;                     % Timestep

% Stopping time 
t_stop = 2*t_p;

% Run solver
[T_centreline_max] =  cohesive_solver_harden(rho_hat,S,lambda,omega,gammadot0,eps,beta1,beta2,beta4,T_p,s_p,t_p,q_flux,N,y,dy,dt,filename,t_stop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Centreline temperature -- asymptotic 

% Stopping temperature T ~ 1 + O(1) i.e. when correction grows to order 1
T_stop = 2;

% Run solver
[tau,T_1,s_1,gamma_hat] = sb_centreline_harden(h_flux,beta1,beta2,beta3,beta4,rho_hat_0,lambda,omega,gammadot0,a,b,Lambda_p,Lambda_k,eps,delta,t_p,T_stop);
save(sprintf('output/asymptotic_%s',filename))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


