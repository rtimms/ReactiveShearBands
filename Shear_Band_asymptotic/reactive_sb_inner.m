% Solution of centreline temperature for reactive shear band using an
% elastic-plastic constitutive law

close all
clear all

%% Model set-up

% Non-dimensional parameters
rho_hat = 1e-4;                           % Density 
omega = 1e-1;                             % Rescaled nominal strain-rate
lambda = 1e-3;                            % Plastic work

% Plastic strain-rate parameters
T_p = 1.01;
s_p = 1.01;
beta1 = 1;
beta2 = 1;
gammadot0 = 1;

% Reaction Parameters
T_R = 1.02;
Omega_hat = 1;
A0_hat = 5e-2;
E_hat = 1e3;

% Small parameter 
eps = T_R^2/E_hat;

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

%% Solve plastic problem 

% Initial time 
eta_0 = -beta3*t_p/eps + log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2)));

% Run solver until T = T_R
[f,g,eta] = fg_solver_EPreactive(h_flux,beta1,beta2,beta3,rho_hat_0,lambda,...
    omega,gammadot0,a,b,Lambda_p,0,0,eps,delta,t_p,T_R,0,eta_0,0);

% Convert back to tau timescale
tau = (eta - log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2))))/beta3;
T_1 = f/beta1;
s_1 = -g/beta2;


%% Reaction stage

% Critical reaction time
tau_R = tau(end);

% Stress at onset of reaction
s_R = s_elastic_fun(omega,t_p + eps*tau_R) + eps*s_1(end);

% Solve ODEs for inner temperature and stress
options=odeset('Stats','on','RelTol',1e-6);
[tau_hat,theta_psi] = ode15s(@ODE_reactive_sb_inner,[-(t_p+eps*tau_R)/(eps^(1/2)) 1/(10*eps)],[0 0],options,...
    beta1,beta2,beta3,gammadot0,lambda,A0_hat,Omega_hat,tau_R,s_R,T_1(end),s_1(end));

% Reaction time
t_R = t_p + eps*tau_R;

save('Approach_II.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%{
%%%%%%%%%%%%%%%%%     --  Plot Temperature/Stress History  --  %%%%%%%%%%%%
figure(1)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau') + eps*T_1','r-.','linewidth',2)
plot(t_p + eps*(tau_R) + eps^(3/2)*(tau_hat), T_R + eps*theta_psi(:,1),'r','linewidth',2)
xlabel('t','interpreter','latex')
ylabel('T','interpreter','latex')
set(gca,'fontsize',18)
box on

figure(2)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau')  + eps*s_1','b-.','linewidth',2)
plot(t_p + eps*(tau_R) + eps^(3/2)*(tau_hat),s_R + eps*theta_psi(:,2),'b','linewidth',2)
xlabel('t','interpreter','latex')
ylabel('s','interpreter','latex')
set(gca,'fontsize',18)
box on

figure(3)
plot(t_p + eps*(tau_R-tau(1)) + eps^(3/2)*(tau_hat-tau_hat(1)), T_R + eps*theta_psi(:,1),'r','linewidth',2)
xlabel('t','interpreter','latex')
ylabel('$T_R + \theta$','interpreter','latex')
set(gca,'fontsize',18)
box on

figure(4)
plot(t_p + eps*(tau_R-tau(1)) + eps^(3/2)*(tau_hat-tau_hat(1)),s_R + eps*theta_psi(:,2),'b','linewidth',2)
xlabel('t','interpreter','latex')
ylabel('$s_R + \psi$','interpreter','latex')
set(gca,'fontsize',18)
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%%

%%%%%%%%%%%%%%%%%     --  Plot Temperature/Stress History  --  %%%%%%%%%%%%

ts = (t_p + eps*(tau(end))) - (t_p + eps*(tau_R) + eps^(3/2)*(tau_hat(1)));

figure(11)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau') + eps*T_1','k-.')
plot(ts + t_p + eps*(tau_R) + eps^(3/2)*(tau_hat), T_R + eps*theta_psi(:,1),'k')
xlabel('$t$')
ylabel('$T$')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('T_match.tikz', 'showInfo', false, ...
     'parseStrings',false,'standalone', false, ...
     'height', '\figureheight', 'width','\figurewidth');

figure(12)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau')  + eps*s_1','k-.')
plot(ts + t_p + eps*(tau_R) + eps^(3/2)*(tau_hat),s_R + eps*theta_psi(:,2),'k')
xlabel('$t$')
ylabel('$s$')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('s_match.tikz', 'showInfo', false, ...
     'parseStrings',false,'standalone', false, ...
     'height', '\figureheight', 'width','\figurewidth');

figure(13)
plot(tau_hat,theta_psi(:,1),'k')
xlabel('$\hat{\tau}$')
ylabel('$\theta$')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('theta_soln.tikz', 'showInfo', false, ...
     'parseStrings',false,'standalone', false, ...
     'height', '\figureheight', 'width','\figurewidth');

figure(14)
plot(tau_hat,theta_psi(:,2),'k')
xlabel('$\hat{\tau}$')
ylabel('$\psi$')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('s_soln.tikz', 'showInfo', false, ...
     'parseStrings',false,'standalone', false, ...
     'height', '\figureheight', 'width','\figurewidth');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% comapre with single stage approach

t_R2 = t_R;


ts = (t_p + eps*(tau(end))) - (t_p + eps*(tau_R) + eps^(3/2)*(tau_hat(1)));
figure(101)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau') + eps*T_1','k-.')
plot(1000,1000,'k-')
legend('Asymptotic (3 stage)','Asymptotic','Location','northwest')
plot(ts + t_p + eps*(tau_R) + eps^(3/2)*(tau_hat), T_R + eps*theta_psi(:,1),'k-.')

figure(201)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau')  + eps*s_1','k-.')
plot(ts + t_p + eps*(tau_R) + eps^(3/2)*(tau_hat),s_R + eps*theta_psi(:,2),'k-.')


name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A5em2_Tp1010_Sp1010_TR1020_N201';
load(sprintf('output/%s.mat',name))
load(sprintf('output/asymptotic_%s',FILENAME))
% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');


figure(101)
plot(t_p + eps*(tau'),T_elastic + eps*T_1','k-')
plot(t_p*ones(10),linspace(1.01,1.025,10),'k--',t_R*ones(10),linspace(1.01,1.025,10),'r--',t_R2*ones(10),linspace(1.01,1.025,10),'r:')
box on
ylim([1 1.025])
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_T_3stagecomparison.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
xlim([0.12 0.18])
ylim([1.01 1.025])
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_T_3stagecomparison_zoom.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

figure(201)
plot(t_p + eps*(tau'),s_elastic + eps*s_1','k-')
plot(t_p*ones(10),linspace(0.995,1.005,10),'k--',t_R*ones(10),linspace(0.995,1.005,10),'r--',t_R2*ones(10),linspace(0.995,1.005,10),'r:')
box on
ylim([0.995 1.005])
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_s_3stagecomparison.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
xlim([0.12 0.18])
ylim([0.995 1.005])
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_s_3stagecomparison_zoom.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
