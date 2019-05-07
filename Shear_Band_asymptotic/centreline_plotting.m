%% Centreline Plots

close all
clear all

% Load numerical data

% GOOD COMPARISON
%name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A1em2_Tp1010_Sp1010_TR1020_N201';
%name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A5em2_Tp1010_Sp1010_TR1020_N201';
%name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A1em1_Tp1010_Sp1010_TR1020_N201';

name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A1em2_Tp1010_Sp1010_TR1020_N201';

load(sprintf('output/%s.mat',name))
t = fread(fopen(sprintf('output/time_%s.bin',name)),'double');
data = fread(fopen(sprintf('output/out_%s.bin',name)),[5 numel(t)*round(N/10)],'double')';
vel = reshape(data(:,1),length(data(:,1))/length(t),length(t));
stress = reshape(data(:,2),length(data(:,2))/length(t),length(t));
gdot = reshape(data(:,3),length(data(:,3))/length(t),length(t));
g = reshape(data(:,4),length(data(:,4))/length(t),length(t));
T = reshape(data(:,5),length(data(:,5))/length(t),length(t));

% load asymptotic data
load(sprintf('output/asymptotic_%s',FILENAME))

% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');

%% Plots

% Increment for numerical plots
[~,idx] = min(abs(t_R - t));
%tt = [1:4:idx idx+2:2:numel(t)];
tt = [1:25:numel(t)];

figure(101);
hold on
plot(t(tt),T(1,tt),'bx',t_p + eps*(tau'), T_elastic + eps*T_1','k-')
plot(linspace(0,0.18,10),...
    arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),linspace(0,0.18,10)),'k-.')
plot(t_p*ones(1,2),linspace(1,1.025,2),'k--',t_R*ones(1,2),linspace(1,1.025,2),'r:')
xlabel('$t$');
ylabel('$T$');
xlim([0 0.18])
ylim([1 1.025])
legend('Numerical','Asymptotic','Elastic','Location','northwest')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_T_comparison.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');



figure(201);
hold on
plot(t(tt),stress(1,tt),'bx',t_p + eps*(tau'),s_elastic + eps*s_1','k-')
plot(linspace(0,0.18,10),...
    arrayfun(@(t) s_elastic_fun(omega,t),linspace(0,0.18,10)),'k-.')
plot(t_p*ones(1,2),linspace(0.995,1.02,2),'k--',t_R*ones(1,2),linspace(0.995,1.02,2),'r:')
xlabel('$t$');
ylabel('$s$');
xlim([0 0.18])
ylim([0.995 1.02])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_s_comparison.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

%% 

figure(301)
hold on
name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A1em2_Tp1010_Sp1010_TR1020_N201';
load(sprintf('output/%s.mat',name))
load(sprintf('output/asymptotic_%s',FILENAME))
plot(tau,T_1,'k-',tau,-s_1,'k-.')
name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A0_Tp1010_Sp1010_TR1020_N201';
load(sprintf('output/%s.mat',name))
load(sprintf('output/asymptotic_%s',FILENAME))
plot(tau,T_1,'b-',tau,-s_1,'b-.')
ylim([0 7])
xlim([-40 60])
xlabel('$\tau$')
legend('$T_1(\tau;\hat{A}_0 = 10^{-2})$',...
    '$-s_1(\tau; \hat{A}_0 = 10^{-2})$',...
    '$T_1(\tau;\hat{A}_0 = 0)$',...
    '$-s_1(\tau; \hat{A}_0 = 0)$','Location','northwest')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('fg_app1.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

%%

name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A1em2_E1e3_Tp1010_Sp1010_TR1020_N201';

load(sprintf('output/%s.mat',name))
t = fread(fopen(sprintf('output/time_%s.bin',name)),'double');
data = fread(fopen(sprintf('output/out_%s.bin',name)),[5 numel(t)*round(N/10)],'double')';
vel = reshape(data(:,1),length(data(:,1))/length(t),length(t));
stress = reshape(data(:,2),length(data(:,2))/length(t),length(t));
gdot = reshape(data(:,3),length(data(:,3))/length(t),length(t));
g = reshape(data(:,4),length(data(:,4))/length(t),length(t));
T = reshape(data(:,5),length(data(:,5))/length(t),length(t));

% load asymptotic data
load(sprintf('output/asymptotic_%s',FILENAME))

t_p1 = t_p;
t_R1 = t_R;

% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');
% Increment for numerical plots
[~,idx] = min(abs(t_R - t));
%tt = [1:4:idx idx+2:2:numel(t)];
tt = [1:4:numel(t)];

figure(501);
hold on
plot(t(tt),T(1,tt),'bx',t_p + eps*(tau'), T_elastic + eps*T_1','k-.')
figure(502);
hold on
plot(t(tt),stress(1,tt),'bx',t_p + eps*(tau'),s_elastic + eps*s_1','k-.')

name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A1em2_E2e3_Tp1010_Sp1010_TR1442_N201';

load(sprintf('output/%s.mat',name))
t = fread(fopen(sprintf('output/time_%s.bin',name)),'double');
data = fread(fopen(sprintf('output/out_%s.bin',name)),[5 numel(t)*round(N/10)],'double')';
vel = reshape(data(:,1),length(data(:,1))/length(t),length(t));
stress = reshape(data(:,2),length(data(:,2))/length(t),length(t));
gdot = reshape(data(:,3),length(data(:,3))/length(t),length(t));
g = reshape(data(:,4),length(data(:,4))/length(t),length(t));
T = reshape(data(:,5),length(data(:,5))/length(t),length(t));

% load asymptotic data
load(sprintf('output/asymptotic_%s',FILENAME))

% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');
% Increment for numerical plots
[~,idx] = min(abs(t_R - t));
%tt = [1:4:idx idx+2:2:numel(t)];
tt = [1:8:numel(t)];


figure(501);
hold on
plot(t(tt),T(1,tt),'b^',t_p + eps*(tau'), T_elastic + eps*T_1','k-.')
plot(t_p1*ones(10),linspace(1.005,1.03,10),'k--',t_R1*ones(10),linspace(1.005,1.03,10),'r:')
xlabel('$t$');
ylabel('$T$');
xlim([0.08 0.22])
ylim([1.005 1.03])
legend('Numerical $\hat{E} = 10^3$','Asymptotic $\hat{E} = 10^3$',...
    'Numerical $\hat{E} = 2 \times 10^3$','Asymptotic $\hat{E} = 2 \times 10^3$',...
    'Location','northwest')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_T_Ecomparison.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

figure(502);
hold on
plot(t(tt),stress(1,tt),'b^',t_p + eps*(tau'),s_elastic + eps*s_1','k-.')
xlabel('$t$');
ylabel('$s$');
xlim([0.08 0.22])
%ylim([0.995 1.005])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_s_Ecomparison.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');



