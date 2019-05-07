%% Centreline Plots

close all
clear all

% Load numerical data
name = 'numerical_SB_rho1em4_omega1em2_lambda1em3_A1em2_Tp1010_Sp1010_TR1020_N201';

load(sprintf('output/%s.mat',name))
t = fread(fopen(sprintf('output/time_%s.bin',name)),'double');
data = fread(fopen(sprintf('output/out_%s.bin',name)),[5 numel(t)*round(N/10)],'double')';
stress = reshape(data(:,2),length(data(:,2))/length(t),length(t));
T = reshape(data(:,5),length(data(:,5))/length(t),length(t));

% load asymptotic data
load(sprintf('output/asymptotic_%s',FILENAME))

% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');

%% Plots

% Increment for numerical plots
% [~,idx] = min(abs(t_p - t));
% tt = [1:10:idx-10 idx-10:1:numel(t)];
tt = 1:5:numel(t);

figure(101);
hold on
plot(t(tt),T(1,tt),'bx',t_p + eps*tau(1:20:end), T_elastic(1:20:end) + eps*T_1(1:20:end)','k-')
plot(linspace(0,0.18,10),...
    arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),linspace(0,0.18,10)),'k-.')
plot(t_p*ones(1,2),linspace(1.01,1.025,2),'k--',t_R*ones(1,2),linspace(1.01,1.025,2),'r:')
xlabel('$t$');
ylabel('$T$');
xlim([0.12 0.18])
ylim([1.01 1.025])
legend('Numerical','Asymptotic','Elastic','Location','northwest')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_T_comparison_A1em2.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');


figure(201);
hold on
plot(t(tt),stress(1,tt),'bx',t_p + eps*tau(1:20:end),s_elastic(1:20:end) + eps*s_1(1:20:end)','k-')
plot(linspace(0,0.18,10),...
    arrayfun(@(t) s_elastic_fun(omega,t),linspace(0,0.18,10)),'k-.')
plot(t_p*ones(1,2),linspace(0.995,1.02,2),'k--',t_R*ones(1,2),linspace(0.995,1.02,2),'r:')
xlabel('$t$');
ylabel('$s$');
xlim([0.12 0.18])
ylim([0.995 1.005])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_s_comparison_A1em2.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');