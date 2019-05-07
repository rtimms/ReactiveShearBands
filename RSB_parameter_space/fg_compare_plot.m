% Plots for reactive shear band results

close all
clear all

% Plot Ap = 10
Ap = 10;
AR = 0;
At = 0.01;

% Load data
filename = sprintf('output_NEW/fgreactive_Ap%.0f_AR%.0f_At%.0f',Ap*100,AR*100,At*100);
load(filename)

figure(1)
hold on
plot(eta,f,'k:','linewidth',2)
figure(2)
hold on
plot(eta,g,'k:','linewidth',2)
figure(3)
hold on
plot(eta,f,'k:','linewidth',2)
figure(4)
hold on
plot(eta,g,'k:','linewidth',2)

% Load and plot results
lineS = {'-.k','--k','-k'};

% parameter space
Ap = 10;
AR_plotvec = [1 10 15];
At = 0.01;

for i = 1:numel(AR_plotvec)
    
    
    % Load data
    filename = sprintf('output_NEW/fgreactive_Ap%.0f_AR%.0f_At%.0f',Ap*100,AR_plotvec(i)*100,At*100);
    load(filename)
    
    % Plot temperature perturbation
    figure(1)
    plot(eta,f,lineS{i},'linewidth',0.5)
    
    % Plot stress perturbation
    figure(2)
    plot(eta,g,lineS{i},'linewidth',0.5)
    
%     [Ap AR At]
%     pause
    
    
end

% parameter space
Ap = 100;
AR_plotvec = [0 1 10 15];
At = 0.01;

% Load and plot results
lineS = {':k','-.k','--k','-k'};

for i = 1:numel(AR_plotvec)
    
    
    % Load data
    filename = sprintf('output_NEW/fgreactive_Ap%.0f_AR%.0f_At%.0f',Ap*100,AR_plotvec(i)*100,At*100);
    load(filename)
    
    % Plot temperature perturbation
    figure(3)
    plot(eta,f,lineS{i},'linewidth',0.5)
    
    % Plot stress perturbation
    figure(4)
    plot(eta,g,lineS{i},'linewidth',0.5)
    
%     [Ap AR At]
%     pause
    
end







figure(1)
xlabel('$\eta$')
ylabel('$f(\eta)$')
xlim([-5 15])
ylim([0 15])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('f1.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

figure(2)
xlabel('$\eta$')
ylabel('$g(\eta)$')
xlim([-5 15])
ylim([0 15])
legend('$\Lambda_{\text{R}} = 0$','$\Lambda_{\text{R}} = 1$','$\Lambda_{\text{R}} = 10$','$\Lambda_{\text{R}} = 15$','Location','southeast')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('g1.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');


figure(3)
xlabel('$\eta$')
ylabel('$f(\eta)$')
xlim([-5 15])
ylim([0 15])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('f2.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

figure(4)
xlabel('$\eta$')
ylabel('$g(\eta)$')
legend('$\Lambda_{\text{R}} = 0$','$\Lambda_{\text{R}} = 0$','$\Lambda_{\text{R}} = 1$','$\Lambda_{\text{R}} = 10$','$\Lambda_{\text{R}} = 15$','Location','southeast')
xlim([-5 15])
ylim([0 15])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('g2.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');


