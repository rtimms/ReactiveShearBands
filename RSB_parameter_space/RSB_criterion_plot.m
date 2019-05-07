% parameter space plot

close all
clear all

% critical eta
eta_c = 7;

Ap_vec = [1 5 10:5:100];
AR_vec = 1:30;

[Ap_plot,AR_plot] =  meshgrid(Ap_vec,AR_vec);

%% Load and plot results

figure(1)
hold on

% figure(2)
% hold on
% 
% figure(3)
% hold on

for i = 1:numel(Ap_plot)
    
    % Set parameters
    Ap = Ap_plot(i);
    AR = AR_plot(i);
    At = 0.1;
    
    % Load data
    filename = sprintf('output_NEW/fgreactive_Ap%.0f_AR%.0f_At%.0f',Ap*100,AR*100,At*100);
    load(filename)
    
    % Find eta = eta_c
    [d, ind] = min(abs(eta-eta_c));
    
    % Plot marker
    if eta(end) < eta_c
        % Runaway/SB before eta_c and computation ceased
        figure(1)
        plot(Ap,AR,'xr')
    else if f(ind) && g(ind) < 10
            % No SB
            figure(1)
            plot(Ap,AR,'ob')
        else
            % SB
            figure(1)
            plot(Ap,AR,'xr')
        end
        
    end
    
%     figure(2)
%     plot(eta,f)
%     xlim([-5 15])
%     ylim([0 10])
%     figure(3)
%     plot(eta,g)
%     xlim([-5 15])
%     ylim([0 10])
%     pause
    
end

figure(1)
box on
xlim([0 100])
ylim([0 30])
xlabel('$\Lambda_{\text{p}}$')
ylabel('$\Lambda_{\text{R}}$')
cleanfigure
matlab2tikz('PR_001.tikz', 'showInfo', false, ...
     'parseStrings',false,'standalone', false, ...
     'height', '\figureheight', 'width','\figurewidth');
%% Fit boundary

% pR = [ 11 0; 13 0.5; 15 1; 30 5; 40 7; 50 8.2; 60 9; 80 10; 90 10.1; 100 10.125];
% x = linspace(0,100,100);
% SBboundary = interp1(pR(:,1),pR(:,2),x,'spline');
% 
% figure(1)
% plot(x,SBboundary,'k-.')
% xlim([0 100])
% ylim([0 15])
% box on
% xlabel('p')
% ylabel('r')
% 
% figure(2)
% plot(x,SBboundary,'k-.')
% xlim([0 100])
% ylim([0 15])
% box on
% xlabel('p')
% ylabel('r')
% box on
% set(0,'defaulttextinterpreter','none');
% cleanfigure
% 
% matlab2tikz('PR.tikz', 'showInfo', false, ...
%      'parseStrings',false,'standalone', false, ...
%      'height', '\figureheight', 'width','\figurewidth');

