function [tau,tau_R,T_1,s_1] = sb_centreline_EPreactive(h,beta1,beta2,beta3,rho0,lambda,omega,gammadot0,a,b,Lambda,LambdaR,Lambdat,eps,delta,t_p,T_R,T_stop)

% Time to start integration
eta0 = -beta3*t_p/eps + log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2)));

%%%%%%%%    -- Step elastic-plastic solution forward in time --     %%%%%%%

% Iterate to find eta_R
iterate = 0;
tol = 1e-3;
eta_R = 1e6;

% Initial guess for eta_R and f_R from inert solution
[f,~,eta] = fg_solver_EPreactive(h,beta1,beta2,beta3,rho0,lambda,omega,gammadot0,a,b,Lambda,0,0,eps,delta,t_p,T_R,0,eta0,0);
f_R_New = f(end);
eta_R_New = eta(end);

% figure
% plot(eta,f,'bx',eta,exp(eta),'k')
% pause

while abs(eta_R_New - eta_R) > tol;
    
    % Break if too many interations
    iterate = iterate + 1;
    if iterate > 5
        sprintf('Maximum number of iterations exceeded in searching for tau_R!')
        %pause
        break
    end
    
    eta_R = eta_R_New;
    
    % Compute up to T = T_R and check new values for f_R and eta_R
    [f,~,eta] = fg_solver_EPreactive(h,beta1,beta2,beta3,rho0,lambda,omega,gammadot0,a,b,Lambda,LambdaR,Lambdat,eps,delta,t_p,T_R,f_R_New,eta0,eta_R_New);
    
    f_R_New = f(end);
    eta_R_New = eta(end);
    
%     abs(eta_R_New - eta_R)
%     pause
end

% Run solution up to T_stop with converged values
[f,g,eta] = fg_solver_EPreactive(h,beta1,beta2,beta3,rho0,lambda,omega,gammadot0,a,b,Lambda,LambdaR,Lambdat,eps,delta,t_p,T_stop,f_R_New,eta0,eta_R_New);
% figure
% plot(eta,f,'ro',eta,g,'bo')

% Convert back to tau timescale
tau = (eta - log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2))))/beta3;
tau_R = (eta_R - log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2))))/beta3;
T_1 = f/beta1;
s_1 = -g/beta2;

end
