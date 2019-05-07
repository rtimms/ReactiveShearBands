close all
clear all

%% load plastic solution

load output/Approach_II.mat

%% Loop to find t_star

% find tau = 0 (t_p)
[~,idx] = min(abs(tau));

for i = idx:numel(tau)
   
    plastic_work = lambda*(s_elastic_fun(omega,t_p + eps*tau(i)) + eps*s_1(i))*gammadot0*exp(beta3*tau(i) + beta1*T_1(i) + beta2*s_1(i));
    reaction =  A0_hat*exp(E_hat/T_R...
        - E_hat/(T_elastic_centreline_fun(delta,t_p + eps*tau(i),h_flux) + eps*T_1(i)));
%     figure(1)
%     hold on
%     plot(tau(i),plastic_work,'bo',tau(i),reaction,'rx')
    %pause
    if reaction > plastic_work
        tau_star = tau(i);
        t_star = t_p + eps*tau_star;
        T_star = T_elastic_centreline_fun(delta,t_p + eps*tau_star,h_flux) + eps*T_1(i);
        s_star = s_elastic_fun(omega,t_p + eps*tau_star) + eps*s_1(i);
        g_star = gammadot0*eps^(-1/2)*exp(-eps^(-1)*(beta1*(T_p - T_star) + beta2*(s_p - s_star)));
        t_idx = i;
        break
    end
end

[time_reaction3,T_reaction3] = ode15s(@(t,T) A0_hat*eps^(-1/2)*exp(eps^(-1)*T_R*(1 - T_R/T)), [t_star 1],T_star);

save('Approach_III.mat')
%%



% 
% for i = 2:numel(tau)
% figure(2)
% hold on
% plot(tau(i),T_elastic_centreline_fun(delta,t_p + eps*tau(i),h_flux) + eps*T_1(i),'x')
% end