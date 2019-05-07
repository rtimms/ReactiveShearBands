function [tau,T_1,s_1,gamma_hat] = sb_centreline_harden(h,beta1,beta2,beta3,beta4,rho0,lambda,omega,gammadot0,a,b,Lambda,Lambdak,eps,delta,t_p,T_stop)

% Time to start integration
eta0 = -beta3*t_p/eps + log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2)));

% Initial f,g and timeVec
i = 1;
f(1) = 0;
g(1) = 0;
k(1) = 0;
eta(1) = eta0;

% Set tolerance
tol = 1e-3;                         % Iteration tolerance
d_eta = 1e-2;                       % Initial timestep
stepTol = 1e-2;                     % Set relative change for step size control
max_it = 20;                        % Maximum allowed iterations
converge = 1;                       % Flag to check if fails to converge in max_it iterations

% Set centreline temperature
T_centreline = T_elastic_centreline_fun(delta,t_p,h);


%%%%%%%%    -- Step elastic-plastic solution forward in time --     %%%%%%%
while T_centreline < T_stop && eta(i) < 100
    
    % Increase timestepping index
    i = i+1;
    
    % Loop to repeat timstep if temperature tolerance not met
    proceed = 0;
    
    while proceed == 0
        
        % Update time
        eta(i) = eta(i-1) + d_eta;
        %eta(i)
        
        % Reset tolerance
        diff_fgk = tol*10;
        iterate = 0;
        
        % Use previous result as initial guess
        fNew = f(i-1);
        gNew = g(i-1);
        kNew = k(i-1);

        while diff_fgk > tol && converge == 1
            
            iterate = iterate + 1;
            
            % Break if maximum iterations exceeded
            if iterate > max_it
                f(i) = fNew;
                g(i) = gNew;
                k(i) = kNew;
                sprintf('Maximum number of iterations exceeded!')
                converge = 0;
            end
            
            % Save old values for computation of diff_f and diff_g
            fOld = fNew;
            gOld = gNew;
            kOld = kNew;

            % Calculate updated f and g
            fNew = f_eqn_harden([f(1:i-1) fNew],[g(1:i-1) gNew],[k(1:i-1) kNew],eta);
            gNew = gOld -...
                (gOld - g_eqn_harden([f(1:i-1) fNew],[g(1:i-1) gNew],[k(1:i-1) kNew],Lambda,eta(1:i)))...
                /(1 + g_eqn_harden([f(1:i-1) fNew],[g(1:i-1) gNew],[k(1:i-1) kNew],Lambda,eta(1:i)));
            kNew = k_eqn_harden([f(1:i-1) fNew],[g(1:i-1) gNew],[k(1:i-1) kNew],Lambdak,eta);
            
            % Calculate change in iterated solution
            diff_fgk = max([ norm(fNew-fOld,inf); norm(gNew-gOld,inf); norm(kNew-kOld,inf)]);
            
        end
        
        % Check if time step needs reducing and update results
        fChange = (fNew - f(i-1));
        gChange = (gNew - g(i-1));
        kChange = (kNew - kOld);
        if any(fChange > stepTol) || any(gChange > stepTol) || any(kChange > stepTol)
            % Reduce timestep
            d_eta = d_eta/2;
        else
            % Update results
            proceed = 1;
            f(i) = fNew;
            g(i) = gNew;
            k(i) = kNew;
        end
        
    end
    
    T_centreline = T_elastic_centreline_fun(delta,t_p +...
        eps*(eta(i) - log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2))))/beta3,h)...
        + eps*fNew/beta1;
    
end

% Convert back to tau timescale
tau = (eta - log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2))))/beta3;
T_1 = f/beta1;
s_1 = -g/beta2;
gamma_hat = k/beta4;

end

%% 

function feq = f_eqn_harden(f,g,k,time)
feq = integral(@(t) f_nest_harden(t,f,g,k,time),time(1),time(end));
end

function expF = f_nest_harden(t,f,g,k,time)
ff = interp1(time,f,t);
gg = interp1(time,g,t);
kk = interp1(time,k,t);
expF = (pi*(time(end)-t)).^(-1/2).*exp(t + ff - gg - kk);
end

function geq = g_eqn_harden(f,g,k,A,time)
 geq = A*exp(time(end) + f(end) - g(end) - k(end));
end

function keq = k_eqn_harden(f,g,k,Ak,time)
keq = Ak*integral(@(t) k_nest_harden(t,f,g,k,time),time(1),time(end));
end

function expK = k_nest_harden(t,f,g,k,time)
ff = interp1(time,f,t);
gg = interp1(time,g,t);
kk = interp1(time,k,t);
expK = exp(t + ff - gg - kk);
end




