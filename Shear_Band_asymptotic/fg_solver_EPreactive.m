function [f,g,eta] = fg_solver_EPreactive(h,beta1,beta2,beta3,rho0,lambda,omega,gammadot0,a,b,Lambda_p,Lambda_R,Lambda_t,eps,delta,t_p,T_stop,fR,eta0,eta_R)

% Initial f,g and timeVec
k = 1;
f(1) = 0;
g(1) = 0;
eta(1) = eta0;

% Set tolerance
tol = 1e-3;                         % Iteration tolerance
d_eta = 1e-2;                       % Initial timestep
stepTol = 1e-2;                     % Set relative change for step size control
max_it = 20;                        % Maximum allowed iterations
converge = 1;                       % Flag to check if fails to converge in max_it iterations

% Set centreline temperature
T_centreline = T_elastic_centreline_fun(delta,t_p,h);

while T_centreline < T_stop && converge == 1 && eta(k) < 100
    
    % Increase timestepping index
    k = k+1;
    
    % Loop to repeat timstep if temperature tolerance not met
    proceed = 0;
    
    while proceed == 0
        
        % Update time
        eta(k) = eta(k-1) + d_eta;
        eta(k)
        
        % Reset tolerance
        diff_f = tol*10;
        diff_g = tol*10;
        iterate = 0;
        
        % Use previous result as initial guess
        fNew = f(k-1);
        gNew = g(k-1);
        
        
        while max(diff_f,diff_g)  > tol && converge == 1;
            
            iterate = iterate + 1;
            
            % Break if maximum iterations exceeded
            if iterate > max_it
                f(k) = fNew;
                g(k) = gNew;
                sprintf('Maximum number of iterations exceeded in fg_solver_EPreactive !')
                converge = 0;
            end
            
            % Save old values for computation of diff_f and diff_g
            fOld = fNew;
            gOld = gNew;
            
            % Calculate updated f and g
            fNew = f_EPreactive([f(1:k-1) fNew],[g(1:k-1) gNew],Lambda_R,Lambda_t,eta_R,fR,eta);
            gNew = gOld...
                - (gOld - g_EPreactive(fOld,gOld,Lambda_p,eta(end)))...
                /(1 + g_EPreactive(fOld,gOld,Lambda_p,eta(end)));
            
            % Calculate change in iterated solution
            diff_f = norm(fNew-fOld,inf);
            diff_g = norm(gNew-gOld,inf);
            
        end
        
        % Check if time step needs reducing and update results
        fChange = (fNew - f(k-1));
        gChange = (gNew - g(k-1));
        
        if any(fChange > stepTol) || any(gChange > stepTol)
            
            % Reduce timestep
            d_eta = d_eta/2;
            
            % Break if d_eta too small
            if d_eta < 1e-6
                sprintf('Minimum timestep!')
                converge = 0;
                f(k) = fNew;
                g(k) = gNew;
                break
            end
            
        else
            % Update results
            proceed = 1;
            f(k) = fNew;
            g(k) = gNew;
        end
        
    end
    
    T_centreline = T_elastic_centreline_fun(delta,t_p +...
        eps*(eta(k) - log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2))))/beta3,h)...
        + eps*fNew/beta1;
    
end

end

%% Functions for Newton iteration

function f = f_EPreactive(f,g,AR,At,tR,fR,time)
f = integral(@(t) f_nest(t,f,g,AR,At,tR,fR,time),time(1),time(end));
end

function g = g_EPreactive(f,g,Ap,time)
g = Ap*exp(time + f - g);
end

function df = dfdg_EPreactive(f,g,time)
df = integral(@(t) f_nest(t,f,g,0,0,0,0,time),time(1),time(end));
end

function expF = f_nest(t,f,g,AR,At,tR,fR,time)
ff = interp1(time,f,t);
gg = interp1(time,g,t);
expF = (pi*(time(end)-t)).^(-1/2)...
    .*(exp(t + ff - gg) + AR*exp(At*(t-tR) + ff - fR));
end
