function [f,g,eta] = fg_solver(Ap,AR,At,eta_R,fR,eta0,T_stop,h,beta1,beta2,beta3,rho0,lambda,omega,gammadot0,a,b,eps,delta,t_p)

% Initial f,g and timeVec
k = 1;
f(1) = 0;
g(1) = 0;
eta(1) = eta0;

% Set tolerance
tol = 1e-6;                         % Iteration tolerance
d_eta = 1e-2;                       % Initial timestep
stepTol = 0.01;                     % Set relative change for step size control

% Set centreline temperature
T_centreline = T_elastic_centreline_fun(delta,t_p,h);

while T_centreline < T_stop
    
    % Increase timestepping index
    k = k+1;
    
    % Loop to repeat timstep if temperature tolerance not met
    proceed = 0;
    
    while proceed == 0
        
        % Update time
        eta(k) = eta(k-1) + d_eta;
        eta(k);
        
        % Reset tolerance
        diff_f = tol*10;
        diff_g = tol*10;
        iterate = 0;
        
        % Use previous result as initial guess
        fNew = f(k-1);
        gNew = g(k-1);
        
        while max(diff_f,diff_g)  > tol;
            
            iterate = iterate + 1;
            
            % Break if maximum iterations exceeded
            if iterate > 10000
                f(k) = fNew;
                g(k) = gNew;
                sprintf('Maximum number of iterations exceeded!')
                return
            end
            
            % Save old values for computation of diff_f and diff_g
            fOld = fNew;
            gOld = gNew;
            
            % Calculate updated f and g
            fNew = fR_eqn([f(1:k-1) fNew],[g(1:k-1) gNew],AR,At,eta_R,fR,eta);
            gNew = gOld -...
                (gOld - g_eqn(fOld,gOld,Ap,eta(end)))/(1 + g_eqn(fOld,gOld,Ap,eta(end)));
            
            % Calculate change in iterated solution
            diff_f = norm(fNew-fOld);
            diff_g = norm(gNew-gOld);
            
        end
        
        % Check if time step needs reducing and update results
        fChange = (fNew - f(k-1));
        gChange = (gNew - g(k-1));
        if any(fChange > stepTol) || any(gChange > stepTol)
            
            % Reduce timestep
            d_eta = d_eta/2;
        else
            % Update results
            proceed = 1;
            f(k) = fNew;
            g(k) = gNew;
        end
        
    end
    T_centreline = T_elastic_centreline_fun(delta,t_p +...
        eps*(eta(k) - log((lambda*gammadot0*(1+omega*t_p))/(b*beta3^(1/2))))/beta3,h)...
        + eps*fNew/beta1
end

end

