function [f,g,eta,converge] = fg_reactive_solver(f0,g0,eta0,eta_Final,eta_R,f_R,Ap,AR,At,fMax,gMax,tol,d_eta,stepTol,it_Max)

% Set initial values
k = 1;
f(1) = f0;
g(1) = g0;
eta(1) = eta0;

% Step forward in time

while eta < eta_Final
    
    % Increase timestepping index
    k = k+1;
    
    % Loop to repeat timstep if temperature tolerance not met
    proceed = 0;
    
    while proceed == 0
        
        % Update time
        eta(k) = eta(k-1) + d_eta;
        
        % Reset tolerance
        diff_f = tol*10;
        diff_g = tol*10;
        iterate = 0;
        
        % Use previous result as initial guess
        fNew = f(k-1);
        gNew = g(k-1);
        
        while max(diff_f,diff_g)  > tol && iterate < it_Max;
            
            iterate = iterate + 1;
            
            % Save old values for computation of diff_f and diff_g
            fOld = fNew;
            gOld = gNew;
            
            % Calculate updated f and g
            fNew = fR_eqn([f(1:k-1) fNew],[g(1:k-1) gNew],AR,At,eta_R,f_R,eta);
            gRHS = g_eqn(fOld,gOld,Ap,eta(end));
            gNew = gOld -...
                (gOld - gRHS)/(1 + gRHS);
            
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
    
    % Break if maximum iterations exceeded
    if iterate == it_Max
        sprintf('Maximum number of iterations exceeded in fg_reactive_solver! stepTol exceeded by %f',abs(fChange-stepTol))
        converge = 0;
        break
    end
    
    % Break if d_eta too small
    if d_eta < 1e-6
        sprintf('Minimum timestep!')
        converge = 0;
        break
    end
        
    % Break if fMax or gMax
    if fNew > fMax || gNew > gMax
        sprintf('fMax or gMax exceeded!')
        converge = 0;
        break
    end
    

        
end

% Converged 
converge = 1;

end

