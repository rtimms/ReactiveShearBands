function fg_solve_test
close all
clear all
clc

% Parameters
Ap = 1;
AR = 0;
At = 0;
eta_R = 0;
f_R = 0;

% Set initial values
k = 1;
f(1) = 0;
g(1) = 0;
eta(1) = -30;
eta_Final = 20;
fMax = 15;
gMax = 15;

% Set tolerance
tol = 1e-3;                        
d_eta = 1e-1;                      
stepTol = 1e-2;                    
it_Max = 10;                     


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
        break
    end
    
    % Break if d_eta too small
    if d_eta < 1e-6
        sprintf('Minimum timestep!')
        break
    end
        
    % Break if fMax or gMax
    if fNew > fMax || gNew > gMax
        sprintf('fMax exceeded!')
        break
    end
    

        
end

figure(1)
plot(eta,f,'r',eta,g,'b',eta,exp(eta),'g:')

end