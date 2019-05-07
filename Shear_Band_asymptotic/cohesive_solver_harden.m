function [T_centreline_max,t_end] = cohesive_solver_harden(rho_hat,S,lambda,omega,gammadot0,eps,beta1,beta2,beta4,T_p,s_p,t_p,q_flux,N,y,dy,dt,FILENAME,t_stop)

% Function to solve elastic-plastic shear band onset using the method of
% characteristics for the mechanical equations and second order finite
% difference for the thermal part of the problem

% Initial values
time = 0;
velNew = omega*y;
stressNew = ones(size(y));
gdotNew = zeros(size(y));
gNew = zeros(size(y));
TNew = ones(size(y));

% Save to file
filename = sprintf('numerical_%s_decoupled_N%.0f',FILENAME,N);
save(sprintf('output/%s',filename))
time_file = fopen(sprintf('output/time_%s.bin',filename),'w');
fwrite(time_file,time,'double');
fclose(time_file);
out = [velNew(1:round(N/10)) stressNew(1:round(N/10)) gdotNew(1:round(N/10)) gNew(1:round(N/10)) TNew(1:round(N/10))];
out_file = fopen(sprintf('output/out_%s.bin',filename),'w');
fwrite(out_file,out','double');
fclose(out_file);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Timestepping routine

% Timestep counter
timestep = 0;

while time < t_stop
    
    % Print timestep
    timestep = timestep + 1;
    time = time + dt;
    %[time max(TNew)]

    
    % Rename previous values
    velOld = velNew;
    stressOld = stressNew;
    TOld = TNew;
    gdotOld = gdotNew;
    gOld = gNew;
    
    
    % Update mechanical equations
    [velNew,stressNew,gdotNew,gNew,converge_m] = ...
        mechanical_solver_harden(velOld,stressOld,gdotOld,gOld,TOld,rho_hat,S,omega,gammadot0,eps,beta1,beta2,beta4,T_p,s_p,y,dy,dt,N);
    
    % Save final timestep if failed to converge
    if converge_m == 0
        sprintf('Mechanical solver failed to converge at timestep %0.f',timestep)
        time_file = fopen(sprintf('output/time_%s.bin',filename),'a+');
        fwrite(time_file,time,'double');
        fclose(time_file);
        out = [velNew(1:round(N/10)) stressNew(1:round(N/10)) gdotNew(1:round(N/10)) gNew(1:round(N/10)) TNew(1:round(N/10))];
        out_file = fopen(sprintf('output/out_%s.bin',filename),'a+');
        fwrite(out_file,out','double');
        fclose(out_file);
        
        % Output maximum temperature on centreline
        T_centreline_max = max(TNew);
        t_end = time;
        return
    end
    
    
    % Update temperature
    TNew = thermal_solver(stressOld,gdotOld,TOld,stressNew,gdotNew,lambda,dy,dt,N,time,q_flux);
    
    
    % Break in inf temp
    if any(isnan(TNew)==1)
        sprintf('Temperature NaN at timestep %0.f',timestep)
        time_file = fopen(sprintf('output/time_%s.bin',filename),'a+');
        fwrite(time_file,time,'double');
        fclose(time_file);
        out = [velNew(1:round(N/10)) stressNew(1:round(N/10)) gdotNew(1:round(N/10)) gNew(1:round(N/10)) TNew(1:round(N/10))];
        out_file = fopen(sprintf('output/out_%s.bin',filename),'a+');
        fwrite(out_file,out','double');
        fclose(out_file);
        
        % Output maximum temperature on centreline at previous timestep
        T_centreline_max = max(TOld);
        t_end = time - dt;
        return
    end
    
    % Save new solutions to file
    time_file = fopen(sprintf('output/time_%s.bin',filename),'a+');
    fwrite(time_file,time,'double');
    fclose(time_file);
    out = [velNew(1:round(N/10)) stressNew(1:round(N/10)) gdotNew(1:round(N/10)) gNew(1:round(N/10)) TNew(1:round(N/10))];
    out_file = fopen(sprintf('output/out_%s.bin',filename),'a+');
    fwrite(out_file,out','double');
    fclose(out_file);
    
    % Plot results
    %if mod(timestep,100) == 1
%                 figure(101)
%                 subplot(2,2,1)
%                 plot(velNew,y)
%                 title('Velocity')
%                 subplot(2,2,2)
%                 plot(stressNew,y)
%                 title('Stress')
%                 subplot(2,2,3)
%                 plot(gdotNew,y)
%                 title('Strain-rate')
%                 subplot(2,2,4)
%                 plot(TNew,y)
%                 title('Temperature')
%                 pause(0.001)
    %end
    
    
end

% Output maximum temperature on centreline
T_centreline_max = max(TNew);
t_end = time;

end


