clear;clc;

% specify coupled system parameters
sigma = 10;
rho = 28;
beta = 8/3;

% Model and experienet setup
model_size = 3;               % specify dimension
delta_t = 0.01;        % specify stepsize 
trans_steps = 600;     % specify transient time
time_steps = 16000;     % specify integration time
x0 = [8; 0; 30];         % specify initial conditions

% observation parameters
obs_error_var = 2.0;           % the variance of observation error
obs_grids = [1; 2; 3];                % observation index (location)
% obs_grids = [1; 3];                % observation index (location)
obs_freq_timestep = 8;         % assimilation/observation interval
nobsgrid = length(obs_grids);   % no. of observation for one assimilation cycle
R  =  eye(nobsgrid).*obs_error_var; % Assume diagonal observation error covariance matrix
H = zeros(nobsgrid,model_size);   % initialize observation operator
for i=1:nobsgrid
    H(i,obs_grids(i))=1.;
end

% remove transient time
x = x0;
for i=1:trans_steps
    step_L63;    
end

% redefine initial conditions
x0 = x;

% generate natural run
x = x0;         
truth(:,1) = x;
for time=1:time_steps
    step_L63;
    truth(:,time+1)=x;
end

% generate observations by adding N(0,R) to natural run
nobstime = time_steps/obs_freq_timestep;
noise = randn(nobsgrid,nobstime).*sqrt(obs_error_var);
for i=1:nobstime
    yobs(:,i)=truth(obs_grids,i*obs_freq_timestep+1)+noise(:,i);
end
save('Observation.mat','yobs')
save('Naturalrun.mat','truth')