%--------------------------------------------------------------
%     Lorenz 3-variable system with data assimilation
%     DA scheme: 3D-Var component
%     created by DA4fun, 2023
%--------------------------------------------------------------
% This is a simple example of 3DVar in MATLAB using steepest descent and the Lorenz 63 model with the 4th order Runge-Kutta method
clear;clc;
addpath ../advance_model 
addpath ../utils % plot tools
addpath ../data % observation file, truth file


% specify coupled system parameters
sigma = 10;
rho = 28;
beta = 8/3;

% Model and experienet setup
model_size = 3;               % specify dimension
delta_t = 0.01;        % specify stepsize 
trans_steps = 600;     % specify transient time
time_steps = 16000;     % specify integration time
x0 = [8; 0; 30];       % specify initial conditions

% observation parameters
obs_error_var = 2.0;           % the variance of observation error
obs_grids = [1; 2; 3];                % observation index (location)
% obs_grids = [1; 3];                % observation index (location) if you
% want to change obs_grids, you need to regenerate observations
obs_freq_timestep = 8;         % assimilation/observation interval
nobsgrid = length(obs_grids);   % no. of observation for one assimilation cycle
R  =  eye(nobsgrid).*obs_error_var; % Assume diagonal observation error covariance matrix
H = zeros(nobsgrid,model_size);   % initialize observation operator
for i=1:nobsgrid
    H(i,obs_grids(i))=1.;
end


% load observation data
load('Observation.mat')

% load natural run
load('Naturalrun.mat')
truth=truth(:,obs_freq_timestep+1:obs_freq_timestep:end);

% load the background error covariance matrix
load('Static_B.mat')


% remove transient time
x = x0;
for i=1:trans_steps
    step_L63;    
end
x0 = x;                
dx0 = [5.; 5.; 5.];             % initial perturbation

%% control run
contr(:,1) = x0+dx0;            % control run initial conditions
Cystart = contr(:,1);
ncycles = time_steps/obs_freq_timestep;
times = 1;
for iassim = 1:ncycles
    % advance until we have observations
    obsstep = (iassim-1)*obs_freq_timestep + 1;
    step_end = iassim*obs_freq_timestep;          
    for istep= obsstep:step_end        
           x = Cystart;
           step_L63;
           Cystart = x;
     % store zens
            contr(:,times+1) = Cystart;
            times = times+1;
    end 
    stateraw(:,iassim) = contr(:,times);
    errraw(iassim) =sqrt( ( (stateraw(1,iassim)-truth(1,iassim)).^2+ (stateraw(2,iassim)-truth(2,iassim)).^2 +(stateraw(3,iassim)-truth(3,iassim)).^2)/3.);
end


%% start DA cycle
analy = zeros(model_size,time_steps+1);
analy(:,1) = x0+dx0;
Cystart = analy(:,1);
ncycles = time_steps/obs_freq_timestep;
times = 1;

% define the iteration parameters
eps = 1e-6;
max_iter = 200;
alpha = 0.1;

for iassim = 1:ncycles
    % advance until we have observations
    obsstep = (iassim-1)*obs_freq_timestep + 1;
    step_end = iassim*obs_freq_timestep;          
    for istep= obsstep:step_end        
           x = Cystart;
           step_L63;
           Cystart = x;
     % store zens
            analy(:,times+1) = Cystart;
            times = times+1;
    end 
    time(iassim)=delta_t*obs_freq_timestep*iassim;
    % obtain observations, assimilate them by 3Dvar
    obs = yobs(:,iassim);
    Xb = analy(:,times); % background: n x 1
    X = Xb; % first guess of X: n x 1  
    prior(:,iassim) = X; % prior state before DA 
    % Loop until convergence or maximum number of iterations reached
    for iter = 1:max_iter
    
         % Calculate the gradient of the cost function
         grad_J = inv(B)*(X-Xb) + H'*inv(R)*(H * X-obs);
    
         % Update the analysis
         X = X - alpha * grad_J;
    
         % Check for convergence
         if norm(grad_J) < eps
             break;
         end
    
    end
    poste(:,iassim) = X;  % posterior state after DA   
    Cystart = X;       
       
    %calculate the root mean square error
    erra(iassim) =sqrt( ( (poste(1,iassim)-truth(1,iassim)).^2+ (poste(2,iassim)-truth(2,iassim)).^2 +(poste(3,iassim)-truth(3,iassim)).^2)/3.);
    errf(iassim) =sqrt( ( (prior(1,iassim)-truth(1,iassim)).^2+ (prior(2,iassim)-truth(2,iassim)).^2 +(prior(3,iassim)-truth(3,iassim)).^2)/3.);
end
% mean RMS analysis/background error
rmse_a=mean(erra(500:ncycles));
rmse_f=mean(errf(500:ncycles));
rmse_raw=mean(errraw(500:ncycles));
fprintf(1,'No assimilation, Mean RMSE= %g \n',rmse_raw)
fprintf(1,'Obseveraion error= %g \n',sqrt(obs_error_var))
fprintf(1,'Number of DA cycles= %g \n',ncycles)
fprintf(1,'Mean background RMSE= %g \n',rmse_f)
fprintf(1,'Mean analysis RMSE= %g \n',rmse_a)

ind=[1:500];
% plot the error
DAerrplt(time(ind),erra(ind),rmse_a,sqrt(obs_error_var),obs_grids);
print -f1 -dpng ../output/output1.png;
% plot the states vs. obs
stateplt(time(ind),yobs(1,ind),prior(1,ind),poste(1,ind),'X')
print -f2 -dpng ../output/output2.png;
