%--------------------------------------------------------------
%     Lorenz 3-variable system with data assimilation
%     DA scheme: EnKF with perturbed observations
%     created by DA4fun, 2023
%--------------------------------------------------------------
% This is a simple example of EnKF in Matlab with perturbed observations and the Lorenz 63 model with the 4th order Runge-Kutta method
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
ensemble_size = 40;

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


% load observation data
load('Observation.mat')

% load natural run
load('Naturalrun.mat')
truth=truth(:,obs_freq_timestep+1:obs_freq_timestep:end);

% remove transient time
x = x0;
for i=1:trans_steps
    step_L63;    
end
x0 = x;                
dx0 = [5.; 5.; 5.];             % initial perturbation
%% start DA cycle
ics = randn(ensemble_size,1);
analy = zeros(ensemble_size,model_size,time_steps+1);
for iens=1:ensemble_size
analy(iens,:,1) = x0+dx0+ics(iens,1);
end
Cystart = squeeze(analy(:,:,1));
ncycles = time_steps/obs_freq_timestep;
times = 1;

for iassim = 1:ncycles
    % advance until we have observations
    obsstep = (iassim-1)*obs_freq_timestep + 1;
    step_end = iassim*obs_freq_timestep;          
    for istep = obsstep:step_end
        for iens = 1:ensemble_size
           x = Cystart(iens,:)';
           step_L63;
           Cystart(iens,:) = x';
        end
     % store zens
            analy(:,:,times+1) = Cystart;
            times = times+1;
    end 
    time(iassim)=delta_t*obs_freq_timestep*iassim;
    % perturbed observations, assimilate them by EnKF
    obs = kron(yobs(:,iassim),ones(1,ensemble_size))+sqrt(obs_error_var).*randn(nobsgrid,ensemble_size);
    Xb = analy(:,:,times); % background: N x n
    prior(:,iassim) = mean(Xb); % prior state before DA
    % EnKF with perturbed observations
    Hx = H*Xb';
    Hx_prime = Hx-mean(Hx,2); % n x N
    Xb_prime = Xb-mean(Xb); % N x n
    BHT = Xb_prime'*Hx_prime'/(ensemble_size-1);
    HBHT=Hx_prime*Hx_prime'/(ensemble_size-1);
    KFgain=BHT/(HBHT+R);
    for iens = 1:ensemble_size     
        Xa(iens,:) = Xb(iens,:)+(KFgain*(obs(:,iens)-Hx(:,iens)))';
    end
    poste(:,iassim) = mean(Xa);  % posterior state after DA   
    Cystart = Xa;       
       
    %calculate the root mean square error
    erra(iassim) =sqrt( ( (poste(1,iassim)-truth(1,iassim)).^2+ (poste(2,iassim)-truth(2,iassim)).^2 +(poste(3,iassim)-truth(3,iassim)).^2)/3.);
    errf(iassim) =sqrt( ( (prior(1,iassim)-truth(1,iassim)).^2+ (prior(2,iassim)-truth(2,iassim)).^2 +(prior(3,iassim)-truth(3,iassim)).^2)/3.);
end
% mean RMS analysis/background error
rmse_a=mean(erra(500:ncycles));
rmse_f=mean(errf(500:ncycles));
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
print -f1 -dpng ../output/output2.png;
