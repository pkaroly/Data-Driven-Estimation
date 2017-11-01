%% ************ ESTIMATION OF NMM *****************
% This script runs the code used to generate estimation results in
% Ahmadizadeh et al. 2016

% AUTHORS
% Dean Freestone, Philippa Karoly 2016

% Further details on the estimation method can be found in 
% the following references:

% [1] Freestone, Dean R., Philippa J. Karoly, Dragan Neši?, 
% Parham Aram, Mark J. Cook, and David B. Grayden. 
% "Estimation of effective connectivity via data-driven neural modeling." 
% Frontiers in neuroscience 8 (2014): 383.

% This code is licensed under GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007

%% initalize parameters
anneal_on = 1;
kappa_0 = 1000;


% LOAD SEIZURE
load('SeizureData');

% scale data
new_min = -10;
new_max = 10;
load('dataRange')
data = Seizure';
dataRange = maxEEG - minEEG;   % based on meta-data from patient seizures
m = mean(data,2);
% range scaled according to seizure range, and zero mean
data(1,:) = (new_max - new_min) * (data(1,:) - m(1)) ./ dataRange(1);
data(2,:) = (new_max - new_min) * (data(2,:) - m(2)) ./ dataRange(2);

% INTERPOLATING DATA
data1 = interp(data(1,:),3);
data2 = interp(data(2,:),3);
data = [data1 ; data2];
clear('data1','data2');

mkdir('Seizure_Estimates');


%% START ESTIMATION
case_ind = 1:3;
% loop through different model cases
for iCase = case_ind
    
    % load the model parmas
    load(['/Model' num2str(iCase) '.mat']);
    
    TimeOfSim = 5;                           % time of simulation (s)
    N_samples = round(TimeOfSim/dt);         % number of time points to simulate
    randomOn = 0;                            % forward sim will be non/random
    
    % generate some data
    generateData;
    
    % truncate the alphas
    xi = xi(1:N_est_states,:);
    
    %%
    % INIT MEAN AND COVARIANCE
    xi_hat_init = mean(xi(:,N_samples/2:end),2);                 % to ignore inital transient take the mean of the second half of the test data
    %     xi_hat_init(21) = -1;
    P_hat_init = 10*blkdiag(cov(xi(1:10,N_samples/2:end)'),cov(xi(11:end,N_samples/2:end)'));
    
    % open up the error in external input
    P_hat_init(2*N_syn+1:2*N_syn+N_inputs,2*N_syn+1:2*N_syn+N_inputs) = diag([1,1e-4*ones(1,N_inputs-1)]); 
        
    xi_hat = zeros(N_est_states,N_samples);                 % initialize for speed
    P_hat = zeros(N_est_states,N_est_states,N_samples);     % initialize for speed
    P_diag = zeros(N_est_states,N_samples);                 % initialize for speed
   
    xi_hat(:,1) = xi_hat_init;
    P_hat(1:N_est_states,1:N_est_states,1) = P_hat_init;
    
    tic
    N_samples = length(data);
    t_end_anneal = N_samples/20;
    
    for t=2:N_samples       % N_samples
        
        % augment x and P with alphas
        xi_0p = [xi_hat(:,t-1) ; Alpha'];
        P_0p = blkdiag(P_hat(:,:,t-1),zeros(N_syn,N_syn)); % no uncertainty on alpha
        
        % predict
        %
        [xi_1m, P_1m] = prop_mean_and_cov(N_syn,N_states,N_inputs,A,B,C,P_0p,xi_0p,varsigma,v0,Q);
        
        if (t<=t_end_anneal) && anneal_on
            kappa = kappa_0^((t_end_anneal-t)/(t_end_anneal-1));
        else
            kappa = 1;
        end
        
        % truncate alphas
        xi_1m = xi_1m(1:N_est_states);
        P_1m = P_1m(1:N_est_states,1:N_est_states);
        
        % Kalman update eqn
        %
        K = P_1m*H'/(H*P_1m*H' + kappa*R*eye(2));
 
        xi_uc = xi_1m + K*(data(:,t) - H*xi_1m);
        P_uc = (eye(N_est_states) - K*H)*P_1m;
        
        % if there are constraints on mean and cov apply here
        % (NA)
        xc = xi_uc;
        Pc = P_uc;
        
        % save
        xi_hat(:,t) = xc;
        P_hat(:,:,t) = Pc;
        P_diag(:,t) = diag(Pc);
        
    end

    save(['Seizure_Estimates/Case' num2str(iCase)],'xi_hat','P_diag','-v7.3');
    
    toc

    
end