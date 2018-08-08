%% Estimation
% runs the state/parameter estimation and plots results for a single
% channel at a time

%% 
% Dean Freestone, Philippa Karoly 2016
% This code is licensed under the MIT License 2018


%%
clear
close all
clc

load('../data/Seizure_1.mat');  % change this path to load alternative data
addpath(genpath('../src/'));

input = 300;
input_offset = [];
% generate some data
%
time = 5;
Fs = 0.4e3;
[A,B,C,N_states,N_syn,N_inputs,N_samples,xi, ...
    v0,varsigma,Q,R,H,y] = set_params(input,input_offset,time,Fs);

if ~isempty(input_offset)
    % reset the offset
    my = mean(y);
    input_offset = -my/0.0325;   % if in mV
    [A,B,C,N_states,N_syn,N_inputs,N_samples,xi, ...
        v0,varsigma,Q,R,H,y] = set_params(input,input_offset);
end

xi_hat_init = mean(xi(:,N_samples/2:end),2);                                % to ignore inital transient take the mean of the second half of the test data
P_hat_init = 10*cov(xi(:,N_samples/2:end)');
P_hat_init(2*N_syn+1:end,2*N_syn+1:end) = eye(N_syn+N_inputs)*10e-2;               % open up the error in parameters

% set inital conditions for the KF
xi_hat = zeros(N_states,N_samples);                 % initialize for speed
P_hat = zeros(N_states,N_states,N_samples);         % initialize for speed
P_diag = zeros(N_states,N_samples);                 % initialize for speed

xi_hat(:,1) = xi_hat_init;                                % to ignore inital transient take the mean of the second half of the test data
P_hat(:,:,1) = P_hat_init;

anneal_on = 1;
kappa_0 = 10000;
t_end_anneal = N_samples/20;

% loop through 16 chans
for iCh = 1:16
    
    fprintf('Channel %02d ...',iCh);
    
    % get one channel at a time
    % NB - portal data is inverted. we need to scale it to some
    % 'reasonable' range for the model, but still capture amplitude
    % differences bw seizures
    y = -0.5*Seizure(:,iCh);
    N_samples = length(y);
    
    for t=2:N_samples       % N_samples
        
        xi_0p = squeeze(xi_hat(:,t-1));
        P_0p = squeeze(P_hat(:,:,t-1));
        
        % predict
        %
        [xi_1m, P_1m] = prop_mean_and_cov(N_syn,N_states,N_inputs,A,B,C,P_0p,xi_0p,varsigma,v0,Q);
        
        if (t<=t_end_anneal) && anneal_on
            kappa = kappa_0^((t_end_anneal-t)/(t_end_anneal-1));
        else
            kappa = 1;
        end
        
        K = P_1m*H'/(H*P_1m*H' + kappa*R);
        
        % correct
        %
        xi_hat(:,t) = xi_1m + K*(y(t) - H*xi_1m);
        P_hat(:,:,t) = (eye(N_states) - K*H)*P_1m;
        P_diag(:,t) = diag(squeeze(P_hat(:,:,t)));
        
        if t > 2
            fprintf('\b\b\b\b');
        end
        fprintf('%03d%%', round(100*t/N_samples));
    end
    
    close all
    figure('name','membrane potential estimates' ,'units','normalized','position',[0 0 1 1] )
    subplot(411),plot(xi_hat(1,:))
    subplot(412),plot(xi_hat(3,:))
    subplot(413),plot(xi_hat(5,:))
    subplot(414),plot(xi_hat(7,:))
    
    figure('name','parameter estimates' ,'units','normalized','position',[0 0 1 1] )
    subplot(511),plot(xi_hat(9,:))
    title('Input');
    subplot(512),plot(xi_hat(10,:))
    title('Inhibitory -> Pyramidal');
    subplot(513),plot(xi_hat(11,:))
    title('Pyramidal -> Inhibitory');
    subplot(514),plot(xi_hat(12,:))
     title('Pyramidal -> Excitatory');
    subplot(515),plot(xi_hat(13,:))
    title('Excitatory -> Pyramidal');
    drawnow;
end