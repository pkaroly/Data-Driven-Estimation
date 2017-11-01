% GENERATE DATA
% ~~~~~~~~~~~~~~~~~~~~~~~
% set up the forward simulation
%

% define initial conditions
%

x = zeros(1,2*N_syn);


xi = zeros(N_states,N_samples);
xi(:,1) = [x ex_input input_offset Alpha]';

rng(1)
w = mvnrnd(zeros(1,N_states),Q,N_samples)';

% run the model forward
%
for n=1:N_samples-1
    
    phi = g(C*xi(:,n),v0,varsigma);
    xi(:,n+1) = A*xi(:,n) + B*xi(:,n).*phi + w(:,n);
    
end

if ~randomOn
    rng(1)
end
v = sqrt(R)*randn(2,N_samples);
y = H*xi(1:N_est_states,:) + v;                           % simulated EEG signal