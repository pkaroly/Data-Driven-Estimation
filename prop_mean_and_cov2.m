function [xi_1m, P_1m] = prop_mean_and_cov2(N_syn,N_states,N_inputs,A,B,C,P_0p,xi_0p,varsigma,v0,Q)

v_indexes = 1:2:2*N_syn+2;              % extra plus 2 for v_in the input!
z_indexes = 2:2:2*N_syn;
alpha_indexes = 2*N_syn+N_inputs+1:N_states;

w_fast = [0.1713244923791705 0.3607615730481384  0.4679139345726904 0.1713244923791705  0.3607615730481384 0.4679139345726904];
y_fast = [0.033765242898424 0.169395306766868 0.380690406958402 0.966234757101576 0.830604693233132 0.619309593041599];

CPC = C(z_indexes,:)*P_0p*C(z_indexes,:)';
dCPB = diag(C(z_indexes,:)*P_0p*B(z_indexes,:)');
Bxi = B(z_indexes,alpha_indexes)*xi_0p(alpha_indexes);
AP = A*P_0p;

%% analytic mean
%
gamma = 1./sqrt(2*(diag(CPC) + varsigma^2));
beta = (C(z_indexes,v_indexes)*xi_0p(v_indexes) - v0).*gamma;       % <- v_0
Xi = (erf(beta) + 1)/2;
Upsilon = exp(-(beta.^2)).*gamma./sqrt(pi);
psi = Bxi.*Xi + dCPB.*Upsilon;                              % E[ Bxi_t o g(Cxi_t) ]

xi_1m = A*xi_0p;                                            % m is for minus (prior (or previous a posteriori)), p for plus (a priori)
xi_1m(2:2:2*N_syn) = xi_1m(2:2:2*N_syn) + psi;

%% exact covariance
%
% cov part 1 (Phi)
%
q2 = Upsilon.*(Bxi - dCPB.*beta.*gamma.^2 *2);
Phi = ones(N_states,1)*q2'.*(AP*C(z_indexes,:)') + ones(N_states,1)*Xi'.*(AP*B(z_indexes,:)');

% note: A term A*xi_t * E[ Bxi_t o g(Cxi_t) ]' = A*xi_0p*psi' is cancelled
% from here from the full covariance expansion.

% cov part 2 (Omega)
%
% NOTE: we can reduce the dimensionality (only compute upper triangle part) of the matrices we turn to vectors to increase speed in larger networks.
%
CPCgammagamma = asin(CPC.*(gamma*gamma')*2);
CPCgammagamma = CPCgammagamma(:);                                           % change to a vector
CPCgammagamma_y = CPCgammagamma*y_fast;                                     % a matrix of row vectors

betabeta = (beta*beta')*2;
betabeta_mat = betabeta(:)*ones(1,length(w_fast));                          % change to a vector and repmat

beta2mat = (beta.^2)*ones(1,N_syn);                                         % sq and repmat
beta2matT = beta2mat';                                                      % the transpose allow for permutations when we sum below
beta2_plus_beta2T = (beta2mat(:) + beta2matT(:))*ones(1,length(w_fast));    % change to a vectors, add, and repmat

% put it together
%
Psi = reshape(sum(CPCgammagamma*w_fast.*exp(-(beta2_plus_beta2T - betabeta_mat.*sin(CPCgammagamma_y))./cos(CPCgammagamma_y).^2),2),4,4)/(4*pi);
%                 ~~~~~~~~~~~~~~~~~~~~  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%             ~~~
Omega = ((Xi*Xi') + Psi).*(Bxi*Bxi' + P_0p(alpha_indexes,alpha_indexes));
%        ~~~~~~~~~~~~~~    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         E[g^2(Cxi)]           (alpha^2 + sigma_alpha^2)

%% now allowing for correlations between params and states
% ~~~~~~~~
% v = sqrt(2*diag(d) + varsigma_sq);            % a strange + + was here: v= sqrt(2*diag(d) + + varsigma_sq); removed a plus
% phi = (erf(l.*varsigma./v) + 1)/2;
% z = e./f .* exp(-c.^2./(2*f)) .* ...
%     ((varsigma*e./(pi*v)) .* exp(1 - varsigma^3*c.^2./(2*sqrt(f).*v)) + phi.* (4*b.*f - 2*c.*e)./sqrt(2*pi*f));
% ~~~~~~~~
%%
% here we construct the cov mat.
%
P_1m = AP*A' + Q;
P_1m(z_indexes,z_indexes) = P_1m(z_indexes,z_indexes)  + Omega - psi*psi';
P_1m(:,z_indexes) = P_1m(:,z_indexes) + Phi;
P_1m(z_indexes,:) = P_1m(z_indexes,:) + Phi';


end

