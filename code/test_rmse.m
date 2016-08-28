%%
% Evolution with respect to n.

%%
% Results : 

% 'coherent', t = 3;
% B = 3.735e-02
% 'shrodinger-cat';  t = 3; 
% B = 2.273e-01

% decay exponent of the rho, controls the value of N
r0 = 2;
B0 = .3;
B0 = .8;
B0 = .5;


if not(exist('name'))
name = 'single-photon'; 
name = 'vacuum';
name = 'shrodinger-cat';  t = 3; 
name = 'coherent'; t = 3; % sqrt(9/2); 
name = 'thermal'; t = 1/10;
name = 'thermal'; t = 1/4;
end

namesvg = [name num2str(round(t*10))];

addpath('toolbox/');
rep = 'results/';
if not(exist(rep))
    mkdir(rep);
end

%%
% Helpers.

mynorm = @(x)norm(x(:));
saveeps = @(ext)saveas(gcf, [rep namesvg ext '.eps'], 'epsc');
set_gca = @()set(gca, {'XTick' 'YTick' 'ZTick' 'FontSize'}, {[] [] [] fs});
fs = 20; % font size
lw = 2; % linewidth

%%
% Parameter of the problem

% relationship between n and N
Nval = @(n)floor( (log(n)/(2*B0)).^(2/r0) );
% resolution of the rho to compute the error
N0 = 30;
% Noise factor, eta=1 <=> no noise. In practice, eta in [.8,.9] 
eta = .9;

%%
% Original rho to recover.

rho0 = compute_rho(name, N0, t);
clf; bar3(rho0); axis tight;
colormap jet(256);

%%
% Estimator by soft thresholding.


soft_thresh = @(x,gamma) max(0, 1-gamma./max(1e-10,abs(x))) .* x;
Estimate = @(rho,N,n,epsilon,lambda,Vinf) soft_thresh(rho, ...
    lambda * Vinf * sqrt( 2*log(N*(N+1)/epsilon) / n ) );

if 0 % OLD
[~,Vinf] = perform_rho_deconvolution([],[],eta, N);
soft_thresh = @(x,gamma) max(0, 1-gamma./max(1e-10,abs(x))) .* x;
Estimate = @(rho,n,epsilon,lambda) soft_thresh(rho, ...
    lambda * Vinf * sqrt( 2*log(N*(N+1)/epsilon) / n ) );
end


%%
% Evolution with n of the error. 

nNbr = 40; nNbr = 20;
nList = round(linspace(5000, 10*100000, nNbr));
lambda_list = linspace(.1,1.5,60);
ntrials = 20; ntrials = 2; % number of runs to compute the standard deviation
epsilon = 1;
Err = [];
for i=1:length(nList)
    progressbar(i, length(nList));
    
    % number of samples.
    n = nList(i);
    % number of estimated atoms
    N = Nval(n);
    % need to compute L^inf norm of regressors
    [~,Vinf] = perform_rho_deconvolution([],[],eta, N);
    
    err = [];
    for k=1:ntrials
        % gen samples        
        [X0,Phi] = perform_sampling(name, n, t);
        X = sqrt(eta)*X0 + sqrt((1-eta)/2)*randn(n,1);
        % deconvolve
        rho = perform_rho_deconvolution(X,Phi,eta, N);
        % estimate    
        for j=1:length(lambda_list)
            rho1 = Estimate(rho,N,n,epsilon,lambda_list(j),Vinf);
            Rho1 = zeros(N0); Rho1(1:N,1:N) = rho1;
            Err(k,i,j) = mynorm(rho0-Rho1)^2;        
        end
    end
end

RMSE = squeeze( sqrt( mean(Err) ) ) / mynorm(rho0);
RMSE_std = squeeze(  std( sqrt(Err) ) ) / mynorm(rho0);

%%
% Log-log fit of evolution with respect to n.
% 	error^2 ~ n^{ -B/(4*gamma+B) }

gamma = (1-eta)/(4*eta);
sel = 3:length(nList);
lambda = 1; [~,i] = min( abs(lambda_list-lambda) );
clf;
plot(log10(nList(sel)),log10(RMSE(sel,i)), '.-');
axis tight;

%%
% Linear regression
% log10(error) = -a*log10(n)+b
% error = c * n^{-a},  c = 10^b
% 2*a = B/(4*gamma+B)

P = polyfit(log10(nList(sel)),log10(RMSE(sel,i))',1);
a = -P(1); c = 10^P(2);
B = 8*a*gamma / ( 1 - 2*a );
fprintf('B = %.3e\n', B);


%%
% Display evolution with respect to n.
lw = 2;

clf; hold on;
lambda = 1; [~,i] = min( abs(lambda_list-lambda) );
h = shadedErrorBar(nList, RMSE(:,i), 3*RMSE_std(:,i), 'b');
set(h.mainLine, 'LineWidth', lw);
%
lambda = .8; [~,i] = min( abs(lambda_list-lambda) );
plot(nList, RMSE(:,i),'r--', 'LineWidth', lw);
axis([min(nList) max(nList) 0 1]);
%
lambda = .5; [~,i] = min( abs(lambda_list-lambda) );
plot(nList, RMSE(:,i),'g--', 'LineWidth', lw);
% display fit
% plot(nList, c*nList.^(-a),'k:', 'LineWidth', lw);
% legend('\lambda=1', '\lambda=0.8', '\lambda=0.5');
axis([min(nList) max(nList) 0 1]);
% xlabel('n'); ylabel('E(|\rho^\eta - \rho|^2)');
box on;
set(gca, 'FontSize', fs);
saveeps( ['-eta' num2str(round(eta*10)) '-rmse'] );

if 0
% plot for a fixed n and varying lambda
n = 20000; [~,i] = min( abs(nList-n) );
clf;
plot(lambda_list, RMSE(i,:), 'LineWidth', lw);
xlabel('\lambda'); ylabel('E(|\rho^\eta - \rho|^2)');
axis([min(lambda_list) max(lambda_list) 0 1]);
% saveeps( ['-eta' num2str(eta) '-evolLambda-n' num2str(round(n/1e3)) 'k'] );
end

% save workspace
save([rep namesvg '.mat']);