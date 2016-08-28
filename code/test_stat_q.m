%%
% Display of the rho matrices and their estimation for a few n.

name = 'single-photon'; 
name = 'vacuum';
name = 'thermal'; t = 1/10;

name = 'shrodinger-cat';  t = 3; 
name = 'thermal'; t = 1/4;
name = 'coherent'; t = 3; % sqrt(9/2); 

% decay exponent of the rho, controls the value of N
r0 = 2;
B0 = .5;
B0 = .3;
B0 = .8;

namesvg = [name num2str(round(t*10))];

addpath('toolbox/');

rep = ['results/' name '/b0-' num2str(round(B0*10)) '/'];
if not(exist(rep))
    mkdir(rep);
end


%%
% Helpers.

fs = 20; % font size
lw = 2; % linewidth
mynorm = @(x)norm(x(:));
saveeps = @(ext)saveas(gcf, [rep name ext '.eps'], 'epsc');
set_gca = @()set(gca, {'XTick' 'YTick' 'ZTick' 'FontSize'}, {[] [] [] fs});

%%
% Parameter of the problem

% relationship between n and N
Nval = @(n)floor( (log(n)/(2*B0)).^(2/r0) );
% resolution of the rho. For real problems, it should be a function of the number n of samples.
N = 15;
% Noise factor, eta=1 <=> no noise. In practice, eta in [.8,.9] 
eta = .9;

%%
% Display clean histograms.

[~,~,p,A] = perform_sampling(name, 0, t);
q = 256;
[Phi,X] = meshgrid(linspace(0,pi,q),linspace(-A,A,q));
H = p(Phi,X);
imwrite(rescale(H), [rep namesvg '-histo.png'], 'png');

if 0
    n = 500000;
    [X0,Phi,p,A] = perform_sampling(name, n, t);
    q = 50;
    H0 = hist3([X0 Phi]', q, q);
    clf; imageplot(H0);
    xlabel('x'); ylabel('\theta'); axis on;
    set_gca();
    saveeps( '-histo' );
end

%%
% Original rho to recover.

rho0 = compute_rho(name, N, t);
clf; bar3(rho0); axis tight;
colormap jet(256);
set_gca();
saveeps( '-rho' );

%%
% Estimator by soft thresholding.

soft_thresh = @(x,gamma) max(0, 1-gamma./max(1e-10,abs(x))) .* x;
Estimate = @(rho,N,n,epsilon,lambda,Vinf) soft_thresh(rho, ...
    lambda * Vinf * sqrt( 2*log(N*(N+1)/epsilon) / n ) );

% relaxation of the threshold
lambda = 1; 
lambda = .5; 
% parameter in the threshold, should be in (0,1)
epsilon = 1;

%%
% Examples of estimation for a few values of n.

nList = [10000 50000 100000 500000];
% nList = [10000 500000];
for i=1:length(nList)
    % number of samples.
    n = nList(i);
    % number of estimated atoms
    N = Nval(n);
    % need to compute L^inf norm of regressors
    [~,Vinf] = perform_rho_deconvolution([],[],eta, N);
    % clean rho
    rho0 = compute_rho(name, N, t);
    % sample
    [X0,Phi] = perform_sampling(name, n, t);
    X = sqrt(eta)*X0 + sqrt((1-eta)/2)*randn(n,1);
    % perform deconvolution
    rho = perform_rho_deconvolution(X,Phi,eta, N);
    % estimate    
    rho1 = Estimate(rho,N,n,epsilon,lambda,Vinf);
    % error
    e = mynorm(rho1-rho0)/mynorm(rho0);
    fprintf('--> n=%d, RelErr=%.3f\n', n, e);
    % display
    clf; bar3(rho1); axis tight;
    colormap jet(256);
    set_gca(); 
    saveeps( ['-n' num2str(round(n/1000)) 'k'] );
end

