function [rho,Vinf] = perform_rho_deconvolution(X,Phi,eta, N)

% perform_rho_deconvolution - estimate rho matrix from noisy samples
%
%   [rho,Vinf] = perform_rho_deconvolution(X,Phi,eta, N);
%
%   Vinf is a (N,N) array where Vinf(j,k) is |f_jk^eta|
%
%   If you only want to estimate Vinf, simply set X=[], Phi=[].
%
%   Copyright (c) 2012 Gabriel Peyre

n = length(X);

Ns = 4096;
Xmax = 10;

% load the pattern functions
persistent x;
persistent f1;
persistent svg_eta; 
persistent svg_Xmax;
if isempty(svg_eta) || svg_eta~=eta || ...
   isempty(svg_Xmax) || svg_Xmax~=Xmax || ...
   length(x)~=Ns || ...
   size(f1,1)~=Ns || ...
   size(f1,2)~=N  || ...
   size(f1,3)~=N 
    [~,f1, x] = load_patterns(N, eta, Ns, Xmax);
    svg_eta = eta;
    svg_Xmax = Xmax;  
end

Vinf = squeeze( max(abs(f1),[],1) );

if isempty(X) || isempty(Phi)
    rho = [];
    return;
end

X = clamp(X, -Xmax*.95,Xmax*.95);
% interpolation
fij = @(xi,j,k)interp1(x, f1(:,j+1,k+1), xi, 'spline');
% estimate rho
rho = zeros(N,N);
for j=0:N-1
    for k=0:N-1
        rho(j+1,k+1) = sum( fij(X/sqrt(eta),j,k) .* exp( -1i*(j-k)*Phi ) ) / n;
    end
end
rho = real(rho);