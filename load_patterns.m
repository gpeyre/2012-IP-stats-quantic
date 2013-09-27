function [f,f1, x] = load_patterns(N, eta, Ns, Xmax)

% load_patterns - load the pattern functions
%
%   [f,f1] = load_patterns(N, eta, [Ns, Xmax]);
%
%   f are the clean pattern function
%   f1 are the deconvolved pattern functions
%
%   f(i) is the evaluation of the pattern at location x(i).
%   Ns is the number of sampling locations.
%   The samples are in [-Xmax,Xmax]
%
%   Copyright (c) 2012 Gabriel Peyre

if nargin<3
    Ns = 2048*2;
end
if nargin<4
    Xmax = 10;
end

gamma = (1-eta)/(4*eta);

% functions callbacks
laguerre = @(t,k,alpha)polyval(LaguerreGen(k, alpha), t);
pattern1 = @(t,j,k)pi*(-1i)^(j-k) .* sqrt( 2^(k-j) * factorial(k) / factorial(j) ) .*  ...
                abs(t) .* t.^(j-k) .* exp(-(t.^2)/4) .* laguerre((t.^2)/2, k,j-k);
pattern = @(t,j,k)pattern1(t,max(j,k),min(j,k));

% sampling locations
Tmax = pi/Xmax*Ns/2;
t = 2*Tmax/Ns * [0:Ns/2-1 -Ns/2:-1]'; % sampling frequencies [-Tmax,Tmax]
x = pi/Tmax * [0:Ns/2-1 -Ns/2:-1]'; % sampling locations
x = fftshift( x ); % correct sampling locations

% generate clean patterns
F = zeros(Ns,N,N);
for j=0:N-1
    for k=0:N-1
        F(:,j+1,k+1) = pattern(t,j,k); 
    end
end
f = 1/(2*pi) * 2*Tmax*fftshift( real(ifft(F)), 1 );

% generate deconvolved patterns
u = min(exp(gamma*t.^2), 1e9);
F1 = F .* repmat( u, [1 N N] );
f1 = 1/(2*pi) * 2*Tmax*fftshift( real(ifft(F1)), 1 );

return;

%% BUG %%
for j=0:N-1
    for k=0:N-1
        if mod(j+k,2)==1%impair
            f(:,j+1,k+1) = f(end:-1:1,j+1,k+1);
            f1(:,j+1,k+1) = f1(end:-1:1,j+1,k+1);
        end
    end
end
