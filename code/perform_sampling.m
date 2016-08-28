function [X,Phi,p,A] = perform_sampling(name, n, t)

% perform_sampling - perform the sampling of the distribution
%
%   [X,Phi,p,A] = perform_sampling(name, n, t)
%
%   (X(i), Phi(i)) are sampled from the distribution indicated by name.
%   t is an (optional) parameter of the distribution)
%
%   p(phi,x) is P(X|Phi), which has approximate range in [-1,1]. 
%
%   Copyright (c) 2012 Gabriel Peyre

switch name
    case 'vacuum'
        p = @(phi,x)exp(-x.^2) / sqrt(pi);
        A = 5; b = 1.2;
    case 'single-photon'
        p = @(phi,x)2*x.^2 .* exp(-x.^2) / sqrt(pi);
        A = 10; b = .6;
    case 'coherent'
        q0 = t;
        p = @(phi,x)exp(-(x-q0*cos(phi)).^2) / sqrt(pi);
        A = 6; b = 1.2;
    case 'squeezed'
        error('Not yet implemented');
    case 'thermal'
        beta = t;
        p = @(phi,x)exp(-x.^2 * tanh(beta/2) ) * sqrt(tanh(beta/2)/pi);
        A = 3 / sqrt(tanh(beta/2)/pi); 
        b = 1.1 * sqrt(tanh(beta/2)/pi);        
    case 'shrodinger-cat'
        q0 = t;
        A = @(phi,x) exp(-(x-q0*cos(phi)).^2) + exp(-(x+q0*cos(phi)).^2);
        B = @(phi,x) 2*exp(-x.^2-q0^2*((cos(phi)).^2)).*cos(2*q0*x.*sin(phi));
        p = @(phi,x) ( A(phi,x) + B(phi,x) ) / (2*sqrt(pi)*(1+exp(-q0^2)));
        A = 10; b = 1.2;      
    otherwise 
        error('Unknown');
end

X = []; Phi = []; 
if isempty(n) || n<1
    return;
end

% check integrability
for iter=1:10
    x = linspace(-A,A,n); phi = rand*pi;
    e = sum(p(phi,x))*2*A/n;
    if e<.95 || max(p(phi,x))>b
        warning('Sampling issue.');
        A = A*2; b = max(p(phi,x))*1.5;
    end
    % plot(x, p(phi,x));
end

% random sampling using rejection
% use mini-batches of length N
s = 0;
while length(X)<=n
    x = A*(2*rand(n,1)-1); y = b*rand(n,1); phi = rand(n,1)*pi;
    I = find( y<=p(phi,x) );
    X   = [X;x(I)];
    Phi = [Phi;phi(I)];
    s = s+1;
    if s>100
        error('Sampling problem');
    end
end
Phi = Phi(1:n);
X = X(1:n);