%%
% Test for the pattern functions.

eta = .9;
N = 20;
Ns = 4096*4;
Xmax = 5;

[f,f1, x] = load_patterns(N, eta, Ns, Xmax);
clf;
q = 4;
A = { [0 0] [4 2]  [2 1] [5 5] };
for i=1:length(A)
    j = A{i}(1); k = A{i}(2);
    subplot(ceil(length(A)/2),2,i);
    plot(x, [f(:,j+1,k+1) f1(:,j+1,+k+1)]/pi); axis tight;
    title(['f_{' num2str(j) ',' num2str(k) '}']);
end
legend('f/\pi', 'f1/\pi');