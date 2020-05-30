clear;close all;

load('Example_data.mat')
n=10^4;
D=rand(n,1);
D=ceil(D*m);
D=LT(D);  % LT stores inverse of cdf - D can be any input data
% r=1.7604;
x=0:0.001:1; % Points at which we approximate PDF
PLOT = 1;

% % compute stopping time via kde code from reference [3]
% [bandwidth,~,~,~]=kde(D,10*n,0,1);
% T=bandwidth^2
T=3*10^(-4); % good guess for this example if user does not have above code
%%
SOL = Cts_samples(D,T,r,x,PLOT);
hold on
plot(x,polyval(P,x),'r','linewidth',2)
legend('Initial Data', 'Density Estimate','True Density','location','best')