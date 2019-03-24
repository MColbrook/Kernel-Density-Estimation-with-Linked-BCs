function [SOL] = Cts_func(f0,T,r,x,PLOT)
%% Applies the continous model given an intial function f0 and r

% f0 initial function
% T stopping time
% r parameter linking the BCs
% x values at which we wish to evaluate the solution
% PLOT option to plot the initial data and solution, takes values TRUE/FALSE

N=max(20,ceil(sqrt( 8*log(10)/(pi^2*T) ))); % sum N+1 terms

C0=zeros(N+1,1);
S0=zeros(N,1);
S1=zeros(N,1);

% integrate to find the constants - of course this can be replaced with
% other numerical integration algorithms

C0(1)=quadgk(@(x) f0(x),0,1,'AbsTol',1e-8,'RelTol',1e-8,'MaxIntervalCount',1e8);

for j=1:N
    kn=2*pi*j;
    C0(j+1) = quadgk(@(x) cos(kn*x).*f0(x),0,1,'AbsTol',1e-8,'RelTol',1e-8,'MaxIntervalCount',1e8);
    S0(j)   = quadgk(@(x) sin(kn*x).*f0(x),0,1,'AbsTol',1e-8,'RelTol',1e-8,'MaxIntervalCount',1e8);
    S1(j)   = quadgk(@(x) x.*sin(kn*x).*f0(x),0,1,'AbsTol',1e-8,'RelTol',1e-8,'MaxIntervalCount',1e8);
end

SOL=solution_k(x,T,r,C0,S0,S1,N);

if PLOT=='TRUE' 
    figure
    x0=0:0.01:1;
    LW = 'LineWidth';
    FS = 'FontSize';
    plot(x0,f0(x0),':')
    hold on
    plot(x,SOL,'--k',LW,2) 
    xlabel('x', FS, 16)
    ylabel('Density', FS, 16)
    legend('Initial Data', 'Density Estimate','location','best')
end
end