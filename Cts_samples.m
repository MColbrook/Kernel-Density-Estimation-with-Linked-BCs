function [SOL] = Cts_samples(D,T,r,x,PLOT)
%% Applies the continous model given a finite number of samples D and r

% D vector of sample
% T stopping time
% r parameter linking the BCs
% x values at which we wish to evaluate the solution
% PLOT option to plot the initial data and solution, takes values 1/0

N=max(20,ceil(sqrt( 8*log(10)/(pi^2*T) ))); % sum N+1 terms

C0=zeros(N+1,1);
S0=zeros(N,1);
S1=zeros(N,1);
C0(1)=1;

for j=1:N
    kn=2*pi*j;
    C0(j+1) = mean(cos(kn*D));
    S0(j)   = mean(sin(kn*D));
    S1(j)   = mean(sin(kn*D).*D);
end

SOL=solution_k(x,T,r,C0,S0,S1,N);

if PLOT==1
    D=D(:);
    m=100; % default, can change according to histogram binning
    dx=1/m; N=length(unique(D));
    xmesh=0:dx:1; 
    data=histc(D,xmesh)/N;
    data=data/sum(data)*(m+1); 
    figure
    LW = 'LineWidth';
    FS = 'FontSize';
    plot(xmesh,data,':')
    hold on
    plot(x,SOL,'--k',LW,2) 
    xlabel('x', FS, 16)
    ylabel('Density', FS, 16)
    legend('Initial Data', 'Density Estimate','location','best')
end
end

