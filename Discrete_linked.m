function [SOL] = Discrete_linked(D,m,T,r,PLOT)
%% Applies the discrete model given a finite number of samples D and r

% D vector of sample
% m binning parameter, also determines the values we compute at
% T stopping time
% r parameter linking the BCs
% x values at which we wish to evaluate the solution
% PLOT option to plot the initial data and solution, takes values TRUE/FALSE

%%% ASSUMES pronbability distribution - makes sure samples appropriately
%%% scaled.

% bin data
D=D(:);
dx=1/m; N=length(unique(D));
xmesh=0:dx:1; 
data=histc(D,xmesh)/N;
data=data(1:m);

% parameters for discrete model
alpha=-r/(r+1);
gamma=alpha;
beta=-1/(r+1);
delta=beta;

h = 1/(m+1); 
xmesh = h:h:1-h; 
S = (data(1) + data(end));
dataExtended = [S*(-alpha)
                 data;
                 S*(-beta)]; 
% Scale the initial data to ensure that the integral of data is one.
% Here we will define the `integral' to be: sum(data)*h.            
dataExtended=dataExtended/sum(dataExtended)/h; 
xmeshExtended = 0:h:1;
data = dataExtended(2:end-1);

%% set up the Four Corners Matrix
tic
M = gallery('tridiag',m,-1,2,-1);
M = full(M);
M(1,1) = M(1,1) + gamma; 
M(1,m) = M(1,m) + alpha;
M(m,1) = M(m,1) + beta;
M(m,m) = M(m,m) + delta;

%% Solve by simply computing the matrix exponential.
p = expm(-T * (1/2) * M / h^2)*data;

%% Extend the solution to the two `ghost nodes' at 0 and at n+1. 
% Define S == (y(node 1) + y(node n)).
% Our method of extension DEFINES 
% node 0   value: to be  r/(r+1) of S
% node n+1 value: to be  1/(r+1) of S
% so of course the exact ratio of the left boundary value 
% to the right boundary value is always r, as desired.
S = (p(1) + p(m));
SOL = [S*(-alpha)
             p;
             S*(-beta)]; 
%% Plot the solution, i.e. the estimated density
if PLOT=='TRUE'
    figure
    LW = 'LineWidth';
    FS = 'FontSize';
    plot(xmesh,data,':')
    hold on
    plot(xmeshExtended,SOL,'--k',LW,2) 
    xlabel('x', FS, 16)
    ylabel('Density', FS, 16)
    legend('Initial Data', 'Density Estimate','location','best')
end
end