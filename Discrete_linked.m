function [SOL2] = Discrete_linked(D,m,T,r,PLT)
%% Applies the discrete model given a finite number of samples D and r

% D vector of sample
% m binning parameter, also determines the values we compute at
% T stopping time
% r parameter linking the BCs
% x values at which we wish to evaluate the solution
% PLT option to plot the initial data and solution, takes values 1/0

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
B = gallery('tridiag',m,-1,2,-1);
% vectors for Sherman–Morrison formula
u=sparse(m,1); v=u;
u(1)=gamma; u(end)=beta;
v(1)=1; v(end)=1;

%% Solve by backwards Euler and Sherman–Morrison formula

K=ceil(T/(2*h^2));
Delta_t=T/K;
p2=data;

B=speye(m,m)+Delta_t*B/(2*h^2);
v=Delta_t*v/(2*h^2);

U=B\u;

for j=1:K
    ZZ=B\p2;
    p2=ZZ-U*(transpose(v)*ZZ)/(1+transpose(v)*U);
end


%% Extend the solution to the two `ghost nodes' at 0 and at n+1. 
% Define S == (y(node 1) + y(node n)).
% Our method of extension DEFINES 
% node 0   value: to be  r/(r+1) of S
% node n+1 value: to be  1/(r+1) of S
% so of course the exact ratio of the left boundary value 
% to the right boundary value is always r, as desired.

S2 = (p2(1) + p2(m));
SOL2 = [S2*(-alpha)
             p2;
             S2*(-beta)];


%% Plot the solution, i.e. the estimated density
if PLT==1
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