function [y] = solution_k(x,t,r,C0,S0,S1,N)

y=2/(r+1)*C0(1)*(r+(1-r)*x);

for j=1:N
    kn=2*pi*j;
    br = C0(j+1)*( (r+(1-r)*x).*cos(kn*x)-kn*t*(1-r)*sin(kn*x) ) + (S0(j)-(1-r)*S1(j))*sin(kn*x);
    y=y+4*exp(-kn^2*t/2)*br/(r+1);
end
end

