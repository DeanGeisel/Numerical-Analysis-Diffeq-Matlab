%%%%%%geiselhw2_1Eulerforward
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to approximate the IVP given in problem 1 of
%%%%%%homework 2 using Eulers forward approximation(explicit), graph the 
%%%%%%approximation, then verify first order error.
%%%%%%It will be accompanied by script for RK2 and RK4.
%--------------------------------------------------------------------------
%calculate U values

h=0.05;                  %%%%step size for t
t0=0;                    %%%%initial value for t
tf=1;                    %%%%final value for t
m=(tf-t0)/h;             %%%%Number of intervals of t
t=linspace(t0,tf,m+1);     %%%%Discretize t 
        %%%%%note:t(1)=t0, t(m+1)=tf
u=zeros(m+1,1);             %%%%Vector for u values
u(1)=1 ;                  %%%%Initial condition for u in the u vector
f=zeros(m,1);              %%%%Vector for function values
for n=1:1:m              %%%%Recursive formula for Euler forward(explicit)
    f(n)=exp(t(n)-u(n));
    u(n+1)=u(n)+h*f(n);
end

%-------------------------------------------------------------------------
%Plot solution curve

disp (u)
figure(1)
plot (t,u)

%-------------------------------------------------------------------------
%Verify first degree accuracy

%%%Modified script from HW 1
%%%The purpose of this script is to caclulate L1 and L2 Norm errors in the
%%%numerical solution for problem 1 of HW set n. The n values are
%%%17,33,65,129

%-------------------------------------------------------------------------
testval=[17,33,65,129];
for i=1:1:4
mesh=testval(i); %number of evaluated points
int=(tf-t0);%given interval of x values
hdep=int/(mesh-1);%calculate step size based on interval number
H(i)=hdep; %Store H values for covergence rate

tdep=linspace(t0,tf,mesh);     %%%%Discretize t 
        %%%%%note:tdep(1)=t0, tdep(mesh)=tf
U=zeros(mesh,1);             %%%%Vector for u values
U(1)=u(1) ;                  %%%%Initial condition for u in the u vector
fdep=zeros(mesh-1,1);              %%%%Vector for function values
for n=1:1:mesh-1              %%%%Recursive formula for Euler forward(explicit)
    fdep(n)=exp(tdep(n)-U(n));
    U(n+1)=U(n)+hdep*fdep(n);
end

exact=@(x)log(exp(x)+exp(1)-1); %given exact solution
L=linspace(t0,tf,mesh);%finding values to be evaluated for exact
clear exactval
for k=1:1:mesh
exactval(k)=exact (L(k));
end
E=abs(exactval'-U);
L1(i)=hdep*norm(E(:),1);
L2(i)=(hdep^(1/2))*norm(E(:));
end
figure (2)
loglog(L1,H,'-.o')
for p=1:1:3
    conrate1(p)=(log10(abs(L1(p+1)-L1(p))))/(log10(abs(H(p+1)-H(p))));
    conrate2(p)=(log10(abs(L2(p+1)-L2(p))))/(log10(abs(H(p+1)-H(p))));
end
L1
L2
conrate1
conrate2

