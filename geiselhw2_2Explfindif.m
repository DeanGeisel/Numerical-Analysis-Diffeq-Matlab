%%%%%%geiselhw2_2Explfindif
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to approximate the IBVP given in problem 2 of
%%%%%%homework 2 at time level 150 using explicit finite difference method 
%%%%%%with 2 different delta t and then graph the approximations.
%%%%%%It will be accompanied by script for an implicit method.
%--------------------------------------------------------------------------
format long
%Establish given values and parameters
h=0.1 ;      %given delta x value
p=[0.0008 0.0011];   %2 given delta t values
tlvl=150;    %requested number of iterations on t "time level"
v=10 ;   %given parameter for the problem
K=4.5 ;  %given paremeter for the problem
x0=0;   %x lower boundary condition
xf=1;   %x upper boundary condition
xmesh=((xf-x0)/h)+1
for g=1:1:2
    %calculate coefficients for use in constructing matrix A
a=(p(g)*K)+(h*p(g)*v);
b=(h^2)-2*p(g)*K-h*p(g)*v;
c=p(g)*K;
    %construct extended tridiagonal A
A=zeros(xmesh-2,xmesh);

for i=1:1:xmesh-2;
    A(i,i:i+2)=[a b c];
end
A=(1/(h^2))*A;
    %construct U vector
U=zeros(xmesh,1);
U(1)=x0;     %given x value at 0 for all time t
U(xmesh)=xf ; %given x value at 1 for all time t
for n=1:1:tlvl
U(2:xmesh-1)=A*U;
end
compareU(:,g)=U
end

x=linspace(0,1,xmesh);
hold on
plot(x,compareU(:,2))
