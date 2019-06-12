%%%Script for 2a and b homework set 1 Due 2/11/19
%%%by Dean Geisel
%%%for Dr. Mohamed Sulman
%%%in MTH 7260

%-------------------------------------------------------------------------

%%%The purpose of this script is to solve the matrix equation found in
%%%problem 2 for the Vector of U values that form the numerical solution.
%%%Also to plot the numerical solution and exact solution on the same
%%%graph.

%-------------------------------------------------------------------------

h=.05;%given h value
m=(1/h)-1;%calculate the number of intervals as m
coef=[-1 2+4*h^2 -1];%The coefficients on U estimates
U0=0;%Boundary condition
U1=2;%Boundary condition
A=zeros(m,m);%Empty A matrix
A(1,1:2)=coef(2:3);%first row of Matrix A
A(m,m-1:m)=coef(1:2);%last row of Matrix A
for n=2:1:m-1
    A(n,n-1:n+1)=coef; %fill the tridiagonal matrix from the coeffi
end
A=(1/h^2)*A;%multiply matrix by 1/h^2
F=zeros(1,m);%%establish vector for right side of matrix equation
for j=1:1:m
F(j)=4*h*j;%evaluate for values of F
end
F(m)=F(m)+(2/h^2); %correction for U0 and U1
[u]=LSGramschmidt(A,F);%solve Matrix equation using previously written code
U=[U0;u;U1]%correction for U0 and U1

exact=@(x)exp(2)*(exp(2*x)-exp(-2*x))/(exp(4)-1)+x; %given exact solution
L=linspace(0,1,m+2);%finding values to be evaluated for exact
hold on
fplot (exact,[0 1])
scatter (L,U)
title('exact solution vs. numerical solution')
xlabel('x')
ylabel('u(x)')
legend ('exact solution', 'numerical solution')

