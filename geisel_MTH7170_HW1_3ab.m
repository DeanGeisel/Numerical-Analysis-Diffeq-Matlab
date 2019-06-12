%%%Script for 3a and b homework set 1 Due 2/11/19
%%%by Dean Geisel
%%%for Dr. Mohamed Sulman
%%%in MTH 7260

%-------------------------------------------------------------------------

%%%The purpose of this script is to solve the matrix equation found in
%%%problem 3 for the Vector of U values that form the numerical solution.
%%%Also to plot the numerical solution and exact solution on the same
%%%graph.

%-------------------------------------------------------------------------

h=pi/20;
M=((pi/2)/h)-1;
m=round(M,0);
coef=[-1-(h/2) 2+2*(h^2) -1+(h/2)];
U0=-0.3;
Ue=-0.1;
A=zeros(m,m);
A(1,1:2)=coef(2:3);
A(m,m-1:m)=coef(1:2);
for n=2:1:m-1
    A(n,n-1:n+1)=coef;
end
A=(1/(h^2))*A;
F=zeros(1,m);
for j=1:1:m
F(j)=-cos(h*j);
end
F(1)=F(1)-(coef(1)*(U0/(h^2)));
F(m)=F(m)-(coef(3)*(Ue/(h^2)));
[u]=LSGramschmidt(A,F);
U=[U0;u;Ue]

exact=@(x)-(1/10)*(sin(x)+3*cos(x));
L=linspace(0,pi/2,m+2);
hold on
fplot (exact,[0 pi/2])
scatter (L,U)
title('exact solution vs. numerical solution')
xlabel('x')
ylabel('u(x)')
legend ('exact solution', 'numerical solution')