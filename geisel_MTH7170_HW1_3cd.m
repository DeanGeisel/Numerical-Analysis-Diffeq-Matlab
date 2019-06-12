%%%Script for 3c homework set 1 Due 2/11/19
%%%by Dean Geisel
%%%for Dr. Mohamed Sulman
%%%in MTH 7260

%-------------------------------------------------------------------------

%%%The purpose of this script is to caclulate L1 and L2 Norm errors in the
%%%numerical solution for problem 3 of HW set 1. The n values are
%%%17,33,65,129

%-------------------------------------------------------------------------
testval=[17,33,65,129];
for i=1:1:4
mesh=testval(i); %number of evaluated points
intinit=0; %beginning x value of interval
intfin=pi/2; %end x value of interval
int=(intfin-intinit);%given interval of x values
h=int/(mesh-1);%given h value
H(i)=h;
M=((pi/2)/h)-1;
m=round(M,0);;%calculate the number of intervals as m
coef=[-1-(h/2) 2+2*(h^2) -1+(h/2)];%The coefficients on U estimates
U0=-0.3;%Boundary condition
Ue=-0.1;%Boundary condition
A=zeros(m,m);%Empty A matrix
A(1,1:2)=coef(2:3);%first row of Matrix A
A(m,m-1:m)=coef(1:2);%last row of Matrix A
for n=2:1:m-1
    A(n,n-1:n+1)=coef; %fill the tridiagonal matrix from the coeffi
end
A=(1/h^2)*A;%multiply matrix by 1/h^2
F=zeros(1,m);%establish vector for right side of matrix equation
for j=1:1:m
F(j)=-cos(h*j);%evaluate for values of F
end
F(1)=F(1)-(coef(1)*(U0/(h^2)));
F(m)=F(m)-(coef(3)*(Ue/(h^2)));

[u]=LSGramschmidt(A,F);%solve Matrix equation using previously written code
U=[U0;u;Ue];%correction for U0 and U1

exact=@(x)-(1/10)*(sin(x)+3*cos(x)); %given exact solution
L=linspace(intinit,intfin,mesh);%finding values to be evaluated for exact
clear exactval
for k=1:1:mesh
exactval(k)=exact (L(k));
end
E=abs(exactval'-U);
L1(i)=h*norm(E(:),1);
L2(i)=(h^(1/2))*norm(E(:));
end
loglog(L1,H,'-.o')
for p=1:1:3
    conrate1(p)=(log10(abs(L1(p+1)-L1(p))))/(log10(abs(H(p+1)-H(p))));
    conrate2(p)=(log10(abs(L2(p+1)-L2(p))))/(log10(abs(H(p+1)-H(p))));
end
L1
L2
conrate1
conrate2