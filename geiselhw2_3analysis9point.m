%%%%%%geiselhw2_3analysis9point
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to analyze the approximation to the Poisson problem given in 
%%%%%%problem 3 of homework 2 using 5-point and 9-point stencils for the Laplacian. 
%%%%%%It will be accompanied by script for the fivepoint and ninepoint methods.
%--------------------------------------------------------------------------
%Establish parameters
testval=[17,33,65,129];
for i=1:1:2
mesh=testval(i); %number of evaluated points
xinit=-1; %beginning x value of interval
xfin=1; %end x value of interval
xint=(xfin-xinit);%given interval of x values
yinit=-1; %beginning x value of interval
yfin=1; %end x value of interval
yint=(yfin-yinit);%given interval of x values
hx=xint/(mesh-1);%given h value
hy=yint/(mesh-1);
h=xint/(mesh-1);       %matching hx and hy
H(i)=h;
n=mesh-2;%calculate the number of intervals as n both intervals same
y0=0;
yf=0;
x0=0;
xf=0;
%--------------------------------------------------------------------------
%create A
%from LeVeque book on making sparse Matrix p.64
I=eye(n);
e=ones(n,1);
T=spdiags([4*e -20*e 4*e],[-1 0 1],n,n);
R=spdiags([e 4*e e],[-1 0 1],n,n);
S=spdiags([e e],[-1 1],n,n);
A=(kron(I,T)+kron(S,R))/(6*(h^2));
figure(1)
spy(A)
%--------------------------------------------------------------------------
%create right side vector Fmod
Fmod=-4*ones(n^2,1);
Fmod(1:n)=Fmod(1:n)-6*x0;
Fmod(n^2-(n-1):n^2)=Fmod(n^2-(n-1):n^2)-6*xf;
for p=0:1:(n-1)
Fmod(1+p*n)=Fmod(1+p*n)-4*y0;
Fmod(n+p*n)=Fmod(n+p*n)-4*yf;
end

%--------------------------------------------------------------------------
%solving for solution vector
U=inv(A)*Fmod;
Uarray=zeros(n+2);

for g=0:1:n-1;
   for k=1:1:n;
         Uarray(k+1,g+2)=U(g*n+k);
   end
end
Uarray(:,1)=y0;
Uarray(:,mesh)=yf;
Uarray(1,:)=x0;
Uarray(mesh,:)=xf;
figure (2)
surf(Uarray)

Lx=linspace(xinit,xfin,mesh);%finding values to be evaluated for exact
Ly=linspace(yinit,yfin,mesh);
clear exactval;
exactval=zeros(mesh,mesh);
for k=1:1:mesh
    x=Lx(k);
for g=1:1:mesh
    y=Ly(g);
    syms m
Sum=symsum(((-8*sin((2*m-1)*pi/2))/(((2*m-1)*pi/2)^3))*cosh(((2*m-1)*pi/2)*y)*cos(((2*m-1)*pi/2)*x),m,1,2);
exact=(2-2*(x^2)+Sum);
exactval(k,g)=exact;
 %given exact solution
end
end

E=abs(exactval-Uarray);
L1(i)=h*norm(E,1);
L2(i)=(h^(1/2))*norm(E);
end
figure (2)
loglog(L1,H,'-.o')
for p=1:1:1
    conrate1(p)=(log10(abs(L1(p+1)-L1(p))))/(log10(abs(H(p+1)-H(p))));
    conrate2(p)=(log10(abs(L2(p+1)-L2(p))))/(log10(abs(H(p+1)-H(p))));
end
L1
L2
conrate1
conrate2