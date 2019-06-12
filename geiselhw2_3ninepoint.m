%%%%%%geiselhw2_3ninepoint
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to approximate the Poisson problem given in 
%%%%%%problem 3 of homework 2 using 9-point stencil for the Laplacian and 
%%%%%%then graph the approximations.
%%%%%%It will be accompanied by script for the ninepoint method.
%--------------------------------------------------------------------------
%Establish parameters
h=0.1;   %delta x and delta y
mesh=(2/h)+1;
n=mesh-2;
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
figure (1)
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