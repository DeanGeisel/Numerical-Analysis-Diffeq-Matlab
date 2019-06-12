%%%%%%geiselhw3_1ccrknikol
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to approximate the IBVP given in 
%%%%%%problem 1c of homework 3 using the Crank-Nikolson method and then plot the
%%%%%%results at t=0.5
%%%%%%It will be accompanied by script for the explicit method.
%--------------------------------------------------------------------------
%Establish parameters
mesh=41;    %41 points in x and y directions
int=1;   %x and y intervals are both of lenght 1
h=int/(mesh-1);   %delta x and delta y
y0=0;    %x and y values at boundaries
yf=0;
x0=0;
xf=0;
k=0.05; %delta t from CFL should be less than 0.00625
eps=0.05; %given as coeff on Laplacian
pts=linspace(0,1,mesh);%discretized x and y are equal
%--------------------------------------------------------------------------
%Create A to correspond to time step n+1
%from LeVeque book on making sparse Matrix p.64
n=mesh-2;
I=eye(n);
e=ones(n,1);
T=spdiags([e ((-2*(h^2)/(eps*k))-4)*e e],[-1 0 1],n,n);
S=spdiags([e e],[-1 1],n,n);
A=-1*(kron(I,T)+kron(S,I));
%--------------------------------------------------------------------------
%Create B to correspond to time step n
%from LeVeque book on making sparse Matrix p.64
n=mesh-2;
I=eye(n);
e=ones(n,1);
P=spdiags([e ((2*(h^2)/(eps*k))-4)*e e],[-1 0 1],n,n);
S=spdiags([e e],[-1 1],n,n);
B=(kron(I,P)+kron(S,I));
%--------------------------------------------------------------------------
%Create Vector Un
for p=1:n%ycoordinate
    for q=1:n%xcoordinate
        U(p,q)=sin(pi*pts(q+1))*sin(2*pi*pts(p+1));
    end
end
Un=reshape(U,[],1);
figure (1)
surf(U)
%--------------------------------------------------------------------------
%iterate to desired time length
tfin=0.5; %desired time for the solution
iter=tfin/k;
for b=1:iter
    Un=B*Un;
    Un=A\Un;
end
%--------------------------------------------------------------------------
%reconstruct to matrix on xy plane
U0=zeros(41,41);
Uinside=reshape(Un,39,39);
U0(2:40,2:40)=Uinside;
figure (3)
surf(U0)