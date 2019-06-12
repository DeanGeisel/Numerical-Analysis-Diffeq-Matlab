%%%%%%geiselhw3_2blaxwend
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to approximate the advection equation given in 
%%%%%%problem 2 of homework 3 by using the Lax-Wendroff method.
%%%%%%It will be accompanied by script for the upwind and Leapfrog method.
%--------------------------------------------------------------------------
format long
%Establish parameters
h=0.01;   %x step
k=0.5*h; %time step from problem
mesh=(1/h)+1; %number of nodes
a=2; %coefficient from equation

%--------------------------------------------------------------------------
%Create Matrix A
d=((a^2)*(k^2))/(2*(h^2))+((a*k)/(2*h));  %coefficients on the u approxim.
e=((-(a^2)*(k^2))/(h^2))+1;
f=((a^2)*(k^2))/(2*(h^2))-((a*k)/(2*h));

A=zeros(mesh-1);
for i=2:mesh-2
    A(i,i-1:i+1)=[d e f];  %tridiagonal missing top and bottom row
end
A(1,1:2)=[e f];
A(mesh-1,mesh-2:mesh-1)=[d e];

%--------------------------------------------------------------------------
%Create U at initial time
xvals=linspace(0,1,mesh);
Un=zeros(mesh-1,1);
for i=1:mesh-1
    Un(i)=1+sin(2*pi*xvals(i+1)); %problem setup for initial time
end
%--------------------------------------------------------------------------
%Iterate to desired time length
t=[0.05,0.1,0.8];
Uhold=zeros(mesh,3); %place for separate solutions
Uhold(1,:)=1 ; %for initial x value of 1 on each solution U to top Un
fvect=zeros(mesh-1,1);
for i=1:3
    iter=t(i)/k;
    for j=1:iter
        fvect(mesh-1)=f*Un(1);
        Un=A*Un+d*eye(mesh-1,1)+fvect; %addition for j-1 term and j+1 term
    end
    Uhold(2:mesh,i)=Un;
end
subplot(3,1,1)
  plot(xvals,Uhold(:,1))
  subplot(3,1,2)
  plot(xvals,Uhold(:,2))
  subplot(3,1,3)
  plot(xvals,Uhold(:,3))