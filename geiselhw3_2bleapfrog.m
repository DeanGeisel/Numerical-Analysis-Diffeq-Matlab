%%%%%%geiselhw3_2bleapfrog
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to approximate the advection equation given in 
%%%%%%problem 2 of homework 3 by using the Leapfrog method.
%%%%%%It will be accompanied by script for the upwind and LaxWendroff method.
%--------------------------------------------------------------------------
format long
%Establish parameters
h=0.01;   %x step
k=0.5*h; %time step from problem
mesh=(1/h)+1; %number of nodes
a=2; %coefficient from equation

%--------------------------------------------------------------------------
%Create Matrix A
sup=1;
sub=-1;
A=zeros(mesh-1);
for i=1:mesh-2
    A(i,i+1)=sup;  %main diagonal
end
for i=1:mesh-2
    A(i+1,i)=sub; %sub diagonal
end
A=(-a*k/h)*A;
%--------------------------------------------------------------------------
%Create U at initial time
xvals=linspace(0,1,mesh);
Un=zeros(mesh-1,1);
Un2=zeros(mesh-1,1);
for i=1:mesh-1
    Un(i)=1+sin(2*pi*xvals(i+1)); %problem setup for initial time
end
Un2(1)=Un(1)-(a*k/h)*(Un(1)-1);
 for i=2:mesh-1
     Un2(i)=Un(i)-(a*k/h)*(Un(i)-Un(i-1));
 end
%--------------------------------------------------------------------------
%Iterate to desired time length
t=[0.05,0.1,0.8];
Uhold=zeros(mesh,3); %place for separate solutions
Uhold(1,:)=1; %for initial x value of 1
fvect=zeros(mesh-1,1);
for i=1:3
    iter=t(i)/k;
    for j=1:iter
        fvect(mesh-1)=-(a*k/h)*Un(1);
        Ustor=Un2
        Un2=A*Un2+(a*k/h)*eye(mesh-1,1)+fvect+Un; %addition for j-1 and j+1 term
        Un=Ustor
    end
    Uhold(2:mesh,i)=Un2;
end
subplot(3,1,1)
  plot(xvals,Uhold(:,1))
  subplot(3,1,2)
  plot(xvals,Uhold(:,2))
  subplot(3,1,3)
  plot(xvals,Uhold(:,3))