%%%Script for 1b homework set 1 Due 2/11/19
%%%by Dean Geisel
%%%for Dr. Mohamed Sulman
%%%in MTH 7260

%--------------------------------------------------------------------------

%%%The purpose of the script is to verify if the derived fourth-order 
%%%accurate formula for u'(x) is fourth order accurate.

%--------------------------------------------------------------------------
format long
x=1/2;            %From Problem setup
h=[0.1;0.01;0.005;0.001;0.0005;0.0001];%choice of various h values for comparison
u=@(y)cos(pi*y);%given test function
up=@(y)-pi*sin(pi*y);%derivative of test function
upaprx=@(x,h)(1/(12*h))*(-3*u(x-h)-10*u(x)+18*u(x+h)-6*u(x+2*h)+u(x+3*h));%approximation
E=zeros(6,1);%Empty Error Vector
for n=1:1:6
    E(n)=abs(up(x)-upaprx(x,h(n)));%Calculate Error for each h
end
  acc=E./(h.^4);%divide error by h^4 the expected order of accuracy
T=table(h,E,acc);%create table of data
T





