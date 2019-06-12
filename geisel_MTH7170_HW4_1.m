%%%%%%geiselhw4_1
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to solve the Matrix equation resulting from
%%%%%%the Finite Element Method on problem 1 of HW 4.
%--------------------------------------------------------------------------
%Establish Matrices 
x=[0; 1/6; 1/2; 4/5; 1]
A=[1 0 0 0 0;
    0 9 -3 0 0;
    0 -3 19/3 -10/3 0;
    0 0 -10/3 25/3 -5;
    0 0 0 0 1]
M=[0 0 0 0 0;
    1/36 1/6 1/18 0 0;
    0 1/18 19/90 1/20 0;
    0 0 1/20 1/6 1/30;
    0 0 0 0 0]
F=[1;5/6;1/2;1/5;0] %f(x)=1-x
W=(A+M)\M*F
for i=1:5
U(i)=W(i)+2*x(i);
end
U
plot(x,U)
title('Approximation of u for -u"+u=1+x') 
xlabel('x')
ylabel('u')
