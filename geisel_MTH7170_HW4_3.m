%%%%%%geiselhw4_3
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to solve the Matrix equation resulting from
%%%%%%the Finite Element Method on problem 3 of HW 4.
%--------------------------------------------------------------------------
%Establish Matrices 
y=linspace(0,1,4);

A=[3 -3 0 0;    %%%%%%for a(
    -3 6 -3 0;
    0 -3 6 -3;
    0 0 -3 3];

D=[1/2 -1/2 0 0;    %%%%for (u',v'x)
    -1/2 2 -3/2 0;
    0 -3/2 4 -5/2;
    0 0 -5/2 5/2];

    
S=A+D;  %%%Left side Matrices No right side boundary adjustments 
 S(1,2:4)=0
 S(1,1)=1

 
 F=[0;0;0;2] %%%%0 plus 2*phi(n)

 A
 D
 S
 F
 
  U=S\F
exact=@(x) 2*log(1+x)
 hold on
 plot(y,U)
 fplot(exact)
