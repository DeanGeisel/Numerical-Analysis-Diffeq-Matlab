%%%%%%geiselhw4_2
%%%%%%by Dean Geisel
%%%%%%for Dr. Mohamed Sulman
%%%%%%in MTH 7170
%--------------------------------------------------------------------------
%%%%%%This script is written to solve the Matrix equation resulting from
%%%%%%the Finite Element Method on problem 2 of HW 4.
%--------------------------------------------------------------------------
%Establish Matrices 
x=linspace(0,pi/2,11);
A=zeros (11);
A(1,1)=20/pi;
A(11,11)=20/pi;
for i=1:10
    A(i+1,i)=-20/pi;
    A(i,i+1)=-20/pi;
end
for i=2:10
    A(i,i)=40/pi;
end
M=zeros (11);
M(1,1)=pi/60;
M(11,11)=pi/60;
for i=1:10
    M(i+1,i)=pi/120;
    M(i,i+1)=pi/120;
end
for i=2:10
    M(i,i)=pi/30;
end
C=zeros (11);
C(1,1)=1/2;
C(11,11)=1/2;
for i=1:10
    C(i+1,i)=-1/2;
    C(i,i+1)=-1/2;
end
for i=2:10
    C(i,i)=1;
end
F=zeros(11,1);
for i=1:11
    F(i)=cos(x(i))-0.6;
end
S=A+C+(2*M);
S(1,2:11)=0;
S(1,1)=1;
M(1,1:11)=0;

Con=zeros(11,1);
Con(11)=0.3;
T=-M*F+Con;
W=S\T

for i=1:11
    U(i)=W(i)-0.3
end
    
exact=@(x) (-1/10)*(sin(x)+3*cos(x))%given exact solution
hold on
plot(x,U)
fplot(exact)

 


for k=1:1:11
exactval(k)=exact (x(k));
end
E=abs(exactval'-U);

L2=((pi/20)^(1/2))*norm(E(:))
A
M
C
S

