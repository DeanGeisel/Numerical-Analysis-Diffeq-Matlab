function [fwd,bck,cen,Dap,err]=geiFDMerror(f,x0,h)

%%%by Dean Geisel
%%%for Dr. Mohamed Sulman
%%%in MTH 7260

%--------------------------------------------------------------------------

%%%The purpose of geiFDMerror is to calculate the
%%%absolute error in the forward, backward, center, and D3 approximations of the
%%%derivatives of a function. The inputs are 'f' (a function on x),'x0' (the
%%%x value chosen for function value comparison), and 'h', the step size.

%%%It is necessary to establish x as a symbolic variable first. Therefore,
%%%before running the function, use the command, 'syms x'.

%--------------------------------------------------------------------------

val=[x0-2*h x0-h x0 x0+h];          %values for evaluation points

for n=1:4
    x=val(n) ;                       %evaluate function at the points and create
    fn(n)=subs(f) ;                   %vector of those values
end

fwd=(fn(4)-fn(3))/h;                  %value of forward approximation
bck=(fn(3)-fn(2))/h;                 %value of backward approximation
cen=(fn(4)-fn(2))/(2*h);              %value of center approximation
Dap=(2*fn(4)+3*fn(3)-6*fn(2)+fn(1))/(6*h);  %value of D3 approximation
x=x0;                                     %set value of x as x0 for eval using sub
fp=diff(f);                            %find derivative of f
fptru=subs(fp);                       %find the value of f' at x0

as=abs(fwd-fptru);                     %absolute error of forward approx.
bs=abs(bck-fptru);                    %absolute error of backwards approx.
cs=abs(cen-fptru);                    %absolute error of center approx.
ds=abs(Dap-fptru);                     %absolute error of D3 approx.

a=double(as);
b=double(bs);
c=double(cs);
d=double(ds);
format long;
err=[a  b  c  d];

