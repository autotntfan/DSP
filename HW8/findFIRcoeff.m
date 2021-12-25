function h = findFIRcoeff(P)

% This function finds the FIR coefficients

W = myChebyPol(length(P)-1);  

u = W\P; 
u = u.';

hLeft = fliplr(u(2:end)/2);
hRight = u(2:end)/2;
h = [hLeft, u(1), hRight];






