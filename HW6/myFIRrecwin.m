function h = myFIRrecwin(omega_c, N)

n = -N:1:N;
h = sin(omega_c*n)./(n*pi);
h(N+1) = omega_c/pi;
h = h';

    