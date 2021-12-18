% EE5630 HW7 solution: Optimal FIR design, Part I
%
% Yi-Wen Liu
% Dec. 15, 2021
clear; close all;
omega_p = 0.3*pi;
omega_s = 0.4*pi;
x_p = cos(omega_p);
x_s = cos(omega_s);
L = 12; % result will be Type I FIR with 2N+1 samples

figure(1);
line([-1 x_s],[0 0],'linewidth',1.5); hold on;
line([x_p, 1],[1 1],'linewidth',1.5);
set(gca,'ylim',[-0.2 1.2]);
xlabel('x');

% Manually select L+2 locations on the x-axis in the ascending order, 
% with the following constraints:
% x(1) = -1 and/or x(L+2) = 1
% one of them, say x(6) = x_s;
% the next one, say x(7) = x_p;
df = pi/(L+1);
ff = pi:-df:0;
xx = cos(ff); xx = xx(:);
% Change the one closest to x_p to exactly x_p
[dis,ind] = min(abs(xx-x_p));
xx(ind) = x_p;
xx(ind-1) = x_s;

%% Polynomial fitting
% Now, find a polynomial p(x) with the following constraints for all x values in the array xx:
% If x = xx(k) <= x_s then p(x) = (-1)^k * delta, because  it is in the STOPBAND
% If x = xx(k) >= x_p then p(x) = (-1)^k * delta + 1, because it is in the PASSBAND
% Reference: Oppenheim & Schafer Eq. (7.113)
A = zeros(L+2,L+2);
y = zeros(L+2,1); % y(1:L+1) represents the polynomial, y(end) is delta

%% STUDENTS: you should set up the A matrix for solving Ay = b
% ...
% ...
% ...

%% setting up the b vector, which represents the desired response
b = zeros(L+2,1);
b(ind:end) = 1; % xx(ind) = x_p marks the starting place for the passband.

y = A\b;

delta = abs(y(end));
P = y(1:end-1);

% Plot the polynomial
figure(1);
line([-1 x_s],[delta delta],'linestyle','--','color','red');
line([x_p 1],[1+delta 1+delta],'linestyle','--','color','red');
line([-1 x_s],[-delta -delta],'linestyle','--','color','red');
line([x_p 1],[1-delta 1-delta],'linestyle','--','color','red');


xplotloc = -1:0.005:1; 
len = length(xplotloc);
p_x = zeros(size(xplotloc)); 
for j = 1:len
    thisterm = 1;
    for m = 1:L+1
        p_x(j) = p_x(j) + P(m)*thisterm;
        thisterm = thisterm * xplotloc(j);
    end
end
plot(xplotloc,p_x);

%% Find the FIR coefficients
% P(x) = sum_{k=0:L} a_k x^k = h[0] + 2 sum_{k=1:L} h[k] cos(k arccos(x)) = ...
%    h[0]+2*sum_{k=1:L} T_k(x), with the summation index going from 1 to L.
% The next few lines should find (h[0]... h[L]) from 
% polynomial coeff P(1:L+1)
h = zeros(2*L+1,1);
W = myChebyPol(L);  %% STUDENTS: please write this function. 
% 
% W = myChebyPol(N) should return an (N+1)x(N+1) matrix, in which the
% (m+1)th column represents the mth-order Chebyshev polynomial: 
% T_m(x) = w(1,m+1) + w(2,m+1)x + w(3,m+1)x^2 +... +w(m+1,m+1)x^N for all
% m = 0 to N.
% ...
% ...
% ...

%% Converting from the polynomial back to the filter coeffcients
u = W\P; % the result should be u = [h[0] 2h[1] 2h[2] ... 2h[L]]^T.
h(1) = u(1); % h[0]
h(2:L+1) = u(2:end)/2;

%% Then enforce even symmetry
h = circshift(h,L);
h(1:L) = h(end:-1:L+2);

figure(2);
freqz(h);

setFontSizeForAll(12)
