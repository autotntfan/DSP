clear
close all

k = 28800;
b = 1500;
m = 800;
f = logspace(-1,9);
s = 1j*2*pi*f;
H = (b*s+k)./(m*s.^2+b*s+k);
mag = abs(H);
phs = angle(H);
semilogx(f,mag)
grid on
[MaxValue,ind] = max(mag);
f(ind)
GroupDelay = -diff(phs);
figure
semilogx(f(1:end-1),GroupDelay)
figure
freqz([b k],[m b k])