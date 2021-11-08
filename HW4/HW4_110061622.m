clear
close all

% parameters
k = 28800; % elastic constant [Nt/m]
b = 1500;  % damping coefficient [kg/s]
m = 800;   % mass [kg]
L = 8;     % one cycle length of bump [m]
D = 0.05;  % height of bump [m]

% (a) plot H(j2pif) v.s. f
f = logspace(-1,4,1e4);
s = 1j*2*pi*f;
H = (b*s + k)./(m*s.^2 + b*s + k);
mag = abs(H);
phs = angle(H);
figure
subplot(211)
semilogx(f,mag2db(mag))
xlabel('frequency (Hz)')
ylabel('Mag (dB)')
subplot(212)
semilogx(f,rad2deg(phs))
xlabel('frequency (Hz)')
ylabel('degree')
% find the max value at which fequency 
[~,ind] = max(mag);
MaxValueAt = f(ind)

GroupDelay = -diff(phs);
figure
semilogx(f(1:end-1),GroupDelay)
title('group delay')
xlabel('frequency (Hz)')
ylabel('degree/s')

% part 2 variables
v = 10;      % velocity [m/s]
T = 0.0001;  % sampling interval [s]
t = 10;t = 0:T:t-T; % time axis
y = v*t;   
figure
plot(t,D*sin(2*pi*y/L),'black')
hold on
plot(t,abs(D*sin(2*pi*y/L)),'b--')
hold on
g = (D*sin(2*pi*y/L) + abs(D*sin(2*pi*y/L))) ./2;
plot(t,g,'r')
legend('g(y)','|g(y)|','(g(y)+|g(y)|)/2')
x1 = g;
% initial condition for x2
x2 = [0 0];
alpha = T*b + 2*m;
for ii = 3:length(g)
   x2(ii) = x1(ii)*T*b/alpha + x1(ii-1)*2*T^2*k/alpha - x1(ii-2)*T*b/alpha - x2(ii-1)*(2*T^2*k-4*m)/alpha - x2(ii-2)*(2*m-T*b)/alpha;
end
figure
plot(t,x1)
hold on
plot(t,x2)
title('displacement')
legend('x1','x2')
xlabel('time (s)')
ylabel('displacement (m)')
MaxDisp = max(x2)
% another way to determine x2
y = filter([T*b 2*T^2*k -T*b],[T*b+2*m 2*T^2*k-4*m 2*m-T*b],x1,[]);
figure
plot(t,x1)
hold on
plot(t,y)
title('displacement')
legend('x1','x2')
% Plot spectrum of g or x1
fs = 1/T;
G = fftshift(fft(g));
freq = linspace(-fs/2,fs/2,length(g));
figure
plot(freq,abs(G))
title('spectrum of g')
xlabel('frequency (Hz)')
ylabel('Mag')
xlim([-20 20])
% set sampling freq
fs = 2e1;
t = 0:1/fs:10;
h = (487*sqrt(111)*sin(9*sqrt(111)/16.*t)+1665*cos(9*sqrt(111)/16.*t)).*exp(-15/16.*t)/888;
figure
plot(t,h)
xlabel('time (s)')
ylabel('amp')
title('h(t)')
figure
% use freqz to draw response
freqz(h,fs,length(t),fs)
[amp,w] = freqz(h,fs,length(t),fs);
figure
semilogx(w,mag2db(abs(amp)))
title('amplitude response')
xlabel('frequency (rad/s)')
ylabel('amplitude (dB)')

% draw the spectrum to check bandwidth
% H = fftshift(fft(h));
% freq = linspace(-fs/2,fs/2,length(h));
% figure
% plot(freq,abs(H))
% xlabel('freq (Hz)')
% xlim([-5 5])