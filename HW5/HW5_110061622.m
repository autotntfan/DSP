clear; close all; clc;
% For Q7, read an audio file
[signal,fs] = audioread('boxing.wav');
% clip 
signal = signal(1:2.5e5,1);
%soundsc(signal,fs)

N = 4:1:20; % order
Omega_c = 0.2*pi; % Cut-Off Frequency ( Big Omega_c )
color = {'red','yellow','green','blue','magenta'};
Legend = repmat({''},1,18);
for ii = 1:length(N)
    
    [Td, DEN, NUM] = FunctionForQ1(N(ii), Omega_c, Omega_c);
    [h,w] = freqz(NUM,DEN,65536);
    
    figure(2)
    hold on
    if mod(N(ii),4) == 0
        plot(w/pi,abs(h),color{N(ii)/4},'DisplayName',['N = ' num2str(N(ii))])
    else
        plot(w/pi,abs(h),'black','HandleVisibility','off')
    end
    if N(ii) == 6
        figure(1)
        plot(w/pi,abs(h),'HandleVisibility','off')
        xlabel('Normalized frequency \times\pi')
        ylabel('|H(j\Omega)|')
        title('6th-order filter with cut-off freq. 0.2\pi')
    end
end
grid on
hold on
plot(0.28,0.01,'r*','DisplayName','Critical Point')
legend show
xlabel('Normalized frequency \times\pi')
ylabel('|H(j\Omega)|')

pb = 0.7071;
sb = 0.01;
lb = 0.2*pi;
ub = 0.28*pi;
% 得到新N
[N, DEN, NUM] = FunctionForQ5(lb,ub,pb,sb,Td);
[h,w] = freqz(NUM,DEN,65536);
figure
plot(w/pi,abs(h),'black')
xlabel('Normalized frequency \times\pi')
ylabel('|H(j\Omega)|')
title(['N = ' num2str(N)])
TestForFilter(N, lb, 1);
legend('Our','Built-in')

y = FunctionForQ6(NUM,DEN,signal);

X = fftshift(fft(signal));
Y = fftshift(fft(y));
freq = linspace(-fs/2,fs/2,length(signal));
figure
plot(freq,abs(X),'black')
hold on
plot(freq,abs(Y),'red')
legend('original','LPF')
xlabel('frequency (Hz)')

fs = 1e3;
t = 0:1/fs:5;
TestSignal = 3*cos(2*pi*80*t)+3*sin(2*pi*90*t)+3*sin(2*pi*100*t)+3*cos(2*pi*110*t)+3*sin(2*pi*120*t)+3*sin(2*pi*130*t);
y = FunctionForQ6(NUM,DEN,TestSignal);
y_prime = filter(NUM,DEN,TestSignal);
figure
plot(t,y,'black',t,y_prime,'red--')
xlabel('time (s)')
legend('ourself','built-in')

X = fftshift(fft(TestSignal));
Y = fftshift(fft(y));
freq = linspace(-fs/2,fs/2,length(TestSignal));
figure
plot(freq,abs(X),'black')
hold on
plot(freq,abs(Y),'red')
legend('original','LPF')
xlabel('frequency (Hz)')


function [Td, DEN, NUM] = FunctionForQ1(N, omega, Omega)
    % N : filter order
    % omega : frequency in discrete-time domain
    % Omega : frequency in continuous-time domain
    Td = 2*tan(omega/2)/Omega;
    % uniformly distributed angles
    angle = pi/2 + (0.5 : 1 : N-0.5).*pi./N;
    s_p = Omega*exp(1j*angle); % look up end of text book
    % map s to z
    z_p = (1+Td*s_p/2)./(1-Td*s_p/2); % Bilinear Transform
    DEN = poly(z_p); % Denominator A(z)

    z_zero = -1*ones(1,N);
    NUM = poly(z_zero); % Numerator B(z)

    DCgain = 2^N/sum(DEN); % B(z)/A(z) evaluated at z = 1, i.e. @ omega = 0.
    NUM = NUM/ DCgain;
end

function [N, DEN, NUM] = FunctionForQ5(lb, ub, pb, sb, Td)
    % passband ≦ |H(e^jw)| ≦ 1,            0  ≦ |w| ≦ lower bound
    %             |H(e^jw)| < stopband, upbound ≦ |w| ≦ π
    % lb : lower bound 
    % ub : upper bound
    % pb : passband
    % sb : stopband
    
    N = ceil(log(((1/sb)^2-1)/((1/pb)^2-1))/(2*log(tan(ub/2)/tan(lb/2))));
    Omega_c = 2*tan(lb/2)/Td/(1/pb^2-1)^(1/2/N);
    [~, DEN, NUM] = FunctionForQ1(N, lb, Omega_c);
end

function y = FunctionForQ6(NUM, DEN, x)
    % NUM : numerator   B(z)
    % DEN : denominator A(z)
    ytermlen = length(DEN(2:end));
    xtermlen = length(NUM);
    x = reshape(x,[1,length(x)]);
    y = zeros(size(x));
    for ii = 1:length(x)
        if ii == 1 
            y(ii) = NUM(1)*x(1)/DEN(1);
        elseif 1<ii && ii<xtermlen
            y(ii) = (sum(NUM(1:ii).*fliplr(x(1:ii)))-sum(DEN(2:ii).*fliplr(y(1:ii-1))))/DEN(1);
        else
            y(ii) = (sum(NUM.*fliplr(x(ii-xtermlen+1:ii)))-sum(DEN(2:end).*fliplr(y(ii-ytermlen:ii-1))))/DEN(1);
        end
    end
end

function [b,a] = TestForFilter(N, omega_c, flag)
    [b,a] = butter(N,omega_c/pi);
    if flag == 0
        figure
    else
        hold on
    end
    [h,w] = freqz(b,a,65536);
    plot(w/pi,abs(h),'r--')
    xlabel('Normalized frequency \times\pi')
    ylabel('|H(j\Omega)|')
    title(['N = ' num2str(N)])
end

