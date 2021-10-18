clear 
close all
% [x,fs] = audioread('ABBA-flac2wav.wav');
% sound(x,fs)
% sound(x(2:2:end,:),fs);
% sound(x(2:2:end,:),fs/2);
% sound(resample(x,fs/2,fs), fs/2)
% y = x(2:2:end,:)-resample(x,fs/2,fs);
% plot(y)
% y = resample(x,fs/4,fs);
% 
% audiowrite('ABBAfs_4.wav',y,fs/4)

n = 0:1:40;
% h = sqrt(3)*(-1).^(n+1)./(2*pi*(n-1/3));
h1 = sin((n-1/3)*pi)./((n-1/3)*pi);
stem(n,h1)
xlabel('n')
ylabel('amplitude')
figure
freqz(h1);



% 
% % figure
% % fs = 20;
% % t = -20:1/fs:20;
% % freq = linspace(-fs/2,fs/2,length(t));
% % h1 = sin((t-1/3)*pi)./((t-1/3)*pi);
% % z = fftshift(fft(h1));
% % plot(freq,abs(z))
% % figure
% % tol = 1e-6;
% % z(abs(z) < tol) = 0;
% % plot(freq,angle(z)/pi)
% % 
% % t = 0:0.01:1;
% % y = sin(2*pi*4*t);
% % plot(t,y);
% % figure
% % plot(t,ifft(ifftshift(fftshift(fft(y)))))

