clear
close all

DIR = './';
[x, fs] = audioread([DIR 'boxing.wav']); 
x = x(1:10*fs); 
x = x(:);

SigLength = length(x); % 

f_cut = 4000; % This Low Pass Filter only reserves parts below 4000Hz
h = fir1(128, f_cut/(fs/2)); % since 0 < Wn < 1, students may use help fir1 to get te idea
h = h(:);
L = length(h); % filter size

tic
y = conv(x,h);
toc

tic
L_zp = SigLength + L - 1; % since the final length is L + SigLengh -1, in advance, we append zero
x_zp = zeros(L_zp,1);
x_zp(1:SigLength) = x;
h_zp = zeros(L_zp,1);
h_zp(1:L) = h;
y_conv = zeros(L_zp,1);
inv_x = flip(x_zp);    % according to formulate definition 
xz = zeros(L_zp,1);
for ii = 1:L_zp
    xz(1:ii) = inv_x(end-(ii-1):end);
    y_conv(ii) = h_zp'*xz;
end
toc

sum(abs(y_conv - y))
%save('yconv.mat','y_conv')