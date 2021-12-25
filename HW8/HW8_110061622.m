clc; clear; close all;

%% parameters
omega_p = 0.2*pi;
omega_s = 0.24*pi;
x_p = cos(omega_p);
x_s = cos(omega_s);
L = 32; % polynomial order

%% init.
figure(1);
line([-1 x_s],[0 0],'linewidth',1.5); hold on;
line([x_p, 1],[1 1],'linewidth',1.5);
set(gca,'ylim',[-0.2 1.2]);
xlabel('x');
df = pi/(L+1);
ff = pi:-df:0;
xx = cos(ff); xx = xx(:);
[dis,ind] = min(abs(xx-x_p));
xx(ind) = x_p;
xx(ind-1) = x_s;

%% main algorithm
notFinish = 1; % decide if optimum is found
count = 0; % count iteraion times
maxIt = 80; % max allowed iteration
while notFinish
    L_plus_what = length(xx) - L; % literally L plus "what"
    count = count + 1;
    if count > maxIt
        fprintf('Iteration stops!\n\n');
        break;
    end
    %% Polynomial fitting
    A = zeros(L+L_plus_what,L+L_plus_what);
    A(:,1) = 1;
    for m = 2:L+L_plus_what-1
        A(:,m) = A(:,m-1).* xx;
    end
    A(:,L+L_plus_what) = (-1).^(1:L+L_plus_what);
    b = zeros(L+L_plus_what,1);
    b(ind:end) = 1; 
    y = A\b;
    delta = abs(y(end));
    P = y(1:end-1);
    xplotloc = -1.0:0.001:1.0; 
    len = length(xplotloc);
    p_x = zeros(size(xplotloc)); 
    for j = 1:len
        thisterm = 1;
        for m = 1:L+L_plus_what-1
            % thisterm:x^0,x^1,x^2,x^3....x^L
            p_x(j) = p_x(j) + P(m)*thisterm;
            thisterm = thisterm * xplotloc(j);
        end
    end
    %% STUDENTS: WRITE THE FUNCTION newArgument()
    newArgument(p_x, xplotloc, L, delta, x_p, x_s);
    try
        [arg, mag] = newArgument(p_x, xplotloc, L, delta, x_p, x_s);
    catch
        clear; close all;
        error('Students have to write the function newArgument()');
    end
    %% check if the polynomial is optimal
    IsOptimal = isOptimal(mag, delta);
    IsntOptimal = ~IsOptimal;
    notFinish = IsntOptimal;
    %% if the new argument is almost identical to 
    %  the old argument, then we say optimum is found
    try
        if max(abs(sort(arg) - xx.')) < 0.001*pi
            notFinish = 0; 
        end
    catch
    end
    %% feed the new argument into the next iteration
    xx = sort(arg); xx = xx.';
end
figure(1)
%% plotting
line([-1 x_s],[delta delta],'linestyle','--','color','red');
line([x_p 1],[1+delta 1+delta],'linestyle','--','color','red');
line([-1 x_s],[-delta -delta],'linestyle','--','color','red');
line([x_p 1],[1-delta 1-delta],'linestyle','--','color','red');
plot(xplotloc,p_x);
% scatter(arg, mag, 'filled');

%% find FIR coefficients
try
    h = findFIRcoeff(P);
catch
    clear; close all;
    error('Students have to include myChebyPol() from HW7');
end
   
figure,
[H,w] = freqz(h, 1, 2048);
plot(w, abs(H)); hold on
scatter([omega_p,omega_s], [1-delta,delta],'filled')
text([1.1*omega_p,1.1*omega_s], [1.05-1.05*delta,1.3*(L/6)*delta],...
     {['|H(\omega_p)|,   \omega_p = ',num2str(omega_p/pi),'\pi'],...
      ['|H(\omega_s)|,   \omega_s = ',num2str(omega_s/pi),'\pi']});
line([omega_s, pi],[delta delta],'linestyle','--','color','red');
line([0 omega_p],[1+delta 1+delta],'linestyle','--','color','red');
line([0 omega_p],[1-delta 1-delta],'linestyle','--','color','red'); hold off

%% delta
fprintf('delta is %.3f\n', delta);