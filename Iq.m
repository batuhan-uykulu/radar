clear;
clc
close all;
%%size
% parameter assignment
f0  = 1e9;                 % center f
bw  = 10e6;                % bandwidth 
pd  = 1e-6;                % PULSE DURATİON
prp = 290e-6;              % pulse repetition period 
prf = 1/prp;               % pulse repetition freq 
fs = 64e6;                 % smplinf freq
c = 3e8;                   % speed of light
d = 1e3; v = 24; n = 2^4; rcs = 1; A = .1 ; pt = 100; lamb = c / f0; Ts = 1 / fs;

%%
[IQ, t] = IQgen(bw, prp, Ts,pd);       

pr = rxpower(pt, rcs, A, d, lamb);                              % recived signal power
dp = dpmatrix(d, v, n, pd, prp, IQ, t, fs, c, lamb, pr);          % doppler bank 
rbin = t.*c/2;                                                  % range bin
[vm, rng] = dp_rang(dp, t, rbin, prf, c)                        % range&speed measure;
a = dp(8,:) + dp(1,:)+ dp(end,:);
% plot(t,real(a))
% mfd = matchedfilter(a, IQ);
plot(20*log(abs(a)))
% aa = abs(IQ);
% X = [aa aa aa aa aa];
% L =  length(X);
% Y = fft(X);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% % P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(L/2))/L;
% subplot(2,1,1)
% plot(f,P1) 
% subplot(2,1,2)
% plot(X)


%%
function [vm, rng] = dp_rang(dp, t, rbin, prf, c)
    [a, I] = max(abs(dp(1,:)));  % { initial range measure
    est_r  = rbin(I)       % }
    ts     = dp(:, I);
    [Pxx,F] = periodogram(ts,[],256,prf,'centered');
%     plot(F,10*log10(Pxx))
%     grid
%     xlabel('Frequency (kHz)')
%     ylabel('Power (dB)')
    [Y,k] = max(Pxx);
    lambda = c/1e9;
    vm =  dop2speed(F(k)/2,lambda);
    rng = est_r - vm * t(I);
end
function pr = rxpower(pt, rcs, A, R, lamb)    
    pr = pt * rcs * A^2 / ( 4 * pi * R^4 * lamb ^ 2 ); % radar eq
end
function dp = dpmatrix(d, v, n, pd, prp, IQ, t, fs, c, lamb, pr)
    fd = 2*v/lamb; % dopler freq 
    a  = zeros(n,length(t));
    IQ = IQ .* exp(-1i*fd*t); % ad dopler freq
    for i = 1:n
        dd = 2*d + pd * v  +  2 * (i -1) * prp * v; %calc distance every pulse
        a(i,:) = (IQ) * exp(-2*1i*pi*dd/lamb);      % add phase
        m      = round(dd*fs/(c - v));
        adz    = zeros(1, m);   
        temp   = cat(2, adz, a(i,:));               % add delay 
        temp   = temp(1:length(t));
        temp   = awgn(temp, 20)*sqrt(pr/100);       % ad noise - adjust power
%         a(i, :) = temp;
        a(i, :) = matchedfilter(temp, IQ);
    end
%     plot(abs(temp))   
    dp = a;
end
function mfd = matchedfilter(rx, lfm)
    r    = xcorr(lfm, rx);
    r    = flip(r);
    lx   = (length(r));
    half = ceil(lx/2);
    mfd  = r(half : end)/1e3;
end
function [IQ, t] = IQgen(bw, prp, Ts,pd)
    t = 0:Ts:1*prp-Ts;  
    % time line for rfo
    c = bw/pd;  
    Ichirp = cos(2*pi*(c/2*t.^2));    % chirp gen
    Qchirp = sin(2*pi*(c/2*t.^2)); 
    lo    = zeros(1, length(t));        
    n     = round(pd/Ts);
    lo(1:n) = 1;                       % lo gen
    I = lo .*  Ichirp;                 % mixer
    Q = lo .*  Qchirp;                 % mixer
    IQ = I + 1i * Q;
end