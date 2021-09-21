clear;
clc
close all;
%%
% parameter assignment
f0  = 1e9;                
bw  = 25e6;               
prf = 10e3;
pd  = 1/prf;
fs  = 64e6;
d   = 1e3; 
c   = 3e8;
rcs = 1; A = .1 ; pt = 100;
lamb = c / f0;
%%
waveform = phased.FMCWWaveform("SampleRate", fs, "SweepTime", pd, "SweepBandwidth", bw,...  % FMCW gen
    "SweepDirection", "triangle", "NumSweeps", 1);
cw       = step(waveform); 
t = 0 : 1/fs: pd -1/fs;   %time line
%%
n  = ceil(d/c*fs);         % {
rx = zeros(1,length(cw));  %    add delay
rx(n+1:end) = cw(1:end-n); % }
rx = transpose(rx);
pr = rxpower(pt, rcs, A, d, lamb);  % rx power
rx  = awgn(rx,10)*sqrt(pr/1000);    % ad noise adjst power

y  = dechirp(rx, cw);   % dechirp

plot(t, y)

[Pxx,F] = periodogram(y,[],length(y),fs,'centered'); %{ get beat frequency
[p, I] = max(Pxx);                                   %}
r  = beat2range(F(I),bw/(2*pd))                      % get range
% plot(F/1000,10*log10(Pxx/1e-3)); grid;
%%
function pr = rxpower(pt, rcs, A, R, lamb)    
    pr = pt * rcs * A^2 / ( 4 * pi * R^4 * lamb ^ 2);
end
