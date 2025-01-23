clearvars
addpath(genpath('util/'))

%% IMPORT

loadPath = '../data/raw/ZePSI-E-013002/ZePSI-E-013002-024001-RM/F0006CH1.CSV';
importData = readtable(loadPath);
x = table2array(importData(:, 4));
y = table2array(importData(:, 5));
nx = numel(x);

y = y - mean(y(1000:1500));

%% PLOT

x1 = x(1:nx/2);
y1 = y(1:nx/2);
x2 = x(nx/2 + 1:end);
y2 = y(nx/2 + 1:end);

figure()
clf
plot(x, y)
hold on
plot(x1, y1)
plot(x2, y2)
yline(0)

%%
% FFT
tStep = x1(2) - x1(1);
fSampl = 1/tStep;
nzf = 0;  % Zero filling
if nzf ~= 0
    ff = fSampl/nzf*(-nzf/2:nzf/2 - 1);
else
    ff = fSampl/(nx/2)*(-(nx/2)/2:(nx/2)/2 - 1);
end

tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
if nzf ~= 0
    y1(ii, nzf) = 0;  % Zero filling
    y2(ii, nzf) = 0;  % Zero filling
end

fy1 = fft(y1);
fy2 = fft(y2);

figure(4)
clf
plot(ff, abs(fftshift(fy1)), 'o-')
hold on
plot(ff, abs(fftshift(fy2)), 'o-')
xlim(0.3*[-1, 1]*1e9)

