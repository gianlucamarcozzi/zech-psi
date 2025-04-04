%%
clearvars

%% IMPORT
filename = "../data/raw/ZePSI-E-013001";
[xRaw, yRaw, param] = eprload(filename);
xRaw{1} = xRaw{1}';
xRaw{2} = xRaw{2}';

[nt, nB] = size(yRaw);

%% SCROLLABLE
clf
h = ScrollableAxes('Index', 1);
plot(h, xRaw{1}, xRaw{2}, yRaw');

%%
clf
plot(tt*1e-9, yRaw(:, 9))
hold on
plot(tt*1e-9, 30*sin(tt*1e-9*4.7e5 - 1))


%% FFT FIRST SLICE
tt = xRaw{1};
yy = yRaw(:, 5);
dtt = (tt(2) - tt(1))*1e-9;  % seconds
df = 1/dtt;  % Hz
f = (-nt/2:nt/2-1)*(df/nt);

yfft0 = fft(yy);
yfft = abs(fftshift(yfft0));

clf
plot(f, yfft)
xlim([-1e6, 1e6])
xline(4.7e5)