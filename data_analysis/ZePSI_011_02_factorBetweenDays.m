%
clearvars

%% Comparison between day 2 (xa, ya) and day 1 (xb, yb). T = 140 K

[xa, ya, parama] = eprload("../data/raw/ZePSI-E-011004-a.DTA");
[xb, yb, paramb] = ...
    eprload("../data/raw/ZePSI-E-011004/ZePSI-E-011004-00001.DTA");

% clf
% h = ScrollableAxes();
% plot(h, 1:4096, 1:301, ya');
% hold on
% plot(h, 1:4096, 1:301, yb')

noiseRange = [3, 40; 210, 273];
[snra, ppa, noisea] = getSNR(1:275, ya(792, 1:275), noiseRange);
[snrb, ppb, noiseb] = getSNR(1:275, yb(792, 1:275), noiseRange);

clf
plot(yb(792, :))
hold on
plot(ya(792, :))
xline(noiseRange(1))
xline(noiseRange(2))
xline(noiseRange(3))
xline(noiseRange(4))
xlim([1, 301])
% SNRa/SNRb: 0.44345   ppa/ppb: 0.52056   noisea/noiseb: 1.1739
% Therefore all the measurements after day 1 should be multiplied by 2
% The reason of this is unknown (different freezing?)
disp("SNRa/SNRb: " + string(snra/snrb) + ...
    "   ppa/ppb: " + string(ppa/ppb) + ...
    "   noisea/noiseb: " + string(noisea/noiseb));

%% Comparison between day 3 (xa, ya) and day 2 (xb, yb). T = 220 K

aFolder = dir("../data/raw/ZePSI-E-011008-a/*.DTA");
[xa, ya, parama] = averagetreprscanselexsys(aFolder, 5, 9.6e9, 50);
xa{2} = xa{2}/10;

[xt, yt, paramt] = eprload("../data/raw/ZePSI-E-011008-a/ZePSI-E-011008-a-00001.DTA");
xa{2} = xa{2}/10;

%% Baseline correction field domain

Opt = {};
Opt.order = 1;
Opt.width = 0.1;
[y1, bl1] = deal(yt(:, 43:end-1));
for it = 1:nt
    [y1(it, :), bl1(it, :), p] = ...
        subtractbaseline(x{2}', ya(it, 43:end-1), Opt);
end

iPlot1 = 4;
clf
h = ScrollableAxes();
plot(h, x{2}, x{1}, ya(:, 43:end-1)');
hold on
plot(h, x{2}, x{1}, bl1');
blRegionWidth = Opt.width*(max(x{2}) - min(x{2}));
xline(min(x{2}) + blRegionWidth)
xline(max(x{2}) - blRegionWidth)

%% Baseline correction time domain

Opt.order = 0;
Opt.range = [7000, 8100];
[y2, bl2] = deal(y1);
[nt, nB] = size(y1);
for ib = 1:nB
    [y2(:, ib), bl2(:, ib)] = ...
        subtractbaseline(x{1}', y1(:, ib), Opt);
end

iPlot1 = 1;
clf
h = ScrollableAxes();
plot(h, x{1}, x{2}, y1');
hold on
plot(h, x{1}, x{2}, bl2');
xline(Opt.range(1))
xline(Opt.range(2))

%%
aa = load("../data/processed/ZePSI-E-011.mat");

figure(3)
% clf
% h = ScrollableAxes();
% plot(h, xa{2}, transpose(1:4096), ya(:, 1:258))
% hold on
% plot(h, aa.x{2}, transpose(1:4096), aa.y{8})

ytfin = y2;
clf
plot(1:258, ytfin(792, :)/max(ytfin(792, :)))
hold on
plot((1:258) + 8, aa.y{8}(792, :)/max(aa.y{8}(792, :)))
xlim([1, 258])
%% 
noiseRange = [3, 40; 210, 273];
[snra, ppa, noisea] = getSNR(1:258, ytfin(792, :), noiseRange);
[snrb, ppb, noiseb] = getSNR(1:258, aa.y{8}(729, :), noiseRange);

% SNRa/SNRb: 2.2958   ppa/ppb: 2.5637   noisea/noiseb: 1.1167
% Therefore all the measurements after day 2 should be additionally 
% divided by ~2.6
% The reason of this is unknown (thawing and refreezing of the sample
% maybe)
disp("SNRa/SNRb: " + string(snra/snrb) + ...
    "   ppa/ppb: " + string(ppa/ppb) + ...
    "   noisea/noiseb: " + string(noisea/noiseb));

for it = 1:nt
    yy = ytfin(it, :);
    yyint(it) = trapz(yy);
end
test = sum(ytfin, 2);
clf
plot(yyint)
hold on
plot(test)