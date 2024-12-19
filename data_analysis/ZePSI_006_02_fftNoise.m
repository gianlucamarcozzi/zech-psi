% Power dependence of the trEPR signal
clearvars

% Import
cd('/home/gianlum33/files/projects/zech_psi/data_analysis/')
expName = '../data/raw/ZePSI-E-006';
measNo = {'003'};
nMeas = numel(measNo);

for jj = 1:nMeas
    filename = [expName measNo{jj} '.DTA'];
    [xRaw{jj}, yRaw{jj}, Param{jj}] = eprload(filename);
end

x = xRaw;

%% Baseline correction field domain

Opt.order = 0;
% rangeNoiseB = [3333, 3342; 3368, 3371];
Opt.width = 0.1;
[y1, bl1] = deal(yRaw);
for jj = 1:nMeas
    [nt, ~] = size(yRaw{jj});
    for it = 1:nt
        [y1{jj}(it, :), bl1{jj}(it, :)] = ...
            subtractbaseline(x{jj}{2}', yRaw{jj}(it, :), Opt);
    end
end

clf
h = ScrollableAxes();
jj = 1;
plot(h, x{jj}{2}, x{jj}{1}, yRaw{jj}');
hold on
plot(h, x{jj}{2}, x{jj}{1}, bl1{jj}');
plot(h, x{jj}{2}, x{jj}{1}, y1{jj}');
blRegionWidth = Opt.width*(max(x{jj}{2}) - min(x{jj}{2}));
xline(min(x{jj}{2}) + blRegionWidth);
xline(max(x{jj}{2}) - blRegionWidth);
legend('2ns, 2048')

%% Baseline correction time domain

clear("Opt")
Opt.order = 0;
Opt.range = [0, x{1}{1}(30)];
[y2, bl2] = deal(yRaw);
for jj = 1:nMeas
    [~, nB] = size(yRaw{jj});
    for ib = 1:nB
        [y2{jj}(:, ib), bl2{jj}(:, ib)] = ...
            subtractbaseline(x{jj}{1}', y1{jj}(:, ib), Opt);
    end
end

clf
h = ScrollableAxes();
jj = 1;
plot(h, x{jj}{1}', x{jj}{2}', y1{jj});
hold on
plot(h, x{jj}{1}', x{jj}{2}', bl2{jj});
plot(h, x{jj}{1}', x{jj}{2}', y2{jj});
xline(Opt.range(2))
legend('2ns, 2048')

%% Filter Butterworth

Opt.Type = 'lowpass';
Opt.Order = 4;
Opt.Cutoff = 2e-2;
y_ = y2{jj}(:, 157);
y3 = filter_butterworth(y_, Opt);
clf
% plot(x{1}{1}, y_)
% hold on
plot(x{1}{1}, y3)
hold on

[bFilt, aFilt] = butter(4, 0.09);
y4 = filter(bFilt, aFilt, y_);
plot(x{1}{1}, y4)

%% FFT

[nt, nB] = size(y2{1});
ntFill = nt*4;
y2Fill = zeros(ntFill, nB);
y2Fill(1:nt, :) = y2{1};

dx = (x{1}{1}(2) - x{1}{1}(1))*1e-9;  % seconds
df = 1/dx;  % Hz
f{1} = (-ntFill/2:ntFill/2 - 1)*(df/ntFill);

for ib = 1:nB
    yfft{1}(:, ib) = fft(y2Fill(:, ib));
    yfftSh{1}(:, ib) = fftshift(yfft{1}(:, ib));
end

cmap = viridis();
clf
imagesc(f{1}, x{1}{2}, abs(yfftSh{1}))
colormap(cmap)
colorbar()

clf
h = ScrollableAxes();
jj = 1;
plot(h, f{jj}', x{jj}{2}', abs(yfftSh{jj}));
legend('2ns, 2048')

[bFilt, aFilt] = butter(1, 1e8/(max(f{1})/2));
clf
freqz(bFilt, aFilt, [], max(f{1}))

dataOut = filter(bFilt, aFilt, yfft{1}(:, 157));
dataOutSh = fftshift(dataOut);
clf
tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact')
nexttile()
plot(f{1}, real(yfftSh{1}(:, 157)))
hold on
plot(f{1}, imag(yfftSh{1}(:, 157)))
nexttile()
plot(f{1}, real(dataOutSh))
hold on
plot(f{1}, imag(dataOutSh))

%%

yFiltered = ifft(dataOut);
clf
plot(x{1}{1}, y2{1}(:, 157))
hold on
plot(x{1}{1}, real(yFiltered))
hold on
plot(x{1}{1}, imag(yFiltered))

%% Boxcar integration

% Determine the range as the window in which signal > max/2
for jj = 1:nMeas
    bMM = 157;  % Index of the maximum in the field domain
    [val(jj), tMM(jj)] = max(y2{jj}(:, bMM));
    thresholdVal = val(jj)*0.7;
    for it = tMM(jj):numel(y2{jj}(:, bMM))
        if y2{jj}(it, bMM) > thresholdVal
            tStop(jj) = it;
        else
            break
        end
    end
    for it = tMM(jj):-1:1
        if y2{jj}(it, bMM) > thresholdVal
            tStart(jj) = it;
        else
            break
        end
    end
end

tiledlayout(3, 1, 'TileSpacing', 'tight', 'Padding', 'compact')
nexttile()
h = ScrollableAxes();
jj = 1;
plot(h, x{jj}{1}', x{jj}{2}', y2{jj});
xline(x{jj}{1}(tStart(jj)))
xline(x{jj}{1}(tStop(jj)))
legend('4ns, 4096')
nexttile()
h = ScrollableAxes();
jj = 2;
plot(h, x{jj}{1}', x{jj}{2}', y2{jj});
xline(x{jj}{1}(tStart(jj)))
xline(x{jj}{1}(tStop(jj)))
legend('2ns, 4096')
nexttile()
h = ScrollableAxes();
jj = 3;
plot(h, x{jj}{1}', x{jj}{2}', y2{jj});
xline(x{jj}{1}(tStart(jj)))
xline(x{jj}{1}(tStop(jj)))
legend('2ns, 2048')

for jj = 1:nMeas
    y{jj} = mean(y2{jj}(tStart(jj):tStop(jj), :), 1);
end

%% SNR ppAmp and sigmaNoise
% sigmaNoise does not sensibly change with power, ppAmp changes

noiseRangeB = [3410, 3419; 3450, 3455];
clf
for jj = 1:nMeas
    [SNR(jj), ppAmp(jj), noiseLev(jj)] = ...
        getSNR(x{jj}{2}, y{jj}, noiseRangeB);
    plot(x{jj}{2}, 1.*y{jj})
    hold on
end
xlim([min(x{1}{2}) max(x{1}{2})])
legend('4ns, 4096', '2ns, 4096', '2ns, 2048')
xline(noiseRangeB(1, 1))
xline(noiseRangeB(1, 2))
xline(noiseRangeB(2, 1))
xline(noiseRangeB(2, 2))
% exportgraphics(gcf, '/home/gianlum33/files/inbox/mwpw_comparison.png')

%% SNR ppAmp and sigmaNoise from traces
% sigmaNoise does not sensibly change with power, ppAmp changes

for jj = nMeas:-1:1
    [~, ppAmp2(jj), ~] = getSNR(x{jj}{1}, y2{jj}(:, bMM), [0, 1]);
    [~, ~, noiseLev2(jj)] = getSNR(x{jj}{1}, y2{jj}(:, end), [0, 16000]);
    SNR2(jj) = ppAmp2(jj)/noiseLev2(jj);
end
% xlim([min(x{1}{2}) max(x{1}{2})])
% legend('20dB, SNR = 220', '25 dB, SNR = 240', '30 dB, SNR = 239', ...
%     '35dB, SNR = 149', '43 dB, SNR = 68')
% exportgraphics(gcf, '/home/gianlum33/files/inbox/mwpw_comparison.png')
clf
plot([43, 35, 30, 25, 20], SNR)
hold on
plot([43, 35, 30, 25, 20], SNR2*5)

