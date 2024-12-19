% Calculate SNR for different laser rep rates and for different integration
% windows and different temperatures.
clearvars

% Import
cd('/home/gianlum33/files/projects/zech_psi/data_analysis/')
expName = '../data/raw/ZePSI-E-003';
measNo = {'012', '018'};  % Isaak, X/W
nMeas = numel(measNo);

filename = [expName measNo{1}];
[xRaw{1}, yRaw{1}, yRuns{1}, Param{1}] = loadtreprisaak(filename);

filename = [expName measNo{2} '.DTA'];
[xRaw{2}, yRaw{2}, Param{2}] = eprload(filename);
yRaw{2} = yRaw{2}/Param{2}.NbScansAcc;

% Data correction
freqCorr = 9.6e9;  % Hz
tCorrIndex = 50;  % Index of light flash
for jj = 1:nMeas
    % Laser flash at t = 0
    xCorr{jj}{1} = ...
        xRaw{jj}{1} - xRaw{jj}{1}(tCorrIndex);
    % Correct frequency to freqCorr
    try
        freqCorr = 9.6e9;
        xCorr{jj}{2} = ...
            xRaw{jj}{2}*(freqCorr./Param{jj}.MWFQ);
    catch
        freqCorr = 9.6;
        xCorr{jj}{2} = ...
            xRaw{jj}{2}*(freqCorr./Param{jj}.mwfreq);
    end
end

x = xCorr;

% 
% % Interpolate wrt the first field-axis
% x = xCorr{1};
% for jj = 1:nMeas
%     for it = 1:numel(yRaw{jj}(:, 1))
%         yCorr{jj}(it, :) = interp1( ...
%             xCorr{jj}{2}, yRaw{jj}(it, :), xCorr{1}{2});
%     end
% end
% 
% % Correct field-axis for different length of the frequncy-shifted axis
% for jj = 1:nMeas
%     minFields(jj) = min(xCorr{jj}{2});
%     maxFields(jj) = max(xCorr{jj}{2});
% end
% lowField = max(minFields);
% highField = min(maxFields);
% [~, lowFieldIdx] = min(abs(x{2} - lowField));
% [~, highFieldIdx] = min(abs(x{2} - highField));
% 
% x{2} = x{2}(lowFieldIdx:highFieldIdx);
% for jj = 1:nMeas
%     yCorr{jj} = yCorr{jj}(:, lowFieldIdx:highFieldIdx);
% end

ttilde = 175;


%% Plot spectra

clf
plot(x{1}{2}, yRaw{1}(ttilde, :))
yyaxis right
plot(x{2}{2}, yRaw{2}(ttilde, :))

% xlim(setaxlim(xCorr{2}, 1))
legend('Isaak', 'X/W')
% clf
% plot(1:nB, data(2).yRuns(:, ttilde, 1))
% hold on
% plot(1:nB, data(2).yRuns(:, ttilde, 2))
% plot(1:nB, data(2).yRuns(:, ttilde, 3))
% plot(1:nB, data(2).yRaw(:, ttilde), 'LineWidth', 2)
% yyaxis right
% plot(1:nB, data.y, 'LineWidth', 2)

%% Scrollable traces

clf
tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'compact')
nexttile()
h = ScrollableAxes();
plot(h, x{1}{1}', x{1}{2}', yRaw{1});
legend('Isaak')
% nexttile
% h = ScrollableAxes();
% plot(h, x{2}{1}', x{2}{2}', yRaw{2});
% legend('X/W')

%% Scrollable spectra

clf
tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'compact')
nexttile()
h = ScrollableAxes();
plot(h, x{1}{2}, x{1}{1}, yRaw{1}');
legend('Isaak')
nexttile
h = ScrollableAxes();
plot(h, x{2}{2}, x{2}{1}, yRaw{2}');
legend('X/W')

%% Baseline correction field domain

Opt.order = 0;
% rangeNoiseB = [3333, 3342; 3368, 3371];
Opt.width = 0.1;
[y1, bl1] = deal(yRaw);
for jj = 1:2
    [nt, ~] = size(yRaw{jj});
    for it = 1:nt
        [y1{jj}(it, :), bl1{jj}(it, :)] = ...
            subtractbaseline(x{jj}{2}', yRaw{jj}(it, :), Opt);
    end
end

clf
tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'compact')
nexttile()
h = ScrollableAxes();
jj = 1;
plot(h, x{jj}{2}, x{jj}{1}, yRaw{jj}');
hold on
plot(h, x{jj}{2}, x{jj}{1}, bl1{jj}');
plot(h, x{jj}{2}, x{jj}{1}, y1{jj}');
blRegionWidth = Opt.width*(max(x{jj}{2}) - min(x{jj}{2}));
xline(min(x{jj}{2}) + blRegionWidth);
xline(max(x{jj}{2}) - blRegionWidth);
legend('Isaak')
nexttile
h = ScrollableAxes();
jj = 1;
plot(h, x{jj}{2}, x{jj}{1}, yRaw{jj}');
hold on
plot(h, x{jj}{2}, x{jj}{1}, bl1{jj}');
plot(h, x{jj}{2}, x{jj}{1}, y1{jj}');
blRegionWidth = Opt.width*(max(x{jj}{2}) - min(x{jj}{2}));
xline(min(x{jj}{2}) + blRegionWidth);
xline(max(x{jj}{2}) - blRegionWidth);
legend('X/W')

%% Baseline correction time domain

clear("Opt")
Opt.order = 0;
rangeIsaak = [0, x{1}{1}(24)];
rangeXw = [0, x{2}{1}(30)];
[y2, bl2] = deal(yRaw);
for jj = 1:2
    if jj == 1
        Opt.range = rangeIsaak;
    else
        Opt.range = rangeXw;
    end
    [~, nB] = size(yRaw{jj});
    for ib = 1:nB
        [y2{jj}(:, ib), bl2{jj}(:, ib)] = ...
            subtractbaseline(x{jj}{1}', y1{jj}(:, ib), Opt);
    end
end

tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'compact')
nexttile()
h = ScrollableAxes();
jj = 1;
plot(h, x{jj}{1}', x{jj}{2}', y1{jj})
hold on
plot(h, x{jj}{1}', x{jj}{2}', bl2{jj})
plot(h, x{jj}{1}', x{jj}{2}', y2{jj})
xline(rangeIsaak(2))
legend('Isaak')
nexttile()
h = ScrollableAxes();
jj = 2;
plot(h, x{jj}{1}', x{jj}{2}', y1{jj})
hold on
plot(h, x{jj}{1}', x{jj}{2}', bl2{jj})
plot(h, x{jj}{1}', x{jj}{2}', y2{jj})
xline(rangeXw(2))
legend('X/W')

%% Different integration windows
%{

t1 = 100;  % 1us
nIntegration = 45;
unitIntegration = 20;  % 0.4us
SNR = zeros(nT, nIntegration, max(nRepRate));
yMean = zeros(nT, nIntegration, nB);
clf
for iT = 1:nT
    nrr_ = nRepRate(iT);
    for jj = 1:nrr_
        for ii = 1:nIntegration
            t2 = t1 + (ii - 1)*unitIntegration;
            yMean(iT, jj, :, ii) = mean(data(iT, jj).yRaw(:, t1:t2), 2);
            [~, ppAmp(iT, jj, ii), sigNoise(iT, jj, ii)] ...
                = getSNR(1:nB, yMean(iT, jj, :, ii), [0, 30; 320, 350]);
        end
    end
end

for iT = 1:nT
    nrr_ = nRepRate(iT);
    for jj = 1:nrr_
        % Normalize for averaging
        ppAmpAvg(iT, jj, :) = ...
            ppAmp(iT, jj, :)./mean(data(iT, jj).Param.laser_mean_power(:));
        % Normalize for laser power
        noiseAvg(iT, jj, :) = ...
            sigNoise(iT, jj, :) * sqrt(data(iT, jj).Param.naverages * ...
            data(iT, jj).Param.nruns);
        snrAvg(iT, jj, :) = ppAmpAvg(iT, jj, :) ./ noiseAvg(iT, jj, :);
    end
end

clf
tL = tiledlayout(2, 1, "TileSpacing", "none", "Padding", "tight");
for iT = 1:nT
    nexttile()
    nrr_ = nRepRate(iT);
    for jj = 1:nrr_
        xAxis = ...
            ( ( 0:(nIntegration - 1) )*unitIntegration + 1)*(...
            data(iT, jj).x2(2) - data(iT, jj).x2(1));
        
        snrPlot = squeeze(snrAvg(iT, :, :));
        yPlot = snrPlot;
        
        if iT == 1
            if jj == 2 || jj == 6
                continue
            elseif jj == 7
                plot(xAxis, 0.43/0.65*yPlot(jj, :), 'o-', ...
                    'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            else
                plot(xAxis, yPlot(jj, :), 'o-', ...
                'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            end
        elseif iT == 2
            if jj == 1 
                continue
            elseif jj == 6
                plot(xAxis, 1.15/0.84*yPlot(jj, :), 'o-', ...
                    'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            else
                plot(xAxis, yPlot(jj, :), 'o-', ...
                    'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            end
        end
        hold on
    end
end
legend('Location', 'southeast')
labelaxesfig(tL, 'Boxcar window width / us', 'SNR')
% xlim([-0.5, max(xAxis) + 0.5])
% ylim([0, 4e3])
%exportfig()
%}

%% Boxcar integration

% Determine the range as the window in which signal > max/2
for jj = 1:nMeas
    bMM = 153;  % Index of the maximum in the field domain
    [val(jj), tMM(jj)] = max(y2{jj}(:, bMM));
    for it = tMM(jj):numel(y2{jj}(:, bMM))
        if y2{jj}(it, bMM) > val(jj)*8/10
            tStop(jj) = it;
        else
            break
        end
    end
    for it = tMM(jj):-1:1
        if y2{jj}(it, bMM) > val(jj)*8/10
            tStart(jj) = it;
        else
            break
        end
    end
end

clf
h = ScrollableAxes();
for jj = [5, 4, 3, 2, 1]
    plot(h, x{jj}{1}, x{jj}{2}, y2{jj}');
    hold on
end
legend(Param{5}.PowerAtten, Param{4}.PowerAtten, ...
    Param{3}.PowerAtten, Param{2}.PowerAtten, Param{1}.PowerAtten)
xline(x{1}{1}(tStop(2)))
xline(x{1}{1}(tStart(2)))

for jj = 1:nMeas
    y{jj} = mean(y2{jj}(tStart(jj):tStop(jj), :), 1);
end

%%
clf
for iT = 2
    for jj = 1:nRepRate(iT)
        y_ = y{iT, jj};
        plot(x{iT}{2}, datasmooth(y_, 5, "savgol"), ...
            'DisplayName', string(Param{iT, jj}.laser_rep_rate))
        hold on
%         plot(x{iT}{2}, y_{iT, jj}, ...
%             'DisplayName', string(Param{iT, jj}.laser_rep_rate))
    end
end
legend()

%% SNR ppAmp noise

for iT = 1:nT
    for jj = 1:nRepRate(iT)
        y_ = y{iT, jj}(1:end-10);
        x_ = x{iT}{2}(1:end-10);
        % y_ = yRaw{iT, jj};
        % y_ = y{iT, jj}/(Param{iT, jj}.naverages * Param{iT, jj}.nruns);
        [~, ~, noiseLev(iT, jj)] = ...
            getSNR(x_, y_, rangeNoiseB);
        % Normalize for averaging
%         SNR(iT, jj) = SNR(iT, jj) / sqrt(Param{iT, jj}.naverages * ...
%             Param{iT, jj}.nruns);
        % Normalize for laser power
%         SNR(iT, jj) = SNR(iT, jj) / mean(Param{iT, jj}.laser_mean_power(:));
        [~, ppAmp(iT, jj), ~] = ...
            getSNR(x_, datasmooth(y_, 5, "savgol"), rangeNoiseB);
        SNR(iT, jj) = ppAmp(iT, jj)/noiseLev(iT, jj);
        noiseLev(iT, jj) = noiseLev(iT, jj) * sqrt(Param{iT, jj}.naverages * ...
            Param{iT, jj}.nruns);
        SNR(iT, jj) = SNR(iT, jj) / sqrt(Param{iT, jj}.naverages * ...
            Param{iT, jj}.nruns);
    end
end
% xlim([min(x{1}{2}) max(x{1}{2})])
% legend('20dB, SNR = 220', '25 dB, SNR = 240', '30 dB, SNR = 239', ...
%     '35dB, SNR = 149', '43 dB, SNR = 68')
% exportgraphics(gcf, '/home/gianlum33/files/inbox/mwpw_comparison.png')

%
tiledlayout(3, 1)
nexttile()
for iT = 1:1
    nrr_ = nRepRate(iT);
    plot([10, 20, 30, 50, 100], noiseLev(iT, [1, 3, 4, 5, 7]), 'o-')
    hold on
end
for iT = 2
    nrr_ = nRepRate(iT);
    plot([10, 10, 20, 30, 50, 100], noiseLev(iT, 1:nrr_), 'o-')
    hold on
end
ylabel("sigma noise")

nexttile()
% [80, 80, 160, 240, 320, 720, 720];
reTuning = ppAmp(1, 1)/ppAmp(1, 2);
ppAmp(1, 7) = ppAmp(1, 7)*reTuning;
SNR(1, 7) = SNR(1, 7)*reTuning;
reTuning = ppAmp(2, 1)/ppAmp(2, 2);
ppAmp(2, 6) = ppAmp(2, 6)/reTuning;
SNR(2, 6) = SNR(2, 6)/reTuning;
for iT = 1:1
    nrr_ = nRepRate(iT);
    plot([10, 20, 30, 50, 100], ppAmp(iT, [1, 3, 4, 5, 7]), 'o-')
    hold on
end
for iT = 2
    nrr_ = nRepRate(iT);
    plot([10, 10, 20, 30, 50, 100], ppAmp(iT, 1:nrr_), 'o-')
    hold on
end
% plot((1:numel(index{2}))*10 + 50, ...
%     SNR(numel(index{1}) + 1:numel(index{1}) + numel(index{2})))
ylabel("ppAmp")

nexttile()
% [80, 80, 160, 240, 320, 720, 720];
for iT = 1:1
    nrr_ = nRepRate(iT);
    plot([10, 20, 30, 50, 100], SNR(iT, [1, 3, 4, 5, 7]), 'o-')
    hold on
end
for iT = 2
    nrr_ = nRepRate(iT);
    plot([10, 10, 20, 30, 50, 100], SNR(iT, 1:nrr_), 'o-')
    hold on
end
% plot((1:numel(index{2}))*10 + 50, ...
%     SNR(numel(index{1}) + 1:numel(index{1}) + numel(index{2})))
ylabel("SNR per average")

%% Plot snr ppAmp and noise

clf
tiledlayout(3, 1, "TileSpacing", "none", "Padding", "tight")
for iT = 2
    nrr_ = nRepRate(iT);
    for jj = 1:5
        snrPlot = squeeze(snrAvg(iT, :, :));
        noisePlot = squeeze(noiseAvg(iT, :, :));
        ppAmpPlot = squeeze(ppAmpAvg(iT, :, :));
        nexttile(1)
        if jj == 2
            yPlot = snrPlot;
            plot(xAxis, 0.84/1.15*yPlot(jj, :), 'o-', ...
            'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            ylabel('SNR')
            nexttile(2)
            yPlot = ppAmpPlot;
            plot(xAxis, 0.84*yPlot(jj, :), 'o-', ...
            'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            ylabel('ppAmp')
            nexttile(3)
            yPlot = noisePlot;
            plot(xAxis, 1.15*yPlot(jj, :), 'o-', ...
            'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            ylabel('sigma noise')
        else
            yPlot = snrPlot;
            plot(xAxis, yPlot(jj, :), 'o-', ...
            'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            ylabel('SNR')
            hold on
            nexttile(2)
            yPlot = ppAmpPlot;
            plot(xAxis, yPlot(jj, :), 'o-', ...
            'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            hold on
            ylabel('ppAmp')
            nexttile(3)
            yPlot = noisePlot;
            plot(xAxis, yPlot(jj, :), 'o-', ...
            'DisplayName', strjoin([string(data(iT, jj).laserRepRate), "Hz"]))
            hold on
            ylabel('sigma noise')        
        end
    end
end

%%
iT = 1;
irr = 7;
iboxcar = 15;

% Spin system
[Sys, Exp, Vary, VaryExp] = deal(struct());
Sys.S = [1/2 1/2];
Sys.initState = 'singlet';
Sys.g = [2.0033, 2.0024, 2.0020; 
    2.0065, 2.0053, 2.0022];
Sys.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;
Sys.lwpp = [0.6 0.];  % mT

Sys.J = unitconvert(-1e-3, 'mT->MHz');  % MHz
Sys.dip = unitconvert(+0.177,'mT->MHz'); % MHz
Sys.eeFrame = [0 90 0]*pi/180;

% Vary.g = [0.005; 0.005];
Vary.lwpp = Sys.lwpp;
% Vary.J = 1/3*abs(Sys.J);
% Vary.dip = 1/3*abs(Sys.dip);

% Experimental parameters
Exp.CenterSweep = [mean(data(iT, irr).x1) ...
    (max(data(iT, irr).x1) - min(data(iT, irr).x1))]/10 - 0.5;  % mT
Exp.mwFreq = 9.4;  % GHz
Exp.nPoints = nB;
Exp.Harmonic = 0;  % no field modulation

VaryExp.CenterSweep = [1 0];

FitOpt.BaseLine = 0;

ydata = squeeze(yMean(iT, irr, :, iboxcar));
ydata = ydata/max(ydata);
esfit(ydata, @pepper, {Sys, Exp}, {Vary, VaryExp}, FitOpt);
