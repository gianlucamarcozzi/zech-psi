% Calculate SNR for different laser rep rates and for different integration
% windows and different temperatures.
clearvars

% Import
cd('/home/gianlum33/files/projects/zech_psi/data/raw')
expName = 'ZePSI-E-003';
measNo = {{'004', '006', '007', '008', '009', '012', '013'}, ...
    {'010', '011', '014', '015', '016', '017'}};
nT = numel(measNo);
for iT = 1:nT
    nRepRate(iT) = numel(measNo{iT});
end

for iT = 1:nT
    for jj = 1:nRepRate(iT)
        filename = [expName measNo{iT}{jj}];
        [xRaw{iT, jj}, yRaw{iT, jj}, yRuns{iT, jj}, Param{iT, jj}] ...
            = loadtreprisaak(filename);
    end
end

% Insert laser rep rate
Param{1, 3}.laser_rep_rate = 50;  % 007
Param{1, 4}.laser_rep_rate = 20;  % 008
Param{1, 5}.laser_rep_rate = 30;  % 009
Param{2, 4}.laser_rep_rate = 50;  % 015
Param{2, 5}.laser_rep_rate = 30;  % 016
Param{2, 6}.laser_rep_rate = 20;  % 017

% Sort by ascending rep rate for every temperature
for iT = 1:nT
    for jj = 1:nRepRate(iT)
        laserRepRate(iT, jj) = Param{iT, jj}.laser_rep_rate;
    end
end

for iT = 1:nT
    nrr_ = nRepRate(iT);
    inputCells = {xRaw(iT, 1:nrr_), yRaw(iT, 1:nrr_), yRuns(iT, 1:nrr_), ...
        Param(iT, 1:nrr_)};
    aa = sortcellbyparameter(laserRepRate(iT, 1:nrr_), inputCells);
    [xRaw(iT, 1:nrr_), yRaw(iT, 1:nrr_), yRuns(iT, 1:nrr_), ...
        Param(iT, 1:nrr_)] = aa{:};
end
[nt, nB] = size(yRaw{1, 1});

% Data correction
freqCorr = 9.4;  % GHz
tCorrIndex = 50;  % Index of light flash
for iT = 1:nT
    nrr_ = nRepRate(iT);
    % Order according to laser rep rate for every temperature
    for jj = 1:nrr_
        % Laser flash at t = 0
        xCorr{iT, jj}{1} = ...
            xRaw{iT, jj}{1} - xRaw{iT, jj}{1}(tCorrIndex);
        % Correct frequency to freqCorr
        xCorr{iT, jj}{2} = ...
            xRaw{iT, jj}{2}*(freqCorr./Param{iT, jj}.mwfreq);
        % Interpolate wrt the x-axis of the first measurement
        for it = 1:nt
            yCorr{iT, jj}(it, :) = interp1( ...
                xCorr{iT, jj}{2}, yRaw{iT, jj}(it, :), xCorr{iT, 1}{2}, ...
                'linear', 'extrap');
        end
    end
    x{iT} = xCorr{iT, 1};
end

% for iT = 1:nT
%     nrr_ = nRepRate(iT);
%     for jj = 1:nrr_
%         data(iT, jj).x1 = data(iT, jj).x1Corr;
%         data(iT, jj).x2 = data(iT, jj).x2Raw;
%         data(iT, jj).y = data(iT, jj).yCorr;
%     end
% end

ttilde = 100;

%% Plot spectra T = 80

clf
iT = 1;
for jj = 1:nRepRate(1)
    plot(x{iT}{2}, yRaw{iT, jj}(ttilde, :), ...
        'DisplayName', strjoin([string(Param{iT, jj}.laser_rep_rate), "Hz"]))
    hold on
end
xlim(setaxlim(x{iT}{2}, 1))
legend()
% clf
% plot(1:nB, data(2).yRuns(:, ttilde, 1))
% hold on
% plot(1:nB, data(2).yRuns(:, ttilde, 2))
% plot(1:nB, data(2).yRuns(:, ttilde, 3))
% plot(1:nB, data(2).yRaw(:, ttilde), 'LineWidth', 2)
% yyaxis right
% plot(1:nB, data.y, 'LineWidth', 2)

%% Scrollable traces

iT = 1;
jj = 5;
clf
h = ScrollableAxes();
plot(h, x{iT}{1}, x{iT}{2}, yCorr{iT, jj})
hold on
jj = 6;
plot(h, x{iT}{1}, x{iT}{2}, yCorr{iT, jj})

%% Scrollable spectra

iT = 1;
jj = 1;
clf
h = ScrollableAxes();
plot(h, x{iT}{2}', x{iT}{1}', yCorr{iT, jj}')
hold on
jj = 2;
plot(h, x{iT}{2}', x{iT}{1}', yCorr{iT, jj}')

%% Baseline correction field domain

Opt.order = 0;
rangeNoiseB = [3333, 3342; 3368, 3371];
Opt.range = rangeNoiseB;
[y1, bl1] = deal(yRaw);
for iT = 1:nT
    for jj = 1:nRepRate(iT)
        for it = 1:nt
            [y1{iT, jj}(it, :), bl1{iT, jj}(it, :)] = ...
                subtractbaseline(x{iT}{2}, yCorr{iT, jj}(it, :), Opt);
        end
    end
end

iT = 1;
jj = 6; % nRepRate(iT);
clf
h = ScrollableAxes();
plot(h, x{iT}{2}', x{iT}{1}', yCorr{iT, jj}');
hold on
plot(h, x{iT}{2}', x{iT}{1}', bl1{iT, jj}');
plot(h, x{iT}{2}', x{iT}{1}', y1{iT, jj}');
% blRegionWidth = Opt.width*(max(x{iPlot1}{2}) - min(x{iPlot1}{2}));
xline(Opt.range(1, 1));
xline(Opt.range(1, 2));
xline(Opt.range(2, 1));
xline(Opt.range(2, 2));

%% Baseline correction time domain

Opt.order = 0;
Opt.range = [-1000, 0];
[y2, bl2] = deal(yRaw);
for iT = 1:nT
    for jj = 1:nRepRate(iT)
        for ib = 1:nB
            [y2{iT, jj}(:, ib), bl2{iT, jj}(:, ib)] = ...
                subtractbaseline(x{iT}{1}, y1{iT, jj}(:, ib), Opt);
        end
    end
end

iT = 1;
jj = nRepRate(iT);
clf
h = ScrollableAxes();
plot(h, x{iT}{1}, x{iT}{2}, y1{iT, jj})
hold on
plot(h, x{iT}{1}, x{iT}{2}, bl2{iT, jj})
plot(h, x{iT}{1}, x{iT}{2}, y2{iT, jj})
xline(Opt.range(2))

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
% I use the same boxcar window for every measurement
for iT = 1:1
    for jj = 3
        bMM = 180;  % Index of the maximum in the field domain
        [val, tMM] = max(y2{iT, jj}(:, bMM));
        for it = tMM:numel(y2{iT, jj}(:, bMM))
            if y2{iT, jj}(it, bMM) > val/2
                tStop = it;
            else
                break
            end
        end
        for it = tMM:-1:1
            if y2{iT, jj}(it, bMM) > val/2
                tStart = it;
            else
                break
            end
        end
    end
end
for iT = 1:nT
    for jj = 1:nRepRate(iT)
        y{iT, jj} = mean(y2{iT, jj}(tStart:tStop, :), 1);
        yy{iT, jj} = mean(yRaw{iT, jj}(tStart:tStop, :), 1);
    end
end

iT = 1;
clf
h = ScrollableAxes();
for jj = 1:nRepRate(iT)
    plot(h, x{iT}{1}', x{iT}{2}', y2{iT, jj}');
    hold on
    legendNames(jj) = string(Param{iT, jj}.laser_rep_rate);
end
ylim([-1 2]*1e-2)
legend(legendNames)
xline(x{1}{1}(tStop))
xline(x{1}{1}(tStart))

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
