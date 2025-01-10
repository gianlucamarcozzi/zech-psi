% Power dependence of the trEPR signal
clearvars

% Import
cd('/home/gianlum33/files/projects/zech_psi/data_analysis/')
expName = '../data/raw/ZePSI-E-003';
measNo = {'018', '019', '020', '021', '023'}; % 43, 25, 20, 30, 35 dB
nMeas = numel(measNo);

for ii = 1:nMeas
    filename = [expName measNo{ii} '.DTA'];
    [xRaw{ii}, yRaw{ii}, Param{ii}] = eprload(filename);
end

% Sort by ascending mw power
for ii = 1:nMeas
    mwPower(ii) = Param{ii}.MWPW;
end
aa = sortcellbyparameter(mwPower, {xRaw, yRaw, Param});
[xRaw, yRaw, Param] = aa{:};

[nt, nB] = size(yRaw{1});

% Data correction
freqCorr = 9.6e9;  % Hz
tCorrIndex = 50;  % Index of light flash
for jj = 1:nMeas
    % Laser flash at t = 0
    xCorr{jj}{1} = ...
        xRaw{jj}{1} - xRaw{jj}{1}(tCorrIndex);
    % Correct frequency to freqCorr
    xCorr{jj}{2} = ...
        xRaw{jj}{2}*(freqCorr./Param{jj}.MWFQ);
end

% Interpolate wrt the first field-axis
x = xCorr{1};
for jj = 1:nMeas
    for it = 1:nt
        yCorr{jj}(it, :) = interp1( ...
            xCorr{jj}{2}, yRaw{jj}(it, :), xCorr{1}{2});
    end
end

% Correct field-axis for different length of the frequncy-shifted axis
for jj = 1:nMeas
    minFields(jj) = min(xCorr{jj}{2});
    maxFields(jj) = max(xCorr{jj}{2});
end
lowField = max(minFields);
highField = min(maxFields);
[~, lowFieldIdx] = min(abs(x{2} - lowField));
[~, highFieldIdx] = min(abs(x{2} - highField));

x{2} = x{2}(lowFieldIdx:highFieldIdx);
for jj = 1:nMeas
    yCorr{jj} = yCorr{jj}(:, lowFieldIdx:highFieldIdx);
end

ttilde = 175;

%% Scrollable field

clf
h = ScrollableAxes();
for ii = [5, 4, 3, 2, 1]
    plot(h, xRaw{ii}{2}, x{1}, yRaw{ii}');
    hold on
end
legend(Param{5}.PowerAtten, Param{4}.PowerAtten, ...
    Param{3}.PowerAtten, Param{2}.PowerAtten, Param{1}.PowerAtten)

%% Scrollable traces

clf
h = ScrollableAxes();
for ii = [5, 4, 3, 2, 1]
    plot(h, x{1}, x{2}, yRaw{ii}');
    hold on
end
legend(Param{5}.PowerAtten, Param{4}.PowerAtten, ...
    Param{3}.PowerAtten, Param{2}.PowerAtten, Param{1}.PowerAtten)

%% Baseline correction field domain

Opt.order = 0;
Opt.width = 0.08;
[y1, bl1] = deal(yRaw);
for ii = 1:nMeas
    [nt, nB] = size(yRaw{ii});
    for it = 1:nt
        [y1{ii}(it, :), bl1{ii}(it, :)] = ...
            subtractbaseline(x{2}, yRaw{ii}(it, :)', Opt);
    end
end

iPlot1 = 1;
clf
h = ScrollableAxes();
plot(h, x{2}, x{1}, yRaw{iPlot1}');
hold on
plot(h, x{2}, x{1}, bl1{iPlot1}');
blRegionWidth = Opt.width*(max(x{2}) - min(x{2}));
xline(min(x{2}) + blRegionWidth)
xline(max(x{2}) - blRegionWidth)

%% Baseline correction time domain

Opt.order = 0;
Opt.range = [-1000, 0];
[y2, bl2] = deal(yRaw);
for ii = 1:nMeas
    [nt, nB] = size(yRaw{ii});
    for ib = 1:nB
        [y2{ii}(:, ib), bl2{ii}(:, ib)] = ...
            subtractbaseline(x{1}', y1{ii}(:, ib), Opt);
    end
end

iPlot1 = 4;
clf
h = ScrollableAxes();
plot(h, x{1}, x{2}, y1{iPlot1}');
hold on
plot(h, x{1}, x{2}, bl2{iPlot1}');
xline(Opt.range(1))
xline(Opt.range(2))

%% Scrollable field

clf
h = ScrollableAxes();
for ii = [5, 4, 3, 2, 1]
    plot(h, x{2}, x{1}, y2{ii}');
    hold on
end
legend(Param{5}.PowerAtten, Param{4}.PowerAtten, ...
    Param{3}.PowerAtten, Param{2}.PowerAtten, Param{1}.PowerAtten)

%% Scrollable traces

clf
h = ScrollableAxes();
for ii = [5, 4, 3, 2, 1]
    plot(h, x{1}, x{2}, y2{ii}');
    hold on
end
legend(Param{5}.PowerAtten, Param{4}.PowerAtten, ...
    Param{3}.PowerAtten, Param{2}.PowerAtten, Param{1}.PowerAtten)

%% Plot traces

load("plotColors.mat")
[~, btilde] = max(y2{4}(ttilde, :));
lineAtZero = zeros(1, nt);

% figure()
clf


load("plotColors.mat")

ttilde2 = round(1.7*ttilde);  % index of around 1us
btilde1 = 124;
btilde2 = 167;

ffig = figure();
ffig.Position(3) = ffig.Position(3)*1.5; 
clf
tL = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'compact');

nexttile
lineAtZero = zeros(1, nB);
xplot = x{2};
plot(xplot, lineAtZero, 'HandleVisibility', 'off', ...
    'Color', plotColors(end))
hold on
for jj = [5, 4, 3, 2, 1]
    plot(xplot, y2{jj}(ttilde, :), ...
        'DisplayName', Param{jj}.PowerAtten, 'Color', plotColors(6 - jj));
end
xline(xplot(btilde), "k-", "HandleVisibility", "off")

xlim(setaxlim(x{2}, 1))
ylim(setaxlim(y2{5}(ttilde, :), 1.05))
xticks([3390:10:3412, x{2}(end)])
xticklabels([3390, 3400, 3410, 3420])
yticks(0)
text(.05, .95, "a)",  "Units", "normalized", "FontSize", 13)
text(.05, .85, append("t = 0.5 ", char(956), "s"), ...
    "Units", "normalized", "FontSize", 13)
legend()


nexttile
lineAtZero = zeros(1, nt);
xplot = x{1}/1000;
plot(xplot, lineAtZero, 'HandleVisibility', 'off', ...
    'Color', plotColors(end))
hold on
for jj = [5, 4, 3, 2, 1]
    plot(x{1}/1000, y2{jj}(:, btilde), ...
        'DisplayName', Param{jj}.PowerAtten, 'Color', plotColors(6 - jj));
end
xline(xplot(ttilde), "k-", "HandleVisibility", "off")

xlim(setaxlim(xplot, 1))
yticks(0)
ylim(setaxlim(y2{4}(:, btilde), 1.25))

text(.35, .9425, "B_0 = 3403 G",  "Units", "normalized", ...
    "FontSize", 13)
text(.05, .95, "b)",  "Units", "normalized", "FontSize", 13)
% legend()
xAxisLabel = "Magnetic field / G";
yAxisLabel = "trEPR signal / a.u.";
labelaxesfig(tL, xAxisLabel, yAxisLabel)
ax = gca;
xlabel(append("Time / ", char(956), "s"));
ax.XLabel.FontSize = 14;

% savefigas(gcf, '../images/paramMwPower_ratioPeaksMwPower_E-003')
% savefigas(gcf, '../images/paramMwPower_ratioPeaksMwPower.png')
% exportgraphics(gcf, '../images/paramMwPower_ratioPeaksMwPower_E-003.pdf')

%% Plot different ratio between peaks at different times

load("plotColors.mat")

ttilde2 = round(1.7*ttilde);  % index of around 1us
btilde1 = 124;
btilde2 = 167;
ffig = figure();
ffig.Position(3) = ffig.Position(3)*1.5; 
clf
tL = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'compact');

nexttile
lineAtZero = zeros(1, nB);
plot(x{2}, lineAtZero, 'HandleVisibility', 'off', ...
    'Color', plotColors(end))
hold on
plot(x{2}, y2{5}(ttilde, :), "k-", ...
    'DisplayName', append("t = 0.5 ", char(956), "s"));
plot(x{2}, y2{5}(ttilde2, :), "k--", ...
    'DisplayName', append("t = 1 ", char(956), "s"));
xline(x{2}(btilde1), "k:", "HandleVisibility", "off")
xline(x{2}(btilde2), "k-.", "HandleVisibility", "off")
xlim(setaxlim(x{2}, 1))
ylim(setaxlim(y2{5}(ttilde, :), 1.05))
xticks([3390:10:3412, x{2}(end)])
xticklabels([3390, 3400, 3410, 3420])
yticks(0)
text(.05, .95, "a)",  "Units", "normalized", "FontSize", 13)
legend()


nexttile
xplot = x{1}/1000;

lineAtZero = zeros(1, nt);
plot(xplot, lineAtZero, 'HandleVisibility', 'off', ...
    'Color', plotColors(end))
hold on
plot(xplot, y2{5}(:, btilde1), "k:", ...
    'DisplayName', "B_0 = 3399.5 G");
plot(xplot, y2{5}(:, btilde2), "k-.",  ...
    'DisplayName', "B_0 = 3404.5 G");
xline(xplot(ttilde), "k-", "HandleVisibility", "off")
xline(xplot(ttilde2), "k--", "HandleVisibility", "off")
xlim([xplot(1), 5])
ylim(setaxlim([max(y2{5}(:, btilde1)), min(y2{5}(:, btilde2))], 1.05))
% xticks([3390:10:3412, x{2}(end)])
% xticklabels([3390, 3400, 3410, 3420])
yticks(0)
text(.05, .95, "b)",  "Units", "normalized", "FontSize", 13)
legend()
xAxisLabel = "Magnetic field / G";
yAxisLabel = "trEPR signal / a.u.";
labelaxesfig(tL, xAxisLabel, yAxisLabel)
ax = gca;
xlabel(append("Time / ", char(956), "s"));
ax.XLabel.FontSize = 14;

% savefigas(gcf, '../images/paramMwPower_ratioPeaksTimeSlices_E-003')
% savefigas(gcf, '../images/paramMwPower_ratioPeaksTimeSlices_E-003.png')
% exportgraphics(gcf, '../images/paramMwPower_ratioPeaksTimeSlices_E-003.pdf')

%% Boxcar integration

% Find a common index close to the maximum of each spectrum
for ii = 1:nMeas
    [~, bMaxs(ii)] = max(y2{ii}(ttilde, :));  % Index of the maximum in the field domain
end
bMM = round(mean(bMaxs));

% Integrate in the range where the signal is larger than valMax*multFactor
% for various values of multFactor between 0.5 and 0.98
multFactorX = 1000;  % multFactor factor of increase of samples points
multFactors = (50*multFactorX:98*multFactorX)/(100*multFactorX);

iMultRange = 1:numel(multFactors);
for ii = 1:nMeas
    [valMax(ii), tMM(ii)] = max(y2{ii}(:, bMM));
    for iMult = iMultRange
        % Find limit in the direction of increasing time
        for it = tMM(ii):nt
            if y2{ii}(it, bMM) > valMax(ii)*multFactors(iMult)
                tStop(iMult, ii) = it;
            else
                break
            end
        end
        % Find limit in the direction of decreasing time
        for it = tMM(ii):-1:1
            if y2{ii}(it, bMM) > valMax(ii)*multFactors(iMult)
                tStart(iMult, ii) = it;
            else
                break
            end
        end
    end
end

clf
h = ScrollableAxes();
for ii = [5, 4, 3, 2, 1]
    plot(h, x{1}, x{2}, y2{ii}');
    hold on
end
legend(Param{5}.PowerAtten, Param{4}.PowerAtten, ...
    Param{3}.PowerAtten, Param{2}.PowerAtten, Param{1}.PowerAtten)
xline(x{1}(tStop(20, 1)), 'Color', plotColors(3), ...
    'HandleVisibility', 'off')
xline(x{1}(tStart(20, 1)), 'Color', plotColors(3), ...
    'HandleVisibility', 'off')
xline(x{1}(tStop(20, 2)), 'Color', plotColors(5), ...
    'HandleVisibility', 'off')
xline(x{1}(tStart(20, 2)), 'Color', plotColors(5), ...
    'HandleVisibility', 'off')
xline(x{1}(tStop(20, 3)), 'Color', plotColors(1), ...
    'HandleVisibility', 'off')
xline(x{1}(tStart(20, 3)), 'Color', plotColors(1), ...
    'HandleVisibility', 'off')
xline(x{1}(tStop(20, 4)), 'Color', plotColors(2), ...
    'HandleVisibility', 'off')
xline(x{1}(tStart(20, 4)), 'Color', plotColors(2), ...
    'HandleVisibility', 'off')
xline(x{1}(tStop(20, 5)), 'Color', plotColors(4), ...
    'HandleVisibility', 'off')
xline(x{1}(tStart(20, 5)), 'Color', plotColors(4), ...
    'HandleVisibility', 'off')

for ii = 1:nMeas
    for iMult = iMultRange
        y{iMult, ii} = ...
            mean(y2{ii}(tStart(iMult, ii):tStop(iMult, ii), :), 1);
    end
end

%% Plot integrated spectrum for different boxcar windows
% paramMwPower_boxcarWindowSpectra_E-003

cmap = flip(viridis());
cmapLinspace = linspaceCmap(viridis, numel(iMultRange)/multFactorX + 1);

figure()
clf
iCmap = 1;
for ii = 4
    for iMult = 1:multFactorX:max(iMultRange) % [3, 20, 40] % iMultRange
        plot(x{2}, y{iMult, ii}, 'Color', cmap(cmapLinspace(iCmap), :))
        iCmap = iCmap + 1;
        hold on
    end
end

xlim(setaxlim(x{2}, 0.95))
ylim(setaxlim(y{iMult, ii}, 1.05))
yticks(0)
xAxisLabel = append("Magnetic field / G");
yAxisLabel = "trEPR signal / a.u.";
labelaxesfig(gca, xAxisLabel, yAxisLabel)
annotation("arrow", [0.575, 0.575], [0.9, 0.75])
text(.6, .875, ...
    append("Increasing", newline, "integration", newline, "window"), ...
    "Units", "normalized", "FontSize", 13)
% savefigas(gcf, '../images/paramMwPower_boxcarWindowSpectra_E-003')
% savefigas(gcf, '../images/paramMwPower_boxcarWindowSpectra_E-003.png')


%% SNR ppAmp and sigmaNoise
% sigmaNoise does not sensibly change with power, ppAmp changes

load("plotColors.mat")
% figure()
clf
noiseRangeField = [3418, 3420];
for ii = 1:nMeas
    for iMult = iMultRange % [3, 20, 40] % iMultRange
        [SNR(iMult, ii), ppAmp(iMult, ii), noiseLev(iMult, ii)] = ...
            getSNR(x{2}, y{iMult, ii}, noiseRangeField);
    end
    [snrMaxValue(ii), snrMM(ii)] = max(SNR(:, ii));
end

%% Plot SNR ppAmp sigmaNoise
% paramMwPower_boxcarWindowSNR_E-003

timeUnit = x{1}(2) - x{1}(1);
yplot = {SNR, ppAmp, noiseLev};
yplot{2} = yplot{2}*1e-4;
yLimits = {[3, 330], [0.2, 1.95], [20, 198]};
xTicks = {0.3:0.3:1.2, 0.3:0.3:1.2, 0:0.3:1.5};
xLabel = append("Boxcar window width / ", char(956), "s");
yLabel = [
    "SNR", "ppAmp / a.u.", "\sigma_{N} / a.u."];
TILEDLAYOUT_LABEL_FONT = 14;

figure();
clf
tL = tiledlayout(3, 1, 'TileSpacing', 'none', 'Padding', 'compact');
for ii = nMeas:-1:1
    xAxis{ii} = timeUnit*(tStop(:, ii) - tStart(:, ii))/1000;  % us
    for iTile = 1:3
        nexttile(iTile)
%         if iTile == 1
%             errorbar(xAxis{ii}, yplot{iTile}(:, ii), yplot{3}(:, ii), ...
%                 'o-', 'Color', plotColors(nMeas - ii + 1), ...
%                 'DisplayName', Param{ii}.PowerAtten)
%             hold on
        plot(xAxis{ii}, yplot{iTile}(:, ii), 'o-', ...
            'Color', plotColors(nMeas - ii + 1), ...
            'DisplayName', Param{ii}.PowerAtten)
        hold on
        xlim([-0.02 1.525])
        ylim(yLimits{iTile})
        xticks(xTicks{iTile})
        yLabelHandle{iTile} = ylabel(yLabel{iTile});
        ax = gca;
        ax.FontSize = TILEDLAYOUT_LABEL_FONT - 3;
        ax.YLabel.FontSize = TILEDLAYOUT_LABEL_FONT;
        ax.XLabel.FontSize = TILEDLAYOUT_LABEL_FONT;
    end
end
xlabel(xLabel)
legend('NumColumns', 3)

% % Vertical lines at the maxima of SNR for each mwPower
% % I do not keep it because the sigmaNoise is between 1 and 0.25 times the
% % SNR, therefore almost any value can be considered in the range of the
% % best SNR value
% nexttile(1)
% for ii = nMeas:-1:1
%     xline(xAxis{ii}(snrMM(ii)), 'Color', plotColors(nMeas - ii + 1), ...
%         'HandleVisibility', 'off')
% end

nexttile(2)
text(-0.08, 0.94, append(char(215), '10^4'), 'Units', 'normalized');
yLabelHandle{2}.Position(1) = yLabelHandle{1}.Position(1);
yLabelHandle{3}.Position(1) = yLabelHandle{2}.Position(1) + 0.03;

% savefigas(gcf, '../images/paramMwPower_boxcarWindowSNR_E-003')
% savefigas(gcf, '../images/paramMwPower_boxcarWindowSNR_E-003.png')
% exportgraphics(gcf, '../images/paramMwPower_boxcarWindowSNR_E-003.pdf', ...
%     'ContentType', 'vector')

%% Plot ratio between peaks vs boxcar window at different mwPower
% paramMwPower_boxcarWindowRatioPeaks_E-003

rangeLow = [3398, 3400.5];  % Low field absorption peak
lowIdx = x{2} > rangeLow(1) & x{2} < rangeLow(2);
rangeHigh = [3401.7, 3404];  % High field absorption peak
highIdx = x{2} > rangeHigh(1) & x{2} < rangeHigh(2);

figure()
clf
for ii = 1:nMeas
    for iMult = iMultRange
        peakLow(iMult, ii) = max(y{iMult, ii}(lowIdx));
        peakHigh(iMult, ii) = max(y{iMult, ii}(highIdx));
    end
end
peakRatio = peakHigh./peakLow;

for ii = nMeas:-1:1
    plot(xAxis{ii}, peakRatio(:, ii), 'o-', ...
        'Color', plotColors(nMeas - ii + 1), ...
        'DisplayName', Param{ii}.PowerAtten)
    hold on
end

xlim([-0.02 1.525])
xAxisLabel = append("Boxcar window width / ", char(956), "s");
yAxisLabel = "Ratio of amplitude of peaks";
labelaxesfig(gca, xAxisLabel, yAxisLabel)
legend('Location', 'southeast')

% savefigas(gcf, '../images/paramMwPower_boxcarWindowRatioPeaks_E-003')
% savefigas(gcf, '../images/paramMwPower_boxcarWindowRatioPeaks_E-003.png')

%% SNR ppAmp and sigmaNoise from traces
% sigmaNoise does not sensibly change with power, ppAmp changes

for ii = nMeas:-1:1
    [~, ppAmp2(ii), ~] = getSNR(x{1}, y{ii}(:, bMM), [0, 1]);
    [~, ~, noiseLev2(ii)] = getSNR(x{1}, y2{ii}(:, end), [0, 16000]);
    SNR2(ii) = ppAmp2(ii)/noiseLev2(ii);
end
% xlim([min(x{1}{2}) max(x{1}{2})])
% legend('20dB, SNR = 220', '25 dB, SNR = 240', '30 dB, SNR = 239', ...
%     '35dB, SNR = 149', '43 dB, SNR = 68')
% exportgraphics(gcf, '/home/gianlum33/files/inbox/mwpw_comparison.png')
figure()
clf
% plot([43, 35, 30, 25, 20], SNR)
hold on
plot([43, 35, 30, 25, 20], SNR2)

