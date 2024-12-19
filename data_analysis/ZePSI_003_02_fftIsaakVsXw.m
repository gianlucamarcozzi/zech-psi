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
xRaw{1}{1} = xRaw{1}{1}*1000;  % ns

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

%% FFT

for jj = 1:2
    [nt, nB] = size(y2{jj});
    ntFill = nt*4;
    y2Fill = zeros(ntFill, nB);
    y2Fill(1:nt, :) = y2{jj}/max(y2{jj}, [], "all");  % Normalization
    
    dx = (x{jj}{1}(2) - x{jj}{1}(1))*1e-9;  % seconds
    df = 1/dx;  % Hz
    f{jj} = (-ntFill/2:ntFill/2 - 1)*(df/ntFill);
    
    for ib = 1:nB
        yfft{jj}(:, ib) = fft(y2Fill(:, ib));
        yfftSh{jj}(:, ib) = fftshift(yfft{jj}(:, ib));
    end
end

cmap = viridis();
clf
imagesc(f{1}, x{1}{2}, abs(yfftSh{1}))
colormap(cmap)
colorbar()


clf
tiledlayout('flow', 'TileSpacing', 'tight', 'Padding', 'compact')
nexttile()
h = ScrollableAxes();
jj = 1;
plot(h, f{jj}', 1:351, abs(yfftSh{jj}));
legend('Isaak')
nexttile()
h = ScrollableAxes();
jj = 2;
plot(h, f{jj}', 1:301, abs(yfftSh{jj}));
legend('X/W')

clf
h = ScrollableAxes();
jj = 2;
plot(h, f{jj}', 1:301, abs(yfftSh{jj}));
hold on
jj = 1;
plot(h, f{jj}', 1:351, abs(yfftSh{jj}));
legend('X/W', 'Isaak')
xlim([-2 2]*1e7)

% [bFilt, aFilt] = butter(1, 1e8/(max(f{1})/2));
% clf
% freqz(bFilt, aFilt, [], max(f{1}))
% 
% dataOut = filter(bFilt, aFilt, yfft{1}(:, 157));
% dataOutSh = fftshift(dataOut);
% clf
% tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'compact')
% nexttile()
% plot(f{1}, real(yfftSh{1}(:, 157)))
% hold on
% plot(f{1}, imag(yfftSh{1}(:, 157)))
% nexttile()
% plot(f{1}, real(dataOutSh))
% hold on
% plot(f{1}, imag(dataOutSh))

%% Plot

load("plotColors.mat")

ffig = figure(); 
clf
tL = tiledlayout(2, 1, 'TileSpacing', 'tight', 'Padding', 'tight');

nexttile
plot(f{2}, abs(yfftSh{2}(:, 1)), 'Color', plotColors(1))

nexttile
plot(f{2}, abs(yfftSh{2}(:, 1)), 'Color', plotColors(1))
hold on
plot(f{1}, abs(yfftSh{1}(:, 1)), 'Color', plotColors(2))
xlim(setaxlim(f{1}, 1))
ylim([0, 40])
xticks((-2:2)*1e7)
% xticklabels([3390, 3400, 3410, 3420])
% yticks(0)
% text(.05, .95, "a)",  "Units", "normalized", "FontSize", 13)
legend('Elexsys 4', 'Isaak')
labelaxesfig(tL, 'Frequency / Hz', 'Absolute FFT')
% savefigas(gcf, '../images/fftIsaakVsXw_offResonance_E-003')
% savefigas(gcf, '../images/fftIsaakVsXw_offResonanceIsaakVsXw_E-003.png')
% exportgraphics(gcf, '../images/fftIsaakVsXw_offResonanceIsaakVsXw_E-003.pdf')

%%
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