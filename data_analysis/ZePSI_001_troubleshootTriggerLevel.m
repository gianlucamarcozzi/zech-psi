% Calculate SNR for different laser rep rates and for different integration
% windows and different temperatures.
clearvars

% Import
cd('/home/gianlum33/files/projects/zech_psi/data_analysis/')
expName = '../data/raw/ZePSI-E-001';
measNo = {'001', '002', '003', '004', '005'};  % T = 80 K
nMeas = numel(measNo);

for jj = 1:nMeas
    filename = [expName measNo{jj}];
    [xRaw{jj}, yRaw{jj}, yRuns{jj}, Param{jj}] ...
        = loadtreprisaak(filename);
end

[nt, nB] = size(yRaw{1});

% Data correction
freqCorr = 9.6;  % GHz
tCorrIndex = 50;  % Index of light flash
for jj = 1:nMeas
    % Laser flash at t = 0
    xCorr{jj}{1} = ...
        xRaw{jj}{1} - xRaw{jj}{1}(tCorrIndex);
    % Correct frequency to freqCorr
    xCorr{jj}{2} = ...
        xRaw{jj}{2}*(freqCorr./Param{jj}.mwfreq);
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

ttilde = 100;

%% Plot raw spectra T = 80

clf
for jj = 1:nMeas
    plot(xRaw{jj}{2}, yRaw{jj}(ttilde, :), ...
        'DisplayName', strjoin([string(Param{jj}.laser_rep_rate), "Hz"]))
    hold on
end
xlim(setaxlim(xRaw{jj}{2}, 1))
legend()

%% Plot Corr spectra

clf
for jj = 1:nMeas
    plot(x{2}, yCorr{jj}(ttilde, :), ...
        'DisplayName', strjoin([string(Param{jj}.laser_rep_rate), "Hz"]))
    hold on
end
xlim(setaxlim(x{2}, 1))
legend()

% clf
% plot(1:nB, data(2).yRuns(:, ttilde, 1))
% hold on
% plot(1:nB, data(2).yRuns(:, ttilde, 2))
% plot(1:nB, data(2).yRuns(:, ttilde, 3))
% plot(1:nB, data(2).yRaw(:, ttilde), 'LineWidth', 2)
% yyaxis right
% plot(1:nB, data.y, 'LineWidth', 2)

%% Plot to clarify error with trigger level

load("plotColors.mat")
idxOutlier = [80, 88, 89, 97, 98, 106, 123, 132, 141];
plotFactor = 1e3;

jj = 2;
yplot = yCorr{jj}*plotFactor;

ffig = figure();
clf
plot(x{2}, yplot(ttilde, :), ...
    'DisplayName', 'Exp. data', ...
    'Color', plotColors(1))
hold on
plot(x{2}(idxOutlier), yplot(ttilde, idxOutlier), ...
    'o', 'DisplayName', 'Outliers', 'Color', plotColors(2))
xlim(setaxlim(x{2}, 1))

xAxisLabel = "Magnetic field / G";
yAxisLabel = "trEPR signal / a.u.";
yticks(0)
labelaxesfig(gca, xAxisLabel, yAxisLabel)
legend()
% savefigas(gcf, '../images/troubleshoot_triggerLevel_E-001.svg')
% savefigas(gcf, '../images/troubleshoot_triggerLevel_E-001.png')
% exportgraphics(gcf, '../images/troubleshoot_triggerLevel_E-001.pdf')
