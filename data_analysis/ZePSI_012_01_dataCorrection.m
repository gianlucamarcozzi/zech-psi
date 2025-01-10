%% CW EPR
clearvars
% addpath(genpath('../util'))

%% IMPORT
expName = "../data/raw/ZePSI-E-012";
% N_SCAN = 5;  % Number of scans to average

measNo = [1:2, 6];
nMeas = numel(measNo);

for ii = 1:nMeas
    filename = append(expName, sprintf('%03d', measNo(ii)));

    [x0{ii}, ytemp, Param{ii}] = eprload(filename);
    y0{ii} = ytemp{1} + 1i*ytemp{2};
    x0{ii} = x0{ii}/10;
end

% figure()
clf
for ii = 1:nMeas
    plot(x0{ii}, imag(y0{ii}))
    hold on
end

%% FREQUENCY CORRECTION
freqCorr = 9.6e9;  % Hz
for ii = 1:nMeas
    % Correct frequency to freqCorr
    x0f{ii} = ...
        x0{ii}*(freqCorr./Param{ii}.MWFQ);
end

% Interpolate wrt the first field-axis
x0i = x0f{1};
for ii = 1:nMeas
    y0i{ii} = interp1(x0f{ii}, y0{ii}, x0i);
end

% Correct field-axis for different length of the frequncy-shifted axis
for ii = 1:nMeas
    minFields(ii) = min(x0f{ii});
    maxFields(ii) = max(x0f{ii});
end
lowField = max(minFields);
highField = min(maxFields);
[~, iLowField] = min(abs(x0i - lowField));
[~, iHighField] = min(abs(x0i - highField));

idxs = iLowField:iHighField;
x = x0i(idxs);
for ii = 1:nMeas
    y00{ii} = y0i{ii}(idxs);
end

% fig = figure();
% fig.Position(4) = 1.5*fig.Position(4);
clf
tL = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for ii = 2
    nexttile(1)
    plot(x0i, real(y0i{ii}))
    hold on
    plot(x, real(y00{ii}))
    nexttile(2)
    plot(x0i, imag(y0i{ii}))
    hold on
    plot(x, imag(y00{ii}))
end

%% BASELINE CORRECTION
blOpt = struct('polyOrder', 0, 'width', 0.2);
[y1, bl1] = deal(y00);
for ii = 1:nMeas
    [y1{ii}, bl1{ii}] = correctbaseline(x, y00{ii}, blOpt);
end

fig = figure();
fig.Position(3:4) = 1.5*fig.Position(3:4);
clf
tL = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
xbl = [min(x) + blOpt.width*(max(x) - min(x)), max(x) - blOpt.width*(max(x) - min(x))];
legendInput = {"Background", "30th Oct 2024", "31st Oct 2024"};
for ii = 1:nMeas
    % Real
    nexttile
    plot(x, real(y00{ii}), "DisplayName", legendInput{ii})
    hold on
    plot(x, real(bl1{ii}), "DisplayName", "Baseline")
    xline(xbl, '--', "HandleVisibility", "off")
    xlim(setaxlim(x, 1))
    ylim(setaxlim(real(y00{ii}), 1.05))
    legend()
    % Imaginary
    nexttile
    plot(x, imag(y00{ii}), "DisplayName", legendInput{ii})
    hold on
    plot(x, imag(bl1{ii}), "DisplayName", "Baseline")
    xline(xbl, '--', "HandleVisibility", "off")
    xlim(setaxlim(x, 1))
    ylim(setaxlim(imag(y00{ii}), 1.05))
    legend()
end

%% PHASE CORRECTION
y2 = deal(y1);
for ii = 1:nMeas
    [y2{ii}, bestPhase(ii)] = correctphase(y1{ii});
end

fig = figure();
% fig.Position(3:4) = 1.5*fig.Position(3:4);
clf
tL = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% xbl = [min(x0{1}) + blOpt.width*(max(x0{1}) - min(x0{1})), max(x0{1}) - blOpt.width*(max(x0{1}) - min(x0{1}))];
legendInput = {"Background", "30th Oct 2024", "31st Oct 2024"};
for ii = 1:nMeas
    nexttile
    plot(x, real(y2{ii}), "DisplayName", append("Re ", legendInput{ii}))
    hold on
    plot(x, imag(y2{ii}), "DisplayName",  append("Imag ", legendInput{ii}))
    xlim(setaxlim(x, 1))
    ylim(setaxlim(real(y2{ii}), 1.05))
    legend()
end

%% SAVE CW EXPERIMENTS
measTitles = legendInput;
y = y2;
% save('../data/processed/ZePSI-E-012_01_cw.mat', 'x', 'y', 'Param', 'measTitles')

%% TRANSIENT EPR
clearvars

%% IMPORT
expName = "../data/raw/ZePSI-E-012";

measNo = [3, 7];
nMeas = numel(measNo);

for ii = 1:2
    if ii == 1
        measFolderName = expName + sprintf("%03d", measNo(ii));
        measFolder = dir(measFolderName + "/*.DTA");
        filename = append(measFolder(1).folder, '/', measFolder(1).name);
    else
        filename = expName + sprintf("%03d", measNo(ii));
    end
    [x0{ii}, y0{ii}, Param{ii}] = eprload(filename);
    x0{ii}{2} = x0{ii}{2}/10;
    x0{ii}{1} = x0{ii}{1}';
    x0{ii}{2} = x0{ii}{2}';
end

% Scrollable traces
% figure()
clf
sax = ScrollableAxes();
for ii = 1:nMeas
    plot(sax, x0{ii}{1}, x0{ii}{2}, y0{ii});
    hold on
end

% Scrollable field
% figure()
clf
sax = ScrollableAxes();
for ii = 1:nMeas
    plot(sax, x0{ii}{2}, x0{ii}{1}, y0{ii});
    hold on
end

%% DATA CORRECTION
[nt, nB] = size(y0{1});
iTimeCorr = 25;
freqCorr = 9.6;

for ii = 1:nMeas
    % Laser flash at t = 0
    x00{ii}{1} = x0{ii}{1} - x0{ii}{1}(iTimeCorr);
    % Correct frequency to freqCorr
    x01{ii}{2} = x00{ii}{2}*(freqCorr./Param{ii}.MWFQ);
end

% Interpolate wrt the field-axis of the first scan
x = x01{1};  % Store the x-axes in xOut
for ii = 1:nMeas
    y00{ii} = interp2( ...
        x01{ii}{2}, x01{ii}{1}, y0{ii}, x{2}, x{1});
end

% Correct field-axis for different length of the frequency-shifted axis
[minFields, maxFields] = deal(zeros(1, nMeas));
for ii = 1:nMeas
    minFields(ii) = min(x01{ii}{2});
    maxFields(ii) = max(x01{ii}{2});
end
lowField = max(minFields);
highField = min(maxFields);
[~, iLowField] = min(abs(x{2} - lowField));
[~, iHighField] = min(abs(x{2} - highField));

x{2} = x{2}(iLowField:iHighField);
for ii = 1:nMeas
    y01{ii} = y00{ii}(:, iLowField:iHighField);
end

clf
tL = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for ii = 2
    nexttile(1)
    plot(x01{ii}, real(y00{ii}))
    hold on
    plot(x, real(y01{ii}))
    nexttile(2)
    plot(x01{ii}, imag(y00{ii}))
    hold on
    plot(x, imag(y01{ii}))
end

