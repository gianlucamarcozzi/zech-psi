%
clearvars

% Import
expName = '../data/raw/ZePSI-E-008';
measNo = 1; 
nMeas = numel(measNo);

for ii = 1:nMeas
    filename = expName + sprintf("%03d", measNo(ii));
    
    [xRaw{ii}, yRaw{ii}, ParamRaw{ii}] = eprload(filename);
    xRaw{ii}{2} = xRaw{ii}{2}/10;
end

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
        xRaw{jj}{2}*(freqCorr./ParamRaw{jj}.MWFQ);
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
for ii = 1
    plot(h, xRaw{ii}{2}, x{1}, yRaw{ii}');
    hold on
end

%% Scrollable traces

clf
h = ScrollableAxes();
for ii = 1
    plot(h, x{1}, x{2}, yRaw{ii}');
    hold on
end

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

iPlot1 = 1;
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
for ii = 1
    plot(h, x{2}, x{1}, y2{ii}');
    hold on
end

%% Scrollable traces

clf
h = ScrollableAxes();
for ii = 1
    plot(h, x{1}, x{2}, y2{ii}');
    hold on
end

%% No boxcar integration for now

y = y2{1};
Param = ParamRaw{1};

%% Save

savePath = '../data/processed/ZePSI-E-008001.mat';
save(savePath, 'x', 'y', 'Param')












