% Power and integration window dependence of the trEPR signal
clearvars

% Import
cd('/home/gianlum33/files/projects/zech_psi/data_analysis/')
expName = '../data/raw/ZePSI-E-004';
for ii = 1:13
    if ii < 10
        measNo{ii} = strjoin(["00", string(ii)], '');
    else
        measNo{ii} = strjoin(["0", string(ii)], '');
    end
end
nMeas = numel(measNo);

for ii = 1:nMeas
    filename = [expName char(measNo{ii}) '.DTA'];
    [x0{ii}, yRaw0{ii}, Param0{ii}]  = eprload(filename);
end

% Divide 43dB measurements from 25dB measurements
index = {[], []};
for ii = 1:nMeas
    switch Param0{ii}.PowerAtten
        case '43 dB'
           index{1}(end + 1) = ii;
        case '25 dB'
           index{2}(end + 1) = ii;
        otherwise
           warning('Unexpected mw attenuation for ii = %d', ii)
    end
end

ix = 1;
for iind = 1:numel(index)
    for ii = 1:numel(index{iind})
        x{ix} = x0{index{iind}(ii)};
        yRaw{ix} = yRaw0{index{iind}(ii)};
        Param{ix} = Param0{index{iind}(ii)};
        ix = ix + 1;
    end
end

%% Baseline correction field domain

Opt.order = 0;
Opt.width = 0.1;
[y1, bl1] = deal(yRaw);
for ii = 1:nMeas
    [nt, nB] = size(yRaw{ii});
    for it = 1:nt
        [y1{ii}(it, :), bl1{ii}(it, :)] = ...
            subtractbaseline(x{ii}{2}, yRaw{ii}(it, :)', Opt);
    end
end

iPlot1 = 4;
clf
h = ScrollableAxes();
plot(h, x{iPlot1}{2}, x{iPlot1}{1}, yRaw{iPlot1}');
hold on
plot(h, x{iPlot1}{2}, x{iPlot1}{1}, bl1{iPlot1}');
blRegionWidth = Opt.width*(max(x{iPlot1}{2}) - min(x{iPlot1}{2}));
xline(min(x{iPlot1}{2}) + blRegionWidth)
xline(max(x{iPlot1}{2}) - blRegionWidth)

%% Baseline correction time domain

Opt.order = 0;
Opt.range = [0, 200];
[y2, bl2] = deal(yRaw);
for ii = 1:nMeas
    [nt, nB] = size(yRaw{ii});
    for ib = 1:nB
        [y2{ii}(:, ib), bl2{ii}(:, ib)] = ...
            subtractbaseline(x{ii}{1}', y1{ii}(:, ib), Opt);
    end
end

iPlot1 = 10;
clf
h = ScrollableAxes();
plot(h, x{iPlot1}{1}, x{iPlot1}{2}, y1{iPlot1}');
hold on
plot(h, x{iPlot1}{1}, x{iPlot1}{2}, bl2{iPlot1}');
xline(Opt.range(1))
xline(Opt.range(2))

%% SNR ppAmp and sigmaNoise from the traces
% sigmaNoise does not sensibly change with power, ppAmp changes

bMM = 2;
for ii = nMeas:-1:1
    [~, ppAmp(ii), ~] = getSNR(x{ii}{1}, y2{ii}(:, bMM), [0, 1]);
    [~, ~, noiseLev(ii)] = getSNR(x{ii}{1}, y2{ii}(:, 1), [0, 4]);
    SNR(ii) = ppAmp(ii)/noiseLev(ii);
end
% xlim([min(x{1}{2}) max(x{1}{2})])
% legend('20dB, SNR = 220', '25 dB, SNR = 240', '30 dB, SNR = 239', ...
%     '35dB, SNR = 149', '43 dB, SNR = 68')
% exportgraphics(gcf, '/home/gianlum33/files/inbox/mwpw_comparison.png')
clf
plot((1:numel(index{1}))*10 + 20, ...
    SNR(1:numel(index{1})))
hold on
plot((1:numel(index{2}))*10 + 50, ...
    SNR(numel(index{1}) + 1:numel(index{1}) + numel(index{2})))
