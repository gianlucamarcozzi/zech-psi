clearvars

%% IMPORT

generalFolder = "../data/raw/";
expName = "ZePSI-E-013009";
measFolder = generalFolder + expName;

[x0, y0, Param0] = loadfolderelexsys(measFolder);
nMeas = numel(Param0);

%% Scrollables
%{
disp("Plot...")
% Traces real
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x0{ii}{1}, x0{ii}{2}, real(y0{ii}));
    hold on
end

% Traces imag
figure(2)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes();
    plot(sax, x0{ii}{1}, x0{ii}{2}, imag(y0{ii}));
    hold on
end
%}
%% SEPARATE ESEEM FROM NUTATIONS

sEseem = size(y0{1});
iEseem = [];
iRabi1 = [];
allMeas = 1:nMeas;
iExclude = 50:62;  % Error in data acquisition
for ii = allMeas
    if size(y0{ii}) == sEseem & ...
        ~ismember(ii - 1, iEseem) & ~ismember(ii, iExclude)
        iEseem(end + 1) = ii;
        iEseem(end + 1) = ii + 1;
        iRabi1(end + 1) = ii + 2;
    end
end
iRabi2 = allMeas(~ismember(allMeas, [iEseem, iRabi1, iExclude]));

xeseem = x0{1};
xrabi = x0{3};

%% ESEEM RINGDOWN

disp("Subtract ringdown...")

for ii = iEseem
    y1{ii} = y0{ii} - y0{ii}(end, :);
end

% Laser ESEEM
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iEseem(1:2:end)
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(y1{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(y1{ii}));
    ylim([-1e4, max([real(y1{ii}), imag(y1{ii})], [], 'all')])
    plotTitle = strsplit(Param0{ii}.TITL, '-');
    title(plotTitle{end})
end

% Dark ESEEM
figure(2)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iEseem(2:2:end)
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(y1{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(y1{ii}));
    ylim([-1e4, max([real(y1{ii}), imag(y1{ii})], [], 'all')])
    plotTitle = strsplit(Param0{ii}.TITL, '-');
    title(plotTitle{end})
end

%% ADJUST PHASE

I_Y1_LASER = nMeas - 6;
shiftphase = @(y, p) y*exp(1i*p*pi/180);
nph = 360;
ylaser = zeros(nph, numel(xeseem{1}));
ydark = zeros(nph, numel(xeseem{1}));
yrabi1 = zeros(nph, numel(xrabi{1}));
for ii = 1:nph
    ylaser(ii, :) = shiftphase(y1{I_Y1_LASER}(150, :), ii);
    ydark(ii, :) = shiftphase(y1{I_Y1_LASER + 1}(150, :), ii);
    yrabi1(ii, :) = shiftphase(y0{I_Y1_LASER + 2}(1, :), ii);
end

% figure()
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xeseem{1}, 1:nph, real(ylaser));
hold on
plot(sax, xeseem{1}, 1:nph, imag(ylaser));
plotTitle = strsplit(Param0{I_Y1_LASER}.TITL, '-');
title(plotTitle{end})
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xeseem{1}, 1:nph, real(ydark));
hold on
plot(sax, xeseem{1}, 1:nph, imag(ydark));
plotTitle = strsplit(Param0{I_Y1_LASER + 1}.TITL, '-');
title(plotTitle{end})
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xrabi{1}, 1:nph, real(yrabi1));
hold on
plot(sax, xrabi{1}, 1:nph, imag(yrabi1));
plotTitle = strsplit(Param0{I_Y1_LASER + 2}.TITL, '-');
title(plotTitle{end})


figure()
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
BEST_PHASE = 63;
for ii = iEseem
    if ii < numel(iEseem) - 7
        phaseShift = 0;
    else
        phaseShift = BEST_PHASE;
    end
    for itau = 1:size(y1{ii}, 1)
        y2{ii}(itau, :) = shiftphase(y1{ii}(itau, :), phaseShift);
    end
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(y2{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(y2{ii}));
end

BEST_PHASE_RABI = 0;
for ii = [iRabi1, iRabi2]
    for itau = 1:size(y0{ii}, 1)
        y2{ii}(itau, :) = shiftphase(y0{ii}(itau, :), BEST_PHASE_RABI);
    end
end

%% ESEEM SUBTRACT DARK

iEseemL = iEseem(1:2:numel(iEseem));
for ii = iEseemL
    y3{ii} = y2{ii} - y2{ii + 1};
end

% ESEEM
figure(3)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iEseemL
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(y3{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(y3{ii}));
    ylim([-1e4, max([real(y3{ii}), imag(y3{ii})], [], 'all')])
    plotTitle = strsplit(Param0{ii}.TITL, '-');
    title(plotTitle{end})
end


%% SAVE

disp("Save...")
pathToProcessed = "../data/processed/";

% ESEEM
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-ESEEM.mat";
x{1} = xeseem{1};
x{2} = xeseem{2}(1:end - 1);
% Cut out the last trace that is zero everywhere
for ii = 1:numel(iEseemL)
    y{ii} = y3{iEseemL(ii)}(1:end - 1, :);
end
Param = Param0(iEseem);
% save(savePath, 'x', 'y', 'Param')

% NUT1
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut1.mat";
x = xrabi;
y = y0(iRabi1);
Param = Param0(iRabi1);
% save(savePath, 'x', 'y', 'Param')

% NUT2
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut2.mat";
x = xrabi;
y = y0(iRabi2);
Param = Param0(iRabi2);
% save(savePath, 'x', 'y', 'Param')

%% PLOTS

iy = iEseem(14);
nSpectrum = 5;
i1 = round(linspace(1, size(y0{iy}, 1), nSpectrum));
i2 = 30:size(y0{iy}, 2);
% cmap = viridis(nSpectrum);
% figure()
clf
tL = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
yplot = real(y0{iy})*1e-4;
for ii = 1:nSpectrum
    plot(x1{1}(i2), yplot(i1(ii), i2),  'DisplayName', append(...
        char(hex2dec('03C4')), sprintf(' = %d ns', 150 + i1(ii) - 1)))
    hold on
end
xlim(setaxlim(x1{1}(i2), 1))
ylim(setaxlim(yplot(i1(2), i2), 1.05))
xticks(100:200:1000)
yticks((-5:2:5))
legend()
nexttile
yplot = imag(y0{iy})*1e-4;
for ii = 1:nSpectrum
    plot(x1{1}(i2), yplot(i1(ii), i2))  % + (ii - 1)*0.6e4
    hold on
end
xlim(setaxlim(x1{1}(i2), 1))
ylim(setaxlim(yplot(i1(2), i2), 1.05))
xticks(100:200:1000)
yticks((-5:0))
labelaxesfig(tL, "Time / ns",...
    ["In-phase / a.u.", "Out-of-phase / a.u."])
% nexttile(1); ylabel("In-phase channel / a.u.")
% nexttile(2); ylabel("Out-of-phase channel / a.u.")
savePath = "../images/" + "ZePSI-013004-01-eseemRaw.png";
% saveas(gcf, savePath)
