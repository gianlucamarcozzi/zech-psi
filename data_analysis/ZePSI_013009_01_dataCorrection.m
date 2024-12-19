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
iEseemL = [];
iEseemD = [];
iRabi1 = [];
allMeas = 1:nMeas;
iExclude = 50:62;  % Error in data acquisition
for ii = allMeas
    if size(y0{ii}) == sEseem & ~ismember(ii, iExclude) & ...
        ~ismember(ii - 1, iEseemL) & ~ismember(ii - 1, iEseemD)
        iEseemL(end + 1) = ii;
        iEseemD(end + 1) = ii + 1;
        iRabi1(end + 1) = ii + 2;
    end
end
iRabi2 = allMeas(~ismember(allMeas, [iEseemL, iEseemD, iRabi1, iExclude]));

xeseem = x0{1};
xrabi = x0{3};

inew = 0;
for ii = iEseemL
    inew = inew + 1;
    yl0{inew} = y0{ii};
    ParamL0{inew} = Param0{ii};
end
inew = 0;
for ii = iEseemD
    inew = inew + 1;
    yd0{inew} = y0{ii};
    ParamD0{inew} = Param0{ii};
end
inew = 0;
for ii = iRabi1
    inew = inew + 1;
    yrone0{inew} = y0{ii};
    ParamRone0{inew} = Param0{ii};
end
inew = 0;
for ii = iRabi2
    inew = inew + 1;
    yrtwo0{inew} = y0{ii};
    ParamRtwo0{inew} = Param0{ii};
end

nEseem = numel(yl0);
nRabi1 = numel(iRabi1);
nRabi2 = numel(iRabi2);

%% ESEEM RINGDOWN

disp("Subtract ringdown...")

for ii = 1:nEseem
    yl1{ii} = yl0{ii} - yl0{ii}(end, :);
    yd1{ii} = yd0{ii} - yd0{ii}(end, :);
end

% Laser ESEEM
f = figure(1);
f.Name = "Laser ESEEM";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nEseem
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yl1{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yl1{ii}));
    ylim([-1e4, max([real(yl1{ii}), imag(yl1{ii})], [], 'all')])
    plotTitle = strsplit(Param0{ii}.TITL, '-');
    title(plotTitle{end})
end

% Dark ESEEM
f = figure(2);
f.Name = "Dark ESEEM";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nEseem
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yd1{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yd1{ii}));
    ylim([-1e4, max([real(yd1{ii}), imag(yd1{ii})], [], 'all')])
    plotTitle = strsplit(Param0{ii}.TITL, '-');
    title(plotTitle{end})
end

%% ADJUST PHASE

I_SHIFT = 17;
I_TIME_SLICE = 100;
shiftphase = @(y, p) y*exp(1i*p*pi/180);
nph = 360;
ylaser = zeros(nph, numel(xeseem{1}));
ydark = zeros(nph, numel(xeseem{1}));
ydarksmooth = zeros(nph, numel(xeseem{1}));
yrabi1 = zeros(nph, numel(xrabi{1}));
for ii = 1:nph
    ylaser(ii, :) = shiftphase(yl1{I_SHIFT}(I_TIME_SLICE, :), ii);
    ydark(ii, :) = shiftphase(yd1{I_SHIFT}(I_TIME_SLICE, :), ii);
    ydarksmooth(ii, :) = datasmooth(ydark(ii, :), 5, 'savgol');
    yrabi1(ii, :) = shiftphase(yrone0{I_SHIFT}(I_TIME_SLICE, :), ii);
end

f = figure(3);
f.Name = "Phase shift";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xeseem{1}, 1:nph, real(ylaser));
hold on
plot(sax, xeseem{1}, 1:nph, imag(ylaser));
plotTitle = strsplit(ParamL0{I_SHIFT}.TITL, '-');
title(plotTitle{end})
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xeseem{1}, 1:nph, real(ydark));
hold on
plot(sax, xeseem{1}, 1:nph, imag(ydark));
plotTitle = strsplit(ParamD0{I_SHIFT}.TITL, '-');
title(plotTitle{end})
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xrabi{1}, 1:nph, real(yrabi1));
hold on
plot(sax, xrabi{1}, 1:nph, imag(yrabi1));
plotTitle = strsplit(ParamRone0{I_SHIFT}.TITL, '-');
title(plotTitle{end})
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xeseem{1}, 1:nph, real(ydarksmooth));
hold on
plot(sax, xeseem{1}, 1:nph, imag(ydarksmooth));
plotTitle = strsplit(ParamD0{I_SHIFT}.TITL, '-');
title(plotTitle{end} + " smooth")

f = figure(1);
f.Name = "Phase shifted data";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
BEST_PHASE = 60;
for ii = 1:nEseem
    if ii < nEseem - 4
        phaseShift = 0;
    else
        phaseShift = BEST_PHASE;
    end
    for itau = 1:size(yl1{ii}, 1)
        yl2{ii}(itau, :) = shiftphase(yl1{ii}(itau, :), phaseShift);
    end
    for itau = 1:size(yd1{ii}, 1)
        yd2{ii}(itau, :) = shiftphase(yd1{ii}(itau, :), phaseShift);
    end
    if ii > nEseem - 5
        nexttile
        sax = ScrollableAxes('Index', 1);
        plot(sax, xeseem{1}, xeseem{2}, real(yl2{ii}));
        hold on
        plot(sax, xeseem{1}, xeseem{2}, imag(yl2{ii}));
        plotTitle = strsplit(ParamL0{ii}.TITL, '-');
        title(plotTitle{end})
    end
end

BEST_PHASE_RABI = 0;
for ii = 1:nRabi1
    for itau = 1:size(y0{ii}, 1)
        yrone2{ii}(itau, :) = shiftphase(yrone0{ii}(itau, :), BEST_PHASE_RABI);
    end
end
for ii = 1:nRabi2
    for itau = 1:size(y0{ii}, 1)
        yrtwo2{ii}(itau, :) = shiftphase(yrtwo0{ii}(itau, :), BEST_PHASE_RABI);
    end
end


%% ESEEM SUBTRACT DARK

for ii = 1:nEseem
    ye{ii} = yl2{ii} - yd2{ii};
end

% ESEEM
figure(3)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = iEseemL
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(ye{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(ye{ii}));
    plot(sax, xeseem{1}, xeseem{2}, real(yl2{ii}), 'Color', 'black');
    ylim([-1e4, max([real(ye{ii}), imag(ye{ii})], [], 'all')])
    plotTitle = strsplit(ParamL0{ii}.TITL, '-');
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
for ii = 1:nEseem
    y{ii} = ye{ii}(1:end - 1, :);
end
Param = ParamL0;
% save(savePath, 'x', 'y', 'Param')

% NUT1
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut1.mat";
x = xrabi;
y = yrone2;
Param = ParamRone0;
% save(savePath, 'x', 'y', 'Param')

% NUT2
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut2.mat";
x = xrabi;
y = yrtwo2;
Param = ParamRtwo0;
% save(savePath, 'x', 'y', 'Param')

%% PLOTS
%{
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
%}
