clearvars

%% IMPORT

generalFolder = "../data/raw/";
expName = "ZePSI-E-013012";
measFolder = generalFolder + expName;

[xel0, yel0, Paramel0] = loadfolderelexsys(measFolder, "*1.DTA");
[xed0, yed0, Paramed0] = loadfolderelexsys(measFolder, "*2.DTA");
[xro0, yro0, Paramro0] = loadfolderelexsys(measFolder, "*3.DTA");
[xrt0, yrt0, Paramrt0] = loadfolderelexsys(measFolder, "*4.DTA");

nMeas = numel(Paramel0);
xeseem = xel0{1};
xrabi = xro0{1};

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

%% ESEEM RINGDOWN

disp("Subtract ringdown...")

for ii = 1:nMeas
    yel1{ii} = yel0{ii} - yel0{ii}(end, :);
    yed1{ii} = yed0{ii} - yed0{ii}(end, :);
end

% Laser ESEEM
f = figure(1);
f.Name = "Laser ESEEM";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yel1{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yel1{ii}));
    ylim([-1e4, max([real(yel1{ii}), imag(yel1{ii})], [], 'all')])
    plotTitle = strsplit(Paramel0{ii}.TITL, '-');
    title(plotTitle{end})
end

% Dark ESEEM
f = figure(2);
f.Name = "Dark ESEEM";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yed1{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yed1{ii}));
    ylim([-1e4, max([real(yed1{ii}), imag(yed1{ii})], [], 'all')])
    plotTitle = strsplit(Paramed0{ii}.TITL, '-');
    title(plotTitle{end})
end

%% ADJUST PHASE

II = 25;
I_TIME_SLICE = 100;
shiftphase = @(y, p) y*exp(1i*p*pi/180);
nph = 360;
ylaser = zeros(nph, numel(xeseem{1}));
ydark = zeros(nph, numel(xeseem{1}));
ydarksmooth = zeros(nph, numel(xeseem{1}));
yrabi1 = zeros(nph, numel(xrabi{1}));
for ii = 1:nph
    ylaser(ii, :) = shiftphase(yel1{II}(I_TIME_SLICE, :), ii);
    ydark(ii, :) = shiftphase(yed1{II}(I_TIME_SLICE, :), ii);
    ydarksmooth(ii, :) = datasmooth(ydark(ii, :), 5, 'savgol');
    yrabi1(ii, :) = shiftphase(yro0{1}(I_TIME_SLICE, :), ii);
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
plotTitle = strsplit(Paramel0{II}.TITL, '-');
title(plotTitle{end})
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xeseem{1}, 1:nph, real(ydark));
hold on
plot(sax, xeseem{1}, 1:nph, imag(ydark));
plotTitle = strsplit(Paramed0{II}.TITL, '-');
title(plotTitle{end})
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xrabi{1}, 1:nph, real(yrabi1));
hold on
plot(sax, xrabi{1}, 1:nph, imag(yrabi1));
plotTitle = strsplit(Paramro0{II}.TITL, '-');
title(plotTitle{end})
nexttile
sax = ScrollableAxes('Index', 1);
plot(sax, xeseem{1}, 1:nph, real(ydarksmooth));
hold on
plot(sax, xeseem{1}, 1:nph, imag(ydarksmooth));
plotTitle = strsplit(Paramed0{II}.TITL, '-');
title(plotTitle{end} + " smooth")

f = figure(1);
f.Name = "Phase shifted data";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
BEST_PHASE = [0, 0, 0]; % [10, 28, 60];
for ii = 1:nMeas
    if ii < nMeas - 2
        phaseShift = 10;
    else
        phaseShift = BEST_PHASE(ii - nMeas + 3);
    end
    for itau = 1:size(yel1{ii}, 1)
        yel2{ii}(itau, :) = shiftphase(yel1{ii}(itau, :), phaseShift);
    end
    for itau = 1:size(yed1{ii}, 1)
        yed2{ii}(itau, :) = shiftphase(yed1{ii}(itau, :), phaseShift);
    end
    if ii > nMeas - 5
        nexttile
        sax = ScrollableAxes('Index', 1);
        plot(sax, xeseem{1}, xeseem{2}, real(yel2{ii}));
        hold on
        plot(sax, xeseem{1}, xeseem{2}, imag(yel2{ii}));
        plotTitle = strsplit(Paramel0{ii}.TITL, '-');
        title(plotTitle{end})
    end
end

%{
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
%}

%% ESEEM SUBTRACT DARK

for ii = 1:nMeas
    ye{ii} = yel2{ii} - yed2{ii};
    % ye{ii} = yel1{ii} - yed1{ii};
end

% ESEEM
figure(3)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(ye{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(ye{ii}));
    % plot(sax, xeseem{1}, xeseem{2}, real(yl2{ii}), 'Color', 'black');
    plot(sax, xeseem{1}, xeseem{2}, real(yel2{ii}), 'Color', 'black');
    ylim([-0.5e4, max([real(ye{ii}), imag(ye{ii})], [], 'all')])
    plotTitle = strsplit(Paramel0{ii}.TITL, '-');
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
for ii = 1:nMeas
    y{ii} = ye{ii}(1:end - 1, :);
end
Param = Paramel0;
save(savePath, 'x', 'y', 'Param')

% NUT1
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut1.mat";
x = xrabi;
y = yro0;
Param = Paramro0;
save(savePath, 'x', 'y', 'Param')

% NUT2
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut2.mat";
x = xrabi;
y = yrt0;
Param = Paramrt0;
save(savePath, 'x', 'y', 'Param')

