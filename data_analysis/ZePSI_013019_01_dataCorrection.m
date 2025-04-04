clearvars

%% IMPORT

generalFolder = "../data/raw/";
expName = "ZePSI-E-013019";
measFolder = generalFolder + expName;

[xel0, yel0, Paramel0] = loadfolderelexsys(measFolder, "*1.DTA");
[xed0, yed0, Paramed0] = loadfolderelexsys(measFolder, "*2.DTA");
[xro0, yro0, Paramro0] = loadfolderelexsys(measFolder, "*3.DTA");
[xrt0, yrt0, Paramrt0] = loadfolderelexsys(measFolder, "*4.DTA");

nMeas = numel(Paramel0);
xeseem = xel0{1};
xdark = xed0{1};
% xrabi = xro0;

%% ADJUST PHASE

II = 25;
I_TIME_SLICE = 1;
shiftphase = @(y, p) y*exp(1i*p*pi/180);
nph = 360;
ylaser = zeros(nph, numel(xeseem{1}));
ydark = zeros(nph, numel(xdark));
% ydarksmooth = zeros(nph, numel(xeseem{1}));
% yrabi1 = zeros(nph, numel(xrabi{1}));
for ii = 1:nph
    ylaser(ii, :) = shiftphase(yel0{II}(I_TIME_SLICE, :), ii);
    ydark(ii, :) = shiftphase(yed0{II}, ii);
    % ydarksmooth(ii, :) = datasmooth(ydark(ii, :), 5, 'savgol');
    % yrabi1(ii, :) = shiftphase(yro0{1}(I_TIME_SLICE, :), ii);
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
plot(sax, xdark, 1:nph, real(ydark));
hold on
plot(sax, xdark, 1:nph, imag(ydark));
plotTitle = strsplit(Paramed0{II}.TITL, '-');
title(plotTitle{end})

figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    [yout, bestPhase(ii)] = correctphase(yed0{ii}, 'integral');
    nexttile
    plot(xdark, real(yout), xdark, imag(yout))
    maxval(ii) = max(real(yout));
end

%%

for ii = 1:nMeas
    for itau = 1:size(yel0{ii}, 1)
        yel1{ii}(itau, :) = shiftphase(yel0{ii}(itau, :), bestPhase(ii));
    end
    % for itau = 1:size(yed1{ii}, 1)
    %     yed2{ii}(itau, :) = shiftphase(yed1{ii}(itau, :), phaseShift);
    % end
end

nPlot = 5;
f = figure(2);
f.Name = "Phase shifted data";
clf
tiledlayout(nPlot, 2, "TileSpacing", "compact", "Padding", "compact")
for ii = round(linspace(1, nMeas, 5))
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yel0{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yel0{ii}));
    ylim(setaxlim([-1e4, 2e4], 1.1))
    plotTitle = strsplit(Paramel0{ii}.TITL, '-');
    title(plotTitle{end})

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yel1{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yel1{ii}));
    ylim(setaxlim([-1e4, 2e4], 1.1))
    plotTitle = strsplit(Paramel0{ii}.TITL, '-');
    title(plotTitle{end})
end

%% SAVE

% --------------------------CHANGE HERE!----------------------------------
xAmp = table2array(readtable(...
    "../data/raw/ZePSI-E-013019/ZePSI-E-013019-param.txt", "Range", "B1"));
% ------------------------------------------------------------------------
% Check for nans
isnanxAmp = isnan(xAmp);
if sum(isnanxAmp)
    for ii = 1:numel(xAmp)
        if isnanxAmp(ii) == 1
            xAmp = xAmp(1:ii - 1);
            break
        end
    end
end

disp("Save...")
pathToProcessed = "../data/processed/";

% ESEEM
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-ESEEM.mat";
x{1} = xeseem{1};
x{2} = xeseem{2};
for ii = 1:nMeas
    Paramel0{ii}.mpfuXAmp = xAmp(ii);
end
y = yel1;
Param = Paramel0;
save(savePath, 'x', 'y', 'Param')

% NUT1
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut1.mat";
x = xro0;
y = yro0;
for ii = 1:nMeas
    Paramro0{ii}.mpfuXAmp = xAmp(ii);
end
Param = Paramro0;
save(savePath, 'x', 'y', 'Param')

% NUT2
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-Nut2.mat";
x = xrt0;
y = yrt0;
Param = Paramrt0;
save(savePath, 'x', 'y', 'Param')

