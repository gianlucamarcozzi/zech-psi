clearvars

%% IMPORT

generalFolder = "../data/raw/";
expName = "ZePSI-E-013029";
measFolder = generalFolder + expName;

[xel0, yel0, Paramel0] = loadfolderelexsys(measFolder, "*1.DTA");
[xed0, yed0, Paramed0] = loadfolderelexsys(measFolder, "*2.DTA");
[xst0, yst0, Paramst0] = loadfolderelexsys(measFolder, "*3.DTA");
[xro0, yro0, Paramro0] = loadfolderelexsys(measFolder, "*4.DTA");
[xrt0, yrt0, Paramrt0] = loadfolderelexsys(measFolder, "*5.DTA");

xAmp = table2array(readtable(append(...
    generalFolder, expName, '/', expName, '-param.txt'), "Range", "B1"));
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

nMeas = numel(Paramel0);
xeseem = xel0{1};
xdark = xed0{1};
xst = xst0{1};
% xrabi = xro0;

%% ADJUST PHASE

shiftphase = @(y, p) y*exp(1i*p*pi/180);

for ii = 1:nMeas
    [yout{ii}, bestPhase(ii)] = correctphase(yst0{ii}(1, :), 'integral');
    maxval(ii) = max(real(yout{ii}));
    imagval(ii) = trapz(abs(imag(yout{ii})));
end

%{
ifit = 1:25;  % [1:11, 15:20];
[p, S, mu] = polyfit(ifit, bestPhase(ifit), 1);
fitPhase = polyval(p, ifit, S, mu);
fitPhase2 = polyval(p, 1:nMeas, S, mu);
%}

f = figure(1);
f.Name = "Phase";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
nexttile([1, 2])
plot(1:nMeas, maxval, '.-')
yyaxis right
plot(1:nMeas, imagval, '.-')
legend("Real", "Imag", "Location", "northwest")
xlim(setaxlim(1:nMeas, 1.1))
nexttile([1, 2])
plot(1:nMeas, bestPhase, '.-')
% xlim([2, 13])
% plot(1:nMeas, fitPhase2, '.-')
% plot(ifit, fitPhase, '.-')
xlim(setaxlim(1:nMeas, 1.1))
for ii = 1:nMeas
    nexttile
    % plot(xst, real(yst0{ii}), xst, imag(yst0{ii}))
    plot(xeseem{1}, real(yed0{ii}(1, :)), xeseem{1}, imag(yed0{ii}(1, :)))
    plotTitle = strsplit(Paramel0{ii}.TITL, '-');
    title(plotTitle{end})
end

%% ADJUST PHASE

for ii = 1:nMeas
    for itau = 1:size(yel0{ii}, 1)
        yel1{ii}(itau, :) = shiftphase(yel0{ii}(itau, :), bestPhase(ii));
        yed1{ii}(itau, :) = shiftphase(yed0{ii}(itau, :), bestPhase(ii));
        yst1{ii} = shiftphase(yst0{ii}, bestPhase(ii));
    end
end

figure(1)
for ii = 1:nMeas
    nexttile(ii + 2)
    plot(xst, real(yst1{ii}), xst, imag(yst1{ii}))
    % plot(xeseem{1}, real(yed1{ii}(1, :)), xeseem{1}, imag(yed1{ii}(1, :)))
    plotTitle = strsplit(Paramel0{ii}.TITL, '-');
    title(plotTitle{end})
end

%%
II = 25;
correctphaseslider(xeseem{1}, yed0{II}(1, :), 0);

%%
% figure(4)
% clf
% sax = ScrollableAxes();
% plot(sax, xst, 1:NPH, real(yph2));
% hold on
% plot(sax, xst, 1:NPH, imag(yph2));


%%

nPlot = 5;
f = figure(2);
f.Name = "Phase shifted data";
clf
tiledlayout(nPlot, 4, "TileSpacing", "compact", "Padding", "compact")
for ii = round(linspace(1, nMeas, 5))
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yel0{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yel0{ii}));
    ylim([-10, 14]*1e3)
    plotTitle = strsplit(Paramel0{ii}.TITL, '-');
    title(append(plotTitle{end}, " BEFORE PHASE SHIFT"))

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yel1{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yel1{ii}));
    ylim([-10, 14]*1e3)
    plotTitle = strsplit(Paramel0{ii}.TITL, '-');
    title(append(plotTitle{end}, " AFTER PHASE SHIFT"))

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yed0{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yed0{ii}));
    ylim([-1.5, 3]*1e3)
    plotTitle = strsplit(Paramed0{ii}.TITL, '-');
    title(append(plotTitle{end}, " BEFORE PHASE SHIFT"))

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yed1{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yed1{ii}));
    ylim([-1.5, 3]*1e3)
    plotTitle = strsplit(Paramed0{ii}.TITL, '-');
    title(append(plotTitle{end}, " AFTER PHASE SHIFT"))
end

%% LIGHT MINUS DARK

for ii = 1:nMeas
    yeseem{ii} = yel1{ii} - yed1{ii};
end

nPlot = 9;
f = figure(2);
f.Name = "Phase shifted data";
clf
tiledlayout('flow', "TileSpacing", "compact", "Padding", "compact")
for ii = round(linspace(1, nMeas, nPlot))
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, xeseem{1}, xeseem{2}, real(yeseem{ii}));
    hold on
    plot(sax, xeseem{1}, xeseem{2}, imag(yeseem{ii}));
    ylim([-10, 14]*1e3)
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
x{2} = xeseem{2};
for ii = 1:nMeas
    Paramel0{ii}.mpfuXAmp = xAmp(ii);
end
y = yeseem;
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

