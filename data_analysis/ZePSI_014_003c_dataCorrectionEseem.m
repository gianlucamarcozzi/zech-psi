clearvars
addpath(genpath('util/'))

%% IMPORT
expName = "ZePSI-E-014003";
loadPath0 = "../data/processed/" + expName;
st = load(loadPath0 + "-Standing-ESE.mat");
el = load(loadPath0 + "-ESEEM-light.mat");
ed = load(loadPath0 + "-ESEEM-dark.mat");

%% ADJUST PHASE
% Get variables
x = st.x;
y = st.y;
param = st.param;
nMeas = numel(y);

% Find best phase
shiftphase = @(y, p) y*exp(1i*p*pi/180);
for ii = 1:nMeas
    [yout{ii}, bestPhase(ii)] = correctphase(y{ii}, 'integral');
    maxval(ii) = max(real(yout{ii}));
    imagval(ii) = trapz(abs(imag(yout{ii})));
end

ifit = 1:nMeas;  % [1:11, 15:20];
[p, S, mu] = polyfit(ifit, bestPhase(ifit), 1);
% fitPhase0 = polyval(p, ifit, S, mu);
fitPhase = polyval(p, 1:nMeas, S, mu);

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
hold on
plot(1:nMeas, fitPhase, '.-')
plot(ifit, fitPhase, '.-')
xlim(setaxlim(1:nMeas, 1.1))
for ii = 1:nMeas
    nexttile
    plot(x, real(y{ii}), x, imag(y{ii}))
    % plot(xdark, real(yout{ii}), xdark, imag(yout{ii}))
    plotTitle = strsplit(param{ii}.TITL, '-');
    title(plotTitle{end})
    xlim(setaxlim(x, 1))
end

for ii = 1:nMeas
    for itau = 1:size(el.y{ii}, 1)
        yel1{ii}(itau, :) = shiftphase(el.y{ii}(itau, :), bestPhase(ii));
        yed1{ii}(itau, :) = shiftphase(ed.y{ii}(itau, :), bestPhase(ii));
        yst1{ii} = shiftphase(y{ii}, bestPhase(ii));
    end
end

pause(5)
figure(1)
for ii = 1:nMeas
    nexttile(ii + 2)
    plot(x, real(yst1{ii}), x, imag(yst1{ii}))
    % plot(xdark, real(yout{ii}), xdark, imag(yout{ii}))
    plotTitle = strsplit(param{ii}.TITL, '-');
    title(plotTitle{end}(1:3))
end

%% LIGHT MINUS DARK
% Get variables
x = struct('x1', el.x.t, 'x2', el.x.b0);
param = el.param;

% Light minus dark (already phase corrected)
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
    plot(sax, x.x1, x.x2, real(yeseem{ii}));
    hold on
    plot(sax, x.x1, x.x2, imag(yeseem{ii}));
    ylim([-10, 14]*1e3)
    plotTitle = strsplit(param{ii}.TITL, '-');
    title(plotTitle{end})
end

nPlot = 10;
f = figure(3);
f.Name = "Phase shifted data";
clf
tiledlayout('flow', "TileSpacing", "compact", "Padding", "compact")
for ii = round(linspace(1, nMeas, nPlot))
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x.x1, x.x2, real(yeseem{ii}));
    hold on
    plot(sax, x.x1, x.x2, imag(yeseem{ii}));
    ylim([-10, 14]*1e3)
    plotTitle = strsplit(param{ii}.TITL, '-');
    title(plotTitle{end})
end

%% SAVE
pathToProcessed = "../data/processed/";

% ESEEM
clear('x', 'y', 'param')
savePath = pathToProcessed + expName + "-ESEEM.mat";
x = struct('x1', el.x.t, 'x2', el.x.b0);
y = yeseem;
param = el.param;
save(savePath, 'x', 'y', 'param')
