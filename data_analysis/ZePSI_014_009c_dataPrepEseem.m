clearvars
addpath(genpath('util/'))

%% IMPORT
expName = "ZePSI-E-014-009";
loadPath0 = "../data/processed/" + expName;
el = load(loadPath0 + "-ESEEM-light.mat");
ed = load(loadPath0 + "-ESEEM-dark.mat");
load(loadPath0 + "-pPhases.mat");
nMeas = numel(el.y{1});

%% ADJUST PHASE
shiftphase = @(y, p) y*exp(1i*p*pi/180);  % degrees
for iset = 1:2
    for ii = 1:nMeas
        el.y1{iset}{ii} = shiftphase(el.y{iset}{ii}, pPhases(ii, iset));
        ed.y1{iset}{ii} = shiftphase(ed.y{iset}{ii}, pPhases(ii, iset));
    end
end

%% LIGHT MINUS DARK
% Get variables
x = struct('x1', el.x.t, ...
    'x2', el.x.b0 + getparampulsespel(el.Param{1}{1}, 'd1 '));
Param = el.Param;

% Light minus dark (already phase corrected)
for ii = 1:nMeas
    yeseemp{ii} = el.y{1}{ii} - ed.y{1}{ii};
    yeseemm{ii} = el.y{2}{ii} - ed.y{2}{ii};
end

fplot = "imag";
f = figure(7);
f.Name = "Phase shifted data";
clf
tiledlayout('flow', "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x.x1, x.x2, feval(fplot, yeseemp{ii}));
    hold on
    plot(sax, x.x1, x.x2, feval(fplot, yeseemm{ii}));
    plot(sax, x.x1, x.x2, feval(fplot, yeseemp{ii} - yeseemm{ii}));
    ylim([-10, 14]*1e3)
    plotTitle = strsplit(Param{1}{ii}.TITL, '-');
    title(plotTitle{end-1})
end

f = figure(8);
f.Name = "Phase shifted data";
clf
tiledlayout('flow', "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x.x1, x.x2, feval(fplot, el.y{1}{ii}));
    hold on
    plot(sax, x.x1, x.x2, feval(fplot, el.y{2}{ii}));
    plot(sax, x.x1, x.x2, feval(fplot, (el.y{1}{ii} - el.y{2}{ii})));
    ylim([-5, 10]*1e3)
    plotTitle = strsplit(Param{1}{ii}.TITL, '-');
    title(plotTitle{end-1})
end

for ii = 1:nMeas
    % yeseem{ii} = el.y{1}{ii} - el.y{2}{ii};
    yeseem{ii} = yeseemp{ii}; % - yeseemm{ii};
end

%% SAVE
pathToProcessed = "../data/processed/";

% ESEEM
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-ESEEM.mat";
x = struct('x1', el.x.t, 'x2', el.x.b0);
y = yeseem;
Param = el.Param{1};
save(savePath, 'x', 'y', 'Param')

%% FUNCTIONS
function plotstandingese(figNum, x, y, Param, bestPhase, maxval, imagval)
    nMeas = numel(bestPhase);
    f = figure(figNum);
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
    % hold on
    % plot(1:nMeas, fitPhase, '.-')
    % plot(ifit, fitPhase, '.-')
    xlim(setaxlim(1:nMeas, 1.1))
    for ii = 1:nMeas
        nexttile
        plot(x, real(y{ii}), x, imag(y{ii}))
        % plot(xdark, real(yout{ii}), xdark, imag(yout{ii}))
        plotTitle = strsplit(Param{ii}.TITL, '-');
        title(plotTitle{end-1})
        xlim(setaxlim(x, 1))
    end
end

function plotsumdiffese(figNum, x, yp, ym)
    nMeas = numel(yp);
    f = figure(figNum);
    f.Name = "Stanging-ESE: echoes sum";
    clf
    hold on
    box on
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    cmap = viridis(nMeas);
    funreim = ["real", "imag"];
    funplmi = ["plus", "minus"];
    for iplot1 = 1:2
        for iplot2 = 1:2
            nexttile(iplot1 + 2*iplot2 - 2)
            hold on
            box on
            for ii = 1:nMeas
                yplotp = feval(funreim(iplot1), yp{ii});
                yplotm = feval(funreim(iplot1), ym{ii});
                yplot = feval(funplmi(iplot2), yplotp, yplotm);
                plot(x, yplot, 'Color', cmap(ii, :))
            end
            title(funreim(iplot1) + " " + funplmi(iplot2))
            % xline([x(jMax), x(jMax - jWidth), x(jMax + jWidth)])
            xlim([100, 350])
        end
    end
end
