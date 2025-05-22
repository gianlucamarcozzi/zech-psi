clearvars
addpath(genpath('util/'))

%% IMPORT
expName = "ZePSI-E-014-007";
loadPath0 = "../data/processed/" + expName;
st{1} = load(loadPath0 + "-Standing-ESE-plus.mat");
st{2} = load(loadPath0 + "-Standing-ESE-minus.mat");
el{1} = load(loadPath0 + "-ESEEM-light-plus.mat");
el{2} = load(loadPath0 + "-ESEEM-light-minus.mat");
ed{1} = load(loadPath0 + "-ESEEM-dark-plus.mat");
ed{2} = load(loadPath0 + "-ESEEM-dark-minus.mat");

%% ADJUST PHASE
OptBl = struct('polyOrder', 1, 'range', [0, 50; 400, 600]);
for iset = 1:2
    % Get variables
    x = st{iset}.x;
    y = st{iset}.y;
    Param = st{iset}.Param;
    nMeas = numel(y);

    % Baseline correction
    for ii = 1:nMeas
        [y1{ii}, bl{ii}] = correctbaseline(x, y{ii}', OptBl);
    end
    % Find best phase
    shiftphase = @(y, p) y*exp(1i*p*pi/180);
    for ii = 1:nMeas
        [yout{ii}, bestPhase(ii)] = correctphase(y1{ii}, 'integral');
        maxval(ii) = max(abs(real(yout{ii})));
        imagval(ii) = trapz(abs(imag(yout{ii})));
    end
    
    for ii = 1:nMeas
        for itau = 1:size(el{iset}.y{ii}, 1)
            el{iset}.y2{ii}(itau, :) = shiftphase(...
                el{iset}.y{ii}(itau, :), bestPhase(ii));
            ed{iset}.y2{ii}(itau, :) = shiftphase(...
                ed{iset}.y{ii}(itau, :), bestPhase(ii));
            st{iset}.y2{ii} = shiftphase(y1{ii}, bestPhase(ii));
        end
    end

    plotstandingese(2*iset - 1, x, y1, Param, bestPhase, maxval, imagval)
    plotstandingese(2*iset, x, st{iset}.y2, Param, bestPhase, maxval, imagval)
end

plotsumdiffese(5, x, st{1}.y, st{2}.y)
plotsumdiffese(6, x, st{1}.y1, st{2}.y2)


%% LIGHT MINUS DARK
% Get variables
x = struct('x1', el{1}.x.t, ...
    'x2', el{1}.x.b0 + getparampulsespel(el{1}.Param{1}, 'd1 '));
Param = el{1}.Param;

% Light minus dark (already phase corrected)
for ii = 1:nMeas
    yeseemp{ii} = el{1}.y{ii} - ed{1}.y{ii};
    yeseemm{ii} = el{2}.y{ii} - ed{2}.y{ii};
end

f = figure(7);
f.Name = "Phase shifted data";
clf
tiledlayout('flow', "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x.x1, x.x2, real(yeseemp{ii}));
    hold on
    plot(sax, x.x1, x.x2, real(yeseemp{ii} - yeseemm{ii}));
    plot(sax, x.x1, x.x2, real(el{1}.y{ii} - el{2}.y{ii}));
    ylim([-10, 14]*1e3)
    plotTitle = strsplit(Param{ii}.TITL, '-');
    title(plotTitle{end-1})
end

f = figure(8);
f.Name = "Phase shifted data";
clf
tiledlayout('flow', "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x.x1, x.x2, real(el{1}.y{ii}));
    hold on
    plot(sax, x.x1, x.x2, real(el{2}.y{ii}));
    plot(sax, x.x1, x.x2, real(el{1}.y{ii} - el{2}.y{ii}));
    ylim([-10, 14]*1e3)
    plotTitle = strsplit(Param{ii}.TITL, '-');
    title(plotTitle{end-1})
end

%% SAVE
pathToProcessed = "../data/processed/";

% ESEEM
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-ESEEM.mat";
x = struct('x1', el.x.t, 'x2', el.x.b0);
y = yeseem;
Param = el.param;
save(savePath, 'x', 'y', 'Param')

%%
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
