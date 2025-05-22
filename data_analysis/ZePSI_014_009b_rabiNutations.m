%%
clearvars
addpath(genpath('util/'))

%% PULSE 2
loadPath = '../data/processed/ZePSI-E-014-009-Nut2.mat';
pul2 = load(loadPath);
nMeas = numel(pul2.x);

% Get variables
x = pul2.x;
y = pul2.y;
Param = pul2.Param;

% Variables relative to fitting
expcos = @(xx, p) [exp(-xx/p(1)).*cos(2*pi*p(2)*xx), ones(numel(xx), 1)];
fitOpt = optimoptions('lsqnonlin','Display','off');
p0 = [40, 30e-3];

for ii = 1:nMeas
    % Data to fit
    ydata = real(y{ii});
    ydata = ydata/max(ydata);
    xdata = x{ii};
    % Define fitmodel using current xdata
    fitmodel = @(p) expcos(xdata, p);

    % Fit
    [yfit{ii}, pfit{ii}, pci{ii}] = lsqnonlin2steps(...
        ydata, fitmodel, p0, fitOpt);

    % Store important fit parameters
    freqRabi(ii) = pfit{ii}(2)*1e3;  % MHz
    piPerf(ii) = 1/2/freqRabi(ii)*1e3;  % ns
    dfreqRabi(ii) = pci{ii}(2)*1e3;  % MHz
end

% -------------------------- PLOT ----------------------------------------
f = figure(1);
f.Name = "Rabi nutations second pulse";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    % Call variables
    ydata = real(y{ii});
    ydata = ydata/max(ydata);
    xdata = x{ii};
    
    % Plot
    nexttile
    plot(xdata, ydata, 'o')
    hold on
    plot(xdata, yfit{ii})
    title(sprintf('%d: %.2f MHz, %.2f ns', ii, freqRabi(ii), piPerf(ii)))
end
% Plot frequencies with uncertainties
nexttile
errorbar(1:nMeas, freqRabi, dfreqRabi, 'o')
xlim([0.9, nMeas + 0.1])
freqRabiAvg = sum(freqRabi./dfreqRabi)/sum(1./dfreqRabi);  % MHz
piPerfAvg = 1/2/freqRabiAvg*1e3;  % ns
title(sprintf("Avg: %.2f MHz, %.2f ns", freqRabiAvg, piPerfAvg))

%% PULSE 1
loadPath = strrep(loadPath, 'Nut2', 'Nut1');
pul1 = load(loadPath);
endings = ["plus", "minus"];

f = figure(2);
f.Name = "Rabi nutations first pulse";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for iph = 1:2  % Phase of the microwave
    % Get variables
    x = pul1.x;
    y = pul1.y{iph};
    Param = pul1.Param{iph};
    nMeas = numel(y);
    
    % Variables relative to fitting
    expcos = @(xx, p) [exp(-xx/p(1)).*cos(2*pi*p(2)*xx), ones(numel(xx), 1)];
    fitOpt = optimoptions('lsqnonlin','Display','off');
    % w0s = linspace(0.5, sqrt(30), nMeas).^2*1e-3;
    w0s = linspace(0.5, 24, nMeas)*1e-3;
    % w0s = flip(w0s);
    p01 = 170;
    
    %------------------------------ FIT --------------------------------------
    for ii = 7:nMeas
        % Data to fit
        ydata = real(y{ii});
        ydata = ydata/max(ydata);
        xdata = x{ii};
        % Current initial parameters
        p0 = [p01, w0s(ii)];
        % Define fitmodel using current xdata
        fitmodel = @(p) expcos(xdata, p);
    
        % Fit
        [yfit{ii}, pfit{ii}, pci{ii}] = lsqnonlin2steps(...
            ydata, fitmodel, p0, fitOpt);
    
        % Store important fit parameters
        freqRabi(ii) = pfit{ii}(2)*1e3;  % MHz
        piPerf(ii) = 1/2/freqRabi(ii)*1e3;  % ns
        dfreqRabi(ii) = pci{ii}(2)*1e3;  % MHz
    end
    
    % -------------------------- PLOT ----------------------------------------
    for ii = 7:nMeas
        % Call variables
        ydata = real(y{ii});
        ydata = ydata/max(ydata);
        xdata = x{ii}';
        
        % Plot
        nexttile(ii)
        hold on; box on;
        plot(xdata, ydata, 'o')
        plot(xdata, yfit{ii})
        title(sprintf('%d: %.2f MHz, %.2f ns', ii, freqRabi(ii), piPerf(ii)))
    end
    % Plot frequencies with uncertainties
    for ii = 7:nMeas
        pAmp(ii) = Param{ii}.pAmpPulse1;
    end
    nexttile(ii + 1)
    hold on; box on;
    errorbar(pAmp(7:end), freqRabi(7:end), dfreqRabi(7:end), 'o')
    xlim(setaxlim(pAmp(7:end), 1.05))
end

% Taken from misc_E_002_003_mpfuBrX.m for the ones where the fit is bad
freqExpected = [2.7219, 3.5798, 4.4428, 5.2997, 6.1490, 7.0064, 7.8794];
freqRabiFin = [freqExpected, freqRabi(8:end)];
piPerfFin = [1/2./freqExpected*1e3, piPerf(8:end)];

%% ------------------------ SAVE INTO ESEEM STRUCTURES ---------------------
%{
loadEseemPath = append(loadPath(1:end-9), '-ESEEM-light.mat');
loadEseem = load(loadEseemPath);
turningAngle = zeros(nMeas, 1);
pulseLength = getparampulsespel(loadEseem.Param{1}, 'p1 ');
for ii = 1:nMeas
    turningAngle(ii) = pulseLength/piPerfFin(ii)*pi*23.11/24; % *(pulseLength/RABI_LENGTH_PI_PULSE);
    loadEseem.Param{ii}.turningAngle = turningAngle(ii);
end
save(loadEseemPath, '-struct', 'loadEseem')
%}

%% RMTM
generalFolder = "../data/raw/";
expName = "ZePSI-E-014-009";
expFolder = generalFolder + expName;

% Pulse amplitude of first pulse
filepathParam = append(expFolder, '/', expName, '-param.txt');
pAmps = getamplitudefirstpulse(filepathParam);

% BrX and BrMinX standing echo data (p plus, m minus)
[yrm1, ytm1, yrm2, ytm2, xsc] = importrmtm(expFolder, expName);
% Number of measurements
nMeas = size(yrm1, 1);

%% OSCILLOSCOPE
N_PULSE = 2;

% Define range of the x-axis where the pulses are expected
XRANGETM{1} = [3.05, 3.20]*1e-7;
XRANGETM{2} = [7.05, 7.20]*1e-7;
XRANGERM{1} = [3.33, 3.51]*1e-7;
XRANGERM{2} = [7.33, 7.51]*1e-7;
% Range of the baseline (outside)
XRANGEBL = [4.5, 6]*1e-7;

% Separate pulse 1 and pulse 2 from the y1 and y2 arrays
for ip = 1:N_PULSE
    xCondition = xsc > XRANGETM{ip}(1) & ...
        xsc < XRANGETM{ip}(2);
    ytm1cut{ip} = ytm1(:, xCondition);
    ytm2cut{ip} = ytm2(:, xCondition);
    xtm{ip} = xsc(xCondition);
    xCondition = xsc > XRANGERM{ip}(1) & ...
        xsc < XRANGERM{ip}(2);
    yrm1cut{ip} = yrm1(:, xCondition);
    yrm2cut{ip} = yrm2(:, xCondition);
    xrm{ip} = xsc(xCondition);
end
xConditionBl = XRANGEBL(1) < xsc & xsc < XRANGEBL(2);
bltm1 = mean(ytm1(:, xConditionBl), 2);
bltm2 = mean(ytm2(:, xConditionBl), 2);
blrm1 = mean(yrm1(:, xConditionBl), 2);
blrm2 = mean(yrm2(:, xConditionBl), 2);

% Average
for ip = 1:N_PULSE
    meantm1{ip} = mean(ytm1cut{ip}, 2) - bltm1;
    meantm2{ip} = mean(ytm2cut{ip}, 2) - bltm2;
    amptm1{ip} = abs(meantm1{ip});
    amptm2{ip} = abs(meantm2{ip});
    phtm1{ip} = angle(meantm1{ip})*180/pi;
    phtm2{ip} = angle(meantm2{ip})*180/pi;
    % Small sign correction
    phtm2{ip}(phtm2{ip} < -150) = ...
        phtm2{ip}(phtm2{ip} < -150) + 360;
end
% Average
for ip = 1:N_PULSE
    meanrm1{ip} = mean(yrm1cut{ip}, 2) - blrm1;
    meanrm2{ip} = mean(yrm2cut{ip}, 2) - blrm2;
    amprm1{ip} = abs(meanrm1{ip});
    amprm2{ip} = abs(meanrm2{ip});
    phrm1{ip} = angle(meanrm1{ip})*180/pi;
    phrm2{ip} = angle(meanrm2{ip})*180/pi;
    % Small sign correction
    phrm2{ip}(phrm2{ip} < -150) = ...
        phrm2{ip}(phrm2{ip} < -150) + 360;
end
idxPlots = 1:4;
plotscopedata(1, "RM BrX", xsc, yrm1, xrm, yrm1cut, blrm1, idxPlots)
plotscopedata(3, "TM BrX", xsc, ytm1, xtm, ytm1cut, bltm1, idxPlots)
plotscopedata(5, "RM BrMinX", xsc, yrm2, xrm, yrm2cut, blrm2, idxPlots)
plotscopedata(7, "TM BrMinX", xsc, ytm2, xtm, ytm2cut, bltm2, idxPlots)
%%
plotscopedatacmp(9, "RM plus", xsc, (yrm2 - blrm2), xsc, yrm1 - blrm1, idxPlots)
plotscopedatacmp(11, "TM plus", xsc, (ytm2 - bltm2), xsc, ytm1 - bltm1, idxPlots)
plotscopedatacmp(13, "RM minus", xsc, -(yrm2 - blrm2), xsc, yrm1 - blrm1, idxPlots)
plotscopedatacmp(15, "TM minus", xsc, -(ytm2 - bltm2), xsc, ytm1 - bltm1, idxPlots)

%% SAVE IMAGE
savepathpdf = '../images/014-009_cmpRMTM.pdf';
delete(savepathpdf);
for ii = [15, 16, 13, 14]
    exportgraphics(figure(ii), savepathpdf, "Append", true)
end

%% SAVE PHASES
pPhases = [phtm1{1}, phtm2{1} - 180];
save("../data/processed/ZePSI-E-014-009-pPhases.mat", "pPhases", "-mat")

%%
plotscopefinal(17, pAmps, ...
    meantm1{1}, meantm2{1}, meantm1{2}, ...
    amptm1{1}, amptm2{1}, amptm1{2}, ...
    phtm1{1}, phtm2{1}, phtm1{2}, [0, 180])
plotscopefinal(18, pAmps, ...
    meanrm1{1}, meanrm2{1}, meanrm1{2}, ...
    amprm1{1}, amprm2{1}, amprm1{2}, ...
    phrm1{1}, phrm2{1}, phrm1{2}, [90, -90])

%% ------------------------- STANDING ESE --------------------------------
% ADJUST PHASE
% Get variables
loadPath0 = "../data/processed/" + expName;
st = load(loadPath0 + "-Standing-ESE.mat");
OptBl = struct('polyOrder', 1, 'width', 0.15);
for iset = 1:2
    x = st.x;
    y = st.y{iset};
    Param = st.Param{iset};
    nMeas = numel(y);

    % Baseline correction
    for ii = 1:nMeas
        [y1{ii}, bl{ii}] = correctbaseline(x, y{ii}, OptBl);
    end
    % Find best phase
    shiftphase = @(y, p) y*exp(1i*p*pi/180);  % degrees
    idxs = 100:300;
    for ii = 1:nMeas
        [yout{ii}, bestPhase{iset}(ii)] = correctphase(y1{ii}, 'pulse');
    end
    
    % SOME CORRECTIONS!!
    if iset == 2
        for ii = 1:nMeas
            if bestPhase{iset}(ii) < 0
                bestPhase{iset}(ii) = bestPhase{iset}(ii) + 360;
            end
        end
        % bestPhase{iset}(1:4) = 150;
    end

    for ii = 1:nMeas
        for itau = 1:size(st.y{iset}{ii}, 1)
            st.y2{iset}{ii} = shiftphase(y1{ii}, bestPhase{iset}(ii));
            intgs(ii) = trapz(st.y2{iset}{ii});
        end
    end

    plotstandingese(20 + 2*iset - 1, x, y1, Param, bestPhase{iset}, intgs)
    plotstandingese(20 + 2*iset, x, st.y2{iset}, Param, bestPhase{iset}, intgs)
end
%%
% plotsumdiffese(5, x, st{1}.y, st{2}.y)
% plotsumdiffese(6, x, st{1}.y1, st{2}.y2)

ypbl = st.y2{1};
ymbl = st.y2{2};
% Integration - Get amplitude and phase
% [valMax, jMax] = max(real(ypbl{1}));
% [~, jWidth] = min(abs(real(ypbl{1}(jMax:end)) - valMax/10));
[valMax, jMax] = max(real(ypbl{1}));
[~, jWidth] = min(abs(real(ypbl{1}(jMax:end)) - valMax/10));
for imeas = 1:nMeas
    yp1(imeas) = trapz(ypbl{imeas});
    ym1(imeas) = trapz(ymbl{imeas});
end
ypamp = abs(yp1);
ypph = angle(yp1)*180/pi;
ymamp = abs(ym1);
ymph = angle(ym1)*180/pi;
ymph(ymph > -150) = ymph(ymph > -150) - 360;

% figure(14)
% nexttile(2)
% hold on
% plot(pulseAmps, asin(ypamp/max(ypamp))/pi*2*0.0476, '--')
% plot(pulseAmps, asin(ymamp/max(ypamp))/pi*2*0.0476, '-.')

% plotblcorrection(11, st{1}.x, yp, ym, blp, blm);
plotsumdiffese(15, st.x, ypbl, ymbl);
plotstandingese(16, pAmps, ypamp, ypph, ymamp, ymph);

%% ------------------------- SAVE TO TXT ---------------------------------
% ------------------------------------------------------------------------
IEXP = 2;
chname = ["mpfu_BrX", "mpfu_BrMinX"];
for ip = 1:2
    savepath = chname(ip) + "_ampPhase_2025-03-28.txt";
    savemat = [mpfuAmps; amplitude{ip}'; phase{ip}'];
    % fid = fopen(savepath, 'w');
    % fprintf(fid, "%f %f %f\n", savemat);
    % fclose(fid);
end

%%

%% ------------------------- FUNCTIONS -----------------------------------
function pAmps = getamplitudefirstpulse(filepathParam)
    % Load the amplitude of the mpfu of the first pulse
    pAmps = table2array(readtable(filepathParam, "Range", "B1"));
    % Check for nans and cut the array as soon as one is found
    isnanxAmp = isnan(pAmps);
    if sum(isnanxAmp)
        for ii = 1:numel(pAmps)
            if isnanxAmp(ii) == 1
                pAmps = pAmps(1:ii - 1);
                break
            end
        end
    end
end

function [yrm1, ytm1, yrm2, ytm2, x] = importrmtm(expFolder, expName)
    pathpart1 = append(expFolder, "/", expName);
    % Pulse 1
    ar = importdata(pathpart1 + "-scope-p1-rm1.txt");
    ai = importdata(pathpart1 + "-scope-p1-rm2.txt");
    br = importdata(pathpart1 + "-scope-p1-tm1.txt");
    bi = importdata(pathpart1 + "-scope-p1-tm2.txt");
    % y values pulse 1
    yrm1 = ar(2:end, :) + 1i*ai(2:end, :);
    ytm1 = br(2:end, :) + 1i*bi(2:end, :);
    % Pulse 2
    ar = importdata(pathpart1 + "-scope-p2-rm1.txt");
    ai = importdata(pathpart1 + "-scope-p2-rm2.txt");
    br = importdata(pathpart1 + "-scope-p2-tm1.txt");
    bi = importdata(pathpart1 + "-scope-p2-tm2.txt");
    % y values pulse 2
    yrm2 = ar(2:end, :) + 1i*ai(2:end, :);
    ytm2 = br(2:end, :) + 1i*bi(2:end, :);
    % x values
    x = ar(1, :);
end

function plotblcorrection(figNum, x, yp, ym, blp, blm)
    nMeas = numel(yp);
    f = figure(figNum);
    f.Name = "Stanging-ESE: echoes";
    clf
    hold on
    box on
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    cmap = viridis(nMeas);
    funplots = ["real", "imag"];
    for iplot = 1:2
        nexttile(iplot)
        hold on
        box on
        for ii = 1:nMeas
            yplot = feval(funplots(iplot), yp{ii});
            plot(x, yplot, 'Color', cmap(ii, :))
            yplot = feval(funplots(iplot), ym{ii});
            plot(x, yplot, 'Color', cmap(ii, :))
        end
        title(funplots(iplot))
        % xline([x(jMax), x(jMax - jWidth), x(jMax + jWidth)])
        xlim([100, 350])
    end
    for ii = nMeas-6:nMeas
        nexttile(ii - nMeas + 9)
        hold on
        box on
        yplotr = real(yp{ii});
        yploti = imag(yp{ii});
        plot(x, yplotr)
        plot(x, yploti)
        yplotr = real(ym{ii});
        yploti = imag(ym{ii});
        plot(x, yplotr)
        plot(x, yploti)
        yplotblp = imag(blp{ii});
        yplotblm = imag(blm{ii});
        plot(x, yplotblp, 'r--')
        plot(x, yplotblm, 'r--')
        yplotblp = real(blp{ii});
        yplotblm = real(blm{ii});
        plot(x, yplotblp, 'r')
        plot(x, yplotblm, 'r')
        xlim(setaxlim(x))
    end
end

function plotscopedata(figNum, figNameStart, x, y, xcut, ycut, bl, idxPlots)
    funplots = ["real", "imag"];
    figNames = figNameStart + funplots;
    for ifig = 1:2
        f = figure(figNum + (ifig - 1));  % 2 to 4
        f.Name = figNames(ifig);
        funplot = funplots(ifig);
        % clf
        tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
        for iplot = idxPlots
            nexttile(iplot)
            hold on; box on;
            xplot = x;
            yplot = y(iplot, :) - bl(iplot, :);
            yplot = feval(funplot, yplot);
            plot(xplot, yplot)
            xlim(setaxlim(xplot, 1))
            hold on
            for ip = 1:numel(ycut)
                yplot = ycut{ip}(iplot, :) - bl(iplot, :);
                xplot = xcut{ip};
                yplot = feval(funplot, yplot);
                plot(xplot, yplot)
            end
        end
    end
end

function plotscopedatacmp(figNum, figNameStart, x, y, x2, y2, idxPlots)
    funplots = ["real", "imag"];
    figNames = figNameStart + " " + funplots;
    for ifig = 1:2
        f = figure(figNum + (ifig - 1));  % 2 to 4
        f.Name = figNames(ifig);
        funplot = funplots(ifig);
        % clf
        tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
        for iplot = idxPlots
            nexttile(iplot)
            hold on; box on;
            xplot = x;
            yplot = y(iplot, :);
            yplot = feval(funplot, yplot);
            plot(xplot, yplot)
            xlim(setaxlim(xplot, 1))
            xplot = x2;
            yplot = y2(iplot, :);
            yplot = feval(funplot, yplot);
            plot(xplot, yplot)            
            xlim(setaxlim(xplot, 1))
        end
    end
end

function plotscopefinal(figNum, pulseAmps, y1, y2, yX, ...
    amps1, amps2, ampsX, phases1, phases2, phasesX, expectPhases)
    figure(figNum)
    axTitles = ["TM real", "TM imag"];
    funplots = split(axTitles(:), ' ');
    funplots = funplots(:, 2);  % real, imag
    clf
    tiledlayout(2, 3, "TileSpacing", "compact", "Padding", "compact")
    % Plots 1 and 3: real and imag part
    numTile = [1, 4];
    for itile = 1:2
        nexttile(numTile(itile))
        hold on
        box on
    
        xplot = pulseAmps;
        plot(xplot, feval(funplots(itile), y1), 'o-')  % <+x>
        plot(xplot, feval(funplots(itile), y2), 'o-')  % <-x>
        plot(xplot, feval(funplots(itile), yX), 'o-')  % +x pi
    
        title(axTitles(itile))
        xlim(setaxlim(xplot, 1))
    end
    
    % Amplitude
    nexttile(2)
    hold on
    box on
    plot(xplot, amps1, 'o-')  % <+x>
    plot(xplot, amps2, 'o-')  % <-x>
    plot(xplot, ampsX, 'o-')  % +x pi
    title('Amplitude')
    xlim(setaxlim(xplot, 1))
    
    % Phase
    nexttile(3)
    hold on
    box on
    plot(xplot, phases1, 'o-')  % <+x>
    plot(xplot, phases2, 'o-')  % <-x>
    plot(xplot, phasesX, 'o-')  % +x pi
    title('phase')
    xlim(setaxlim(xplot, 1))
    
    % Relative amplitude difference
    nexttile(5)
    hold on
    box on
    yplot = (amps1 - amps2)./amps1;
    plot(xplot, yplot, 'ko-')
    title('Relative amp diff')
    xlim(setaxlim(xplot, 1))
    
    % Phase difference from expected
    nexttile(6)
    hold on
    box on
    yplot = (phases1 - expectPhases(1)) - (phases2 - expectPhases(2));
    plot(xplot, yplot, 'ko-')
    title('Phase difference from expected')
    xlim(setaxlim(xplot, 1))
end

function plotstandingese(figNum, x, y, Param, bestPhase, intgs)
    nMeas = numel(bestPhase);
    f = figure(figNum);
    f.Name = "Phase";
    clf
    tiledlayout("flow", "TileSpacing", "none", "Padding", "tight")
    nexttile([1, 2])
    plot(1:nMeas, real(intgs), '.-')
    % yyaxis right
    hold on
    plot(1:nMeas, imag(intgs), '.-')
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
        yticks([])
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
