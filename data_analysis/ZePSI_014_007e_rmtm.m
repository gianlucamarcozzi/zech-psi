%%
clearvars
addpath(genpath('util/'))

%% PULSE 2
generalFolder = "../data/raw/";
expName = "ZePSI-E-014-007";

% x-axis
pulseAmps = linspace(0.06, 0.003, 5)./0.06;
pulseAmps = pulseAmps(1:end-1);

expFolder = generalFolder + expName;
% BrX and BrMinX standing echo data (p plus, m minus)
[yrm1, ytm1, yrm2, ytm2, xsc] = ...
    importrmtm(expFolder, expName);
% Number of measurements
nMeas = size(yrm1, 1);

%% OSCILLOSCOPE
N_PULSE = 2;

% Define range of the x-axis where the pulses are expected
XRANGETM{1} = [3.05, 3.10]*1e-7;
XRANGETM{2} = [7.05, 7.20]*1e-7;
XRANGERM{1} = [3.16, 3.54]*1e-7;
XRANGERM{2} = [7.16, 7.72]*1e-7;
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
plotscopedata(1, "RM BrX", xsc, yrm1 - blrm1, xrm, yrm1cut, idxPlots)
plotscopedata(3, "TM BrX", xsc, ytm1 - bltm1, xtm, ytm1cut, idxPlots)
plotscopedata(5, "RM BrMinX", xsc, yrm2 - blrm2, xrm, yrm2cut, idxPlots)
plotscopedata(7, "TM BrMinX", xsc, ytm2 - bltm2, xtm, ytm2cut, idxPlots)
%%
plotscopedatacmp(9, "RM BrMinX", xsc, (yrm2 - blrm2), xsc, yrm1 - blrm1, idxPlots)
plotscopedatacmp(11, "TM BrMinX", xsc, (ytm2 - bltm2), xsc, ytm1 - bltm1, idxPlots)
%%
plotscopefinal(13, pulseAmps, ...
    meantm1{1}, meantm2{1}, meantm1{2}, ...
    amptm1{1}, amptm2{1}, amptm1{2}, ...
    phtm1{1}, phtm2{1}, phtm1{2}, [0, 180])
plotscopefinal(14, pulseAmps, ...
    meanrm1{1}, meanrm2{1}, meanrm1{2}, ...
    amprm1{1}, amprm2{1}, amprm1{2}, ...
    phrm1{1}, phrm2{1}, phrm1{2}, [0, 180])

%% ------------------------- STANDING ESE --------------------------------
% ADJUST PHASE
% Get variables
loadPath0 = "../data/processed/" + expName;
st{1} = load(loadPath0 + "-Standing-ESE-plus.mat");
st{2} = load(loadPath0 + "-Standing-ESE-minus.mat");
OptBl = struct('polyOrder', 1, 'range', [0, 50; 400, 600]);
for iset = 1:2
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
        for itau = 1:size(st{iset}.y{ii}, 1)
            st{iset}.y2{ii} = shiftphase(y1{ii}, bestPhase(ii));
        end
    end

    % plotstandingese(2*iset - 1, x, y1, Param, bestPhase, maxval, imagval)
    % plotstandingese(2*iset, x, st{iset}.y2, Param, bestPhase, maxval, imagval)
end

% plotsumdiffese(5, x, st{1}.y, st{2}.y)
% plotsumdiffese(6, x, st{1}.y1, st{2}.y2)

ypbl = st{1}.y2;
ymbl = st{2}.y2;
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
plotsumdiffese(15, st{1}.x, ypbl, ymbl);
plotstandingese(16, pulseAmps, ypamp, ypph, ymamp, ymph);

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

%% ------------------------- FUNCTIONS -----------------------------------
function [yrm1, ytm1, yrm2, ytm2, x] = importrmtm(expFolder, expName)
    pathpart1 = append(expFolder, "/", expName);
    % Pulse 1
    ar = importdata(pathpart1 + "-scope-p1rm1.txt");
    ai = importdata(pathpart1 + "-scope-p1rm2.txt");
    br = importdata(pathpart1 + "-scope-p1tm1.txt");
    bi = importdata(pathpart1 + "-scope-p1tm2.txt");
    % y values pulse 1
    yrm1 = ar(2:end, :) + 1i*ai(2:end, :);
    ytm1 = br(2:end, :) + 1i*bi(2:end, :);
    % Pulse 2
    ar = importdata(pathpart1 + "-scope-p2rm1.txt");
    ai = importdata(pathpart1 + "-scope-p2rm2.txt");
    br = importdata(pathpart1 + "-scope-p2tm1.txt");
    bi = importdata(pathpart1 + "-scope-p2tm2.txt");
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

function plotscopedata(figNum, figNameStart, x, y, xcut, ycut, idxPlots)
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
            yplot = y(iplot, :);
            yplot = feval(funplot, yplot);
            plot(xplot, yplot)
            hold on
            for ip = 1:numel(ycut)
                yplot = ycut{ip}(iplot, :);
                xplot = xcut{ip};
                yplot = feval(funplot, yplot);
                plot(xplot, yplot)
            end
        end
    end
end

function plotscopedatacmp(figNum, figNameStart, x, y, x2, y2, idxPlots)
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
            yplot = y(iplot, :);
            yplot = feval(funplot, yplot);
            plot(xplot, yplot)
            xplot = x2;
            yplot = y2(iplot, :);
            yplot = feval(funplot, yplot);
            plot(xplot, yplot)            
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

function plotstandingese(figNum, x, ypamp, ypph, ymamp, ymph)
    f = figure(figNum);
    f.Name = "Stanging-ESE: Amplitude and phase";
    tL = tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact");
    % Amplitude
    nexttile(1)
    hold on
    box on
    plot(x, ypamp, 'o-')
    plot(x, ymamp, 'o-')
    xlim(setaxlim(x, 1.05))
    title("Amplitude")
    
    % Relative amplitude diff
    nexttile(3)
    yplot = (ypamp - ymamp)./ypamp;
    plot(x, yplot, 'ko-')
    xlim(setaxlim(x, 1.05))
    title("Relative amplitude difference")
    
    % Phase
    nexttile(2)
    hold on
    box on
    plot(x, ypph, 'o-')
    plot(x, ymph + 180, 'o-')
    xlim(setaxlim(x, 1.05))
    legend("phase 1", "phase 2 + 180")
    title("Phase")
    
    % Phase diff from expected
    nexttile(4)
    yplot = (ypph - 0) - (ymph + 180);
    plot(x, yplot, 'ko-')
    xlim(setaxlim(x, 1.05))
    title("Phase difference from expected")
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