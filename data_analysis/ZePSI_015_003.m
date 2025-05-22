% clearvars
addpath(genpath('util/'))

%% GENERAL
generalFolder = "../data/raw/";
pathToProcessed = "../data/processed/"; 
expNumbers = ["3", "7", "5", "9"];  % 150, 200, 230, 250 Kelvin

% -------------------------------------------------------------------------
% ------------------- SELECT NUMBER OF EXPERIMENT -------------------------
% -------------------------------------------------------------------------
INUM = 2;

%% ------------------- STORE TO MAT FILE ----------------------------------
expName = "ZePSI-E-015-00" + expNumbers(INUM);
fig0 = 1000*INUM;
measFolder = generalFolder + expName;
filepathParam = append(generalFolder, expName, '/', expName, '-param.txt');

% Store mpfu amplitude of first pulse in param structure
xAmp = getamplitudefirstpulse(filepathParam);

titles = ["-ESEEM-light", "-ESEEM-dark", "-Standing-ESE", ...
 "-Nut1", "-Nut2"];

% THIS SHOULD BE RUN ONLY ONCE
if true
    for ii = 1:numel(titles)
        ending = sprintf("*%i.DTA", ii);
        [x0, y, Param] = loadfolderelexsys(measFolder, ending);
        % Save pulse amplitudes
        for iamp = 1:numel(Param)
            Param{iamp}.pAmpPulse1 = xAmp(iamp);
        end
        % Save only one x-axis for ESEEM and Standing
        if contains(titles(ii), "ESEEM") || contains(titles(ii), "Standing")
            x = x0{1};
        else
            x = x0;
        end
        % Save
        savePath = pathToProcessed + expName + titles(ii) + ".mat";
        save(savePath, "x", "y", "Param");
    end
end

%% ------------------------- RABI PULSE 2 ---------------------------------
loadPath = "../data/processed/" + expName + "-Nut2.mat";
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

[freqRabi, piPerf, dfreqRabi] = deal(zeros(1, nMeas));
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
f = figure(fig0 + 1);
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

%% ---------------------------- PULSE 1 ----------------------------------
loadPath = "../data/processed/" + expName + "-Nut1.mat";
pul1 = load(loadPath);

f = figure(fig0 + 2);
f.Name = "Rabi nutations first pulse";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
% Get variables
x = pul1.x;
y = pul1.y;
Param = pul1.Param;
nMeas = numel(y);

% Variables relative to fitting
expcos = @(xx, p) [exp(-xx/p(1)).*cos(2*pi*p(2)*xx), ones(numel(xx), 1)];
fitOpt = optimoptions('lsqnonlin','Display','off');
% w0s = linspace(0.5, sqrt(30), nMeas).^2*1e-3;
w0s = linspace(0.5, 24, nMeas)*1e-3;
% w0s = flip(w0s);
p01 = 170;

%------------------------------ FIT -----------------------------------
[freqRabi, piPerf, dfreqRabi] = deal(zeros(1, nMeas));
for ii = 1:nMeas
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

% -------------------------- PLOT -------------------------------------
for ii = 1:nMeas
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
for ii = 1:nMeas
    pAmp(ii) = Param{ii}.pAmpPulse1;
end
nexttile(ii + 1)
hold on; box on;
START_IDX = 1;
errorbar(pAmp(START_IDX:end), freqRabi(START_IDX:end), ...
    dfreqRabi(START_IDX:end), 'o')
xlim(setaxlim(pAmp(START_IDX:end), 1.05))

%% RMTM
generalFolder = "../data/raw/";
% expName = "ZePSI-E-015-005";
expFolder = generalFolder + expName;

% Pulse amplitude of first pulse
filepathParam = append(expFolder, '/', expName, '-param.txt');
pAmps = getamplitudefirstpulse(filepathParam);

% BrX and BrMinX standing echo data (p plus, m minus)
[yrm1, ytm1, yrm2, ytm2, xsc] = importrmtm(expFolder, expName);
% Number of measurements
nMeas = size(yrm1, 1);

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
idxPlots = [1, 2, 10, 11];
plotscopedata(fig0 + 3, xsc, yrm1, xrm, yrm1cut, blrm1, idxPlots)
plotscopedata(fig0 + 4, xsc, ytm1, xtm, ytm1cut, bltm1, idxPlots)
% plotscopedata(fig0 + 5, xsc, yrm2, xrm, yrm2cut, blrm2, idxPlots)
% plotscopedata(fig0 + 6, xsc, ytm2, xtm, ytm2cut, bltm2, idxPlots)

% plotscopedatacmp(fig0 + 11, "RM plus", xsc, (yrm2 - blrm2), xsc, yrm1 - blrm1, idxPlots)
% plotscopedatacmp(fig0 + 13, "TM plus", xsc, (ytm2 - bltm2), xsc, ytm1 - bltm1, idxPlots)
% plotscopedatacmp(fig0 + 15, "RM minus", xsc, -(yrm2 - blrm2), xsc, yrm1 - blrm1, idxPlots)
% plotscopedatacmp(fig0 + 17, "TM minus", xsc, -(ytm2 - bltm2), xsc, ytm1 - bltm1, idxPlots)

plotscopefinal(fig0 + 7, pAmps, ...
    meantm1{1}, meantm2{1}, meantm1{2}, ...
    amptm1{1}, amptm2{1}, amptm1{2}, ...
    phtm1{1}, phtm2{1}, phtm1{2}, [0, 180])
plotscopefinal(fig0 + 8, pAmps, ...
    meanrm1{1}, meanrm2{1}, meanrm1{2}, ...
    amprm1{1}, amprm2{1}, amprm1{2}, ...
    phrm1{1}, phrm2{1}, phrm1{2}, [0, 90])

%{ 
% SAVE IMAGE
savepathpdf = '../images/014-009_cmpRMTM.pdf';
delete(savepathpdf);
for ii = [15, 16, 13, 14]
    % exportgraphics(figure(ii), savepathpdf, "Append", true)
end
%}

%% ------------------------- STANDING ESE --------------------------------
loadPath0 = "../data/processed/" + expName;
st = load(loadPath0 + "-Standing-ESE.mat");
OptBl = struct('polyOrder', 1, 'width', 0.15);

x = st.x;
y = st.y;
Param = st.Param;
nMeas = numel(y);

% Baseline correction
for ii = 1:nMeas
    [y1{ii}, bl{ii}] = correctbaseline(x, y{ii}, OptBl);
end
% Find best phase
shiftphase = @(y, p) y*exp(1i*p*pi/180);  % degrees
idxs = 150:240;
for ii = 1:nMeas
    [yout{ii}, bestPhase(ii)] = correctphase(y1{ii}(idxs), 'pulse');
end

% SOME CORRECTIONS!!
% bestPhase(1) = bestPhase(2);

for ii = 1:nMeas
    for itau = 1:size(st.y{ii}, 1)
        st.y2{ii} = shiftphase(y1{ii}, bestPhase(ii));
        intgs(ii) = trapz(st.y2{ii});
    end
end

plotstandingese(fig0 + 11, x, y1, Param, bestPhase, intgs)
plotstandingese(fig0 + 12, x, st.y2, Param, bestPhase, intgs)

figure(fig0 + 13)
clf
tL = tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");
nexttile; hold on; box on;
IDX = 8;  % For normalization
plot(pAmps, amprm1{1}/amprm1{1}(IDX), 'o-')
plot(pAmps, (amptm1{1}/amptm1{1}(IDX)), 'o-')
plot(pAmps, freqRabi/freqRabi(IDX), 'o-')
legend('RM', 'TM', 'Rabi freq', 'Location', 'northwest')
title("Magnitude")
nexttile; hold on; box on;
IDX_FOR_MEAN = 10:25;
deltaphase = mean(bestPhase(IDX_FOR_MEAN) + phrm1{1}(IDX_FOR_MEAN)');
plot(pAmps, -phrm1{1} + deltaphase, 'o-')
plot(pAmps, bestPhase, 'ko-')
% plot(1:nMeas, pPhases + deltaphase, 'ko-')
legend('RM', 'Standing ESE', 'Location', 'northwest')
title("Phase")
labelaxesfig(tL, "TM amp pulse 1 / V", "")
%------------------------ SAVE PHASES ------------------------------------
% pPhases = [phrm1{1}, phrm1{2}];
pPhases = -phrm1{1} + deltaphase;
savePath = pathToProcessed + expName + "-pPhases.mat";
save(savePath, "pPhases", "-mat")

%%
clf
plotscopedata(1, xsc, yrm1, xrm, yrm1cut, blrm1, [1, 2, 24, 25])
% plotscopedata(1, xsc, yrm1*3, xrm, yrm1cut, blrm1, 1:4)r

%% ESEEM PREPARATION
loadPath0 = "../data/processed/" + expName;
el = load(loadPath0 + "-ESEEM-light.mat");
ed = load(loadPath0 + "-ESEEM-dark.mat");
load(loadPath0 + "-pPhases.mat");
nMeas = numel(el.y);

% ADJUST PHASE
% deltaphase = bestPhase(end) - pPhases(end);
shiftphase = @(y, p) y*exp(1i*p*pi/180);  % degrees
for ii = 1:nMeas
    el.y1{ii} = shiftphase(el.y{ii}, pPhases(ii));
    ed.y1{ii} = shiftphase(ed.y{ii}, pPhases(ii));
    % el.y1{ii} = shiftphase(el.y{ii}, bestPhase(ii));
    % ed.y1{ii} = shiftphase(ed.y{ii}, bestPhase(ii));
end

% LIGHT MINUS DARK
% Get variables
x = struct('x1', el.x.t, ...
    'x2', el.x.b0 + getparampulsespel(el.Param{1}, 'd1 '));
Param = el.Param;

% Light minus dark (already phase corrected)
for ii = 1:nMeas
    yeseem{ii} = el.y1{ii} - ed.y1{ii};
end

fplot = "imag";
f = figure(fig0 + 14);
f.Name = "Phase shifted data";
clf
tiledlayout('flow', "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    hold on; box on;
    sax = ScrollableAxes('Index', 1);
    % plot(sax, x.x1, x.x2, feval(fplot, el.y{ii}));
    plot(sax, x.x1, x.x2, feval("real", el.y1{ii}));
    plot(sax, x.x1, x.x2, feval(fplot, el.y1{ii}));
    % plot(sax, x.x1, x.x2, feval(fplot, yeseem{ii}));
    ylim([-10, 14]*1e3)
    plotTitle = strsplit(Param{ii}.TITL, '-');
    title(plotTitle{end-1})
end

% ----------------------------- SAVE DATA --------------------------------
pathToProcessed = "../data/processed/";

% ESEEM
clear('x', 'y', 'Param')
savePath = pathToProcessed + expName + "-ESEEM.mat";
x = struct('x1', el.x.t, 'x2', el.x.b0);
y = yeseem;
Param = el.Param;
save(savePath, 'x', 'y', 'Param')

%% ------------------------ ESEEM ----------------------------------------
loadPath = "../data/processed/" + expName + "-ESEEM.mat";
load(loadPath)
nMeas = numel(y);

% Find and plot integration window
integWidth = 20;
iMax = 196;
deltaTau = getparampulsespel(Param{1}, 'd30');  % ns
nMeas = numel(y);

for ii = 1:nMeas
    [integWindow{ii}, integWindowPlot{ii}] = deal(zeros(size(y{ii})));
    nTau = size(y{ii}, 1);
    for itau = 1:nTau
        valMax = max(real(y{ii}(itau, :)));
        iInteg1 = iMax - integWidth/2 + (itau - 1)*deltaTau;
        iInteg2 = iMax + integWidth/2 + (itau - 1)*deltaTau;
        iInteg = iInteg1:iInteg2;
        integWindow{ii}(itau, iInteg) = ones(1, integWidth + 1);
        integWindowPlot{ii}(itau, iInteg) = ...
            valMax*integWindow{ii}(itau, iInteg);
    end
    y2{ii} = integWindow{ii}.*y{ii};  % Signal in the window
end

figure(fig0 + 15)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:6
    nexttile
    sax = ScrollableAxes();
    plot(sax, x.x1, x.x2, real(y{ii}));
    hold on
    plot(sax, x.x1, x.x2, imag(y{ii}));
    plot(sax, x.x1, x.x2, imag(y2{ii}));
    xline(iMax, '--')
    ylim(setaxlim([-1e4, 2e4], 1.1))
    plotTitle = strsplit(Param{ii}.TITL, '-');
    title(plotTitle{end - 1}(1:3))
end

% INTEGRATE TO GET ESEEM
eseem = zeros(nMeas, nTau);
d1 = getparampulsespel(Param{1}, 'd1 ');
xeseem = d1 + (0:deltaTau:(nTau - 1)*deltaTau);
for ii = 1:nMeas
    eseem(ii, :) = sum(y2{ii}(1:nTau, :), 2);
end

%%
% I_BEST = [1:5];
I_BEST = [18:22];
% I_BEST = [40:43];
figure(fig0 + 16)
clf

fplot = "imag";
cmap = viridis(nMeas);
for ii = 1:nMeas
    yplot = feval(fplot, eseem(ii, :));
    % displayName = sprintf("%d deg", round(Param{ii}.pAmpPulse1/0.046*180));
    displayName = sprintf("%d deg", round(amprm1{1}(ii)/amprm1{2}(ii)*180));
    plot(xeseem, yplot, '-',  'Color', cmap(ii, :), ...
        'DisplayName', displayName)
    hold on
end
xline(xeseem(I_BEST), "HandleVisibility", "off")
xlim(setaxlim(xeseem, 1))
ylim(setaxlim( ...
    [min(feval(fplot, eseem(:))), max(feval(fplot, eseem(:)))], 1.05))
yline(0, "HandleVisibility", "off")
% plot(xeseem, imag(eseemcorr), 'Color', 'r', 'DisplayName', 'Corr')
labelaxesfig(gca, "Time / ns", "Intensity / a.u.")
legend('NumColumns', 3, "FontSize", 9)

saveName = strsplit(loadPath, '-');
% saveas(gcf, append('../images/', saveName{3}, '-04-01-oop.png'))

ybeta = zeros(1, nMeas);
xbeta = zeros(1, nMeas);
for ii = 1:nMeas
    for jj = 1:numel(I_BEST)
        ybeta(ii, jj) = imag(eseem(ii, I_BEST(jj)));
        % ybeta(ii, jj) = yfin(ii, I_BEST(jj))
    end
    xbeta(ii) = Param{ii}.pAmpPulse1/0.049*pi;
    xbeta(ii) = amprm1{1}(ii)/amprm1{2}(ii)*pi;
end

figure(fig0 + 17)
clf
for ii = 1:numel(I_BEST)
    yplot = ybeta(:, ii)/max(abs(ybeta(:, ii)));
    yplot = -sign(yplot(ii))*yplot;
    xplot = xbeta*180/pi;
    displayName = sprintf("%d ns", xeseem(I_BEST(ii)));
    plot(xplot, yplot, 'o-', 'DisplayName', displayName)
    hold on
end
yline(0, 'HandleVisibility', 'off')
% xline(180, '--', 'HandleVisibility', 'off')
legend('Location', 'northwest')

aa = load("/home/gianluca/files/projects/oop-ciss-calculations/data/digitized/zech_p46_oopEseem_expData.csv");
bb = load("/home/gianluca/files/projects/oop-ciss-calculations/data/digitized/zech_p46_oopEseem_fit.csv");
plot(aa(:, 1), aa(:, 2)/max(abs(aa(:, 2))), 'kx-', 'DisplayName', 'Zech')
plot(bb(:, 1), bb(:, 2)/max(abs(aa(:, 2))), 'r.-', 'DisplayName', 'Zech fit')

xlim([0, 200])
ylim(setaxlim(aa(:, 2)/max(abs(aa(:, 2))), 1.05))
labelaxesfig(gca, "Beta / deg", "Intensity / a.u.")
saveName = strsplit(loadPath, '-');
% saveas(gcf, append('../images/', saveName{3}, '-04-02-beta.png'))

%% ESEEM TRACES AND FFT
% Create window function
nTau = size(eseem, 2);
winham = windowhamming(nTau, 40);
figure(fig0 + 18)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttile
    plot(xeseem, real(eseem(ii, :)), 'o-')
    hold on
    plot(xeseem, imag(eseem(ii, :)), 'o-')
    plot(xeseem, winham*max(real(eseem(ii, :))))
    plot(xeseem, winham.*imag(eseem(ii, :)))
    yline(0)
    title(Param{ii}.TITL(end - 4:end - 2))
end

% FFT
tStep = xeseem(2) - xeseem(1);
fSampl = 1/tStep;
nzf = 1024;  % Zero filling
wineseem = winham.*eseem;  % Apply window function
if nzf ~= 0
    fxeseem = fSampl/nzf*(-nzf/2:nzf/2 - 1);  % Freq axis
    wineseem(:, nzf) = zeros(nMeas, 1);  % Zero filling
    feseem = zeros(nMeas, nzf);  % Initialize fft arrays
else
    fxeseem = fSampl/nTau*(-nTau/2:nTau/2 - 1);  % Freq axis
    feseem = zeros(nMeas, nTau);  % Initialize fft arrays
end

figure(fig0 + 19)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    feseem(ii, :) = fft(imag(wineseem(ii, :)));

    nexttile(ii)
    plot(fxeseem, real(fftshift(feseem(ii, :))), '.-')
    hold on
    plot(fxeseem, imag(fftshift(feseem(ii, :))), '.-')
    yline(0)
    xlim(setaxlim(fxeseem))
    xlim([-0.015, 0.015])
end

%% FUNCTIONS
function xAmp = getamplitudefirstpulse(filepathParam)
    % Load the amplitude of the mpfu of the first pulse
    xAmp = table2array(readtable(filepathParam, "Range", "B1"));
    % Check for nans and cut the array as soon as one is found
    isnanxAmp = isnan(xAmp);
    if sum(isnanxAmp)
        for ii = 1:numel(xAmp)
            if isnanxAmp(ii) == 1
                xAmp = xAmp(1:ii - 1);
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

function plotscopedata(figNum, x, y, xcut, ycut, bl, idxPlots)
    funplots = ["real", "imag"];
    % figNames = figNameStart + funplots;
    f = figure(figNum);  % 2 to 4
    % f.Name = figNames(ifig);
    % funplot = funplots(ifig);
    clf
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
    for iplot = idxPlots
        nexttile(iplot)
        hold on; box on;
        xplot = x;
        for ii = 1:numel(funplots)
            yplot = y(iplot, :) - bl(iplot, :);
            yplot = feval(funplots(ii), yplot);
            plot(xplot, yplot)
        end
        xlim(setaxlim(xplot, 1))
        hold on
        for ii = 1:numel(funplots)
            for ip = 1:numel(ycut)
                yplot = ycut{ip}(iplot, :) - bl(iplot, :);
                xplot = xcut{ip};
                yplot = feval(funplots(ii), yplot);
                plot(xplot, yplot, 'r')
            end
        end
    end
    legend("Real ch", "Imag ch")
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
    axTitles = ["real", "imag"];
    funplots = axTitles;
    clf
    tL = tiledlayout(2, 3, "TileSpacing", "compact", "Padding", "compact");
    % Plots 1 and 3: real and imag part
    numTile = [1, 4];
    for itile = 1:2
        nexttile(numTile(itile))
        hold on
        box on
    
        xplot = pulseAmps;
        plot(xplot, feval(funplots(itile), y1), 'o-')  % <+x>
        % plot(xplot, feval(funplots(itile), y2), 'o-')  % <-x>
        plot(xplot, feval(funplots(itile), yX), 'ro-')  % +x pi
    
        title(sprintf("Scope ch %i", itile))
        xlim(setaxlim(xplot, 1))
    end
    
    % Amplitude
    nexttile(2)
    hold on
    box on
    plot(xplot, amps1, 'o-')  % <+x>
    % plot(xplot, amps2, 'o-')  % <-x>
    plot(xplot, ampsX, 'ro-')  % +x pi
    title('Magnitude')
    xlim(setaxlim(xplot, 1))
    
    % Phase
    nexttile(3)
    hold on
    box on
    plot(xplot, phases1, 'o-')  % <+x>
    % plot(xplot, phases2, 'o-')  % <-x>
    plot(xplot, phasesX, 'ro-')  % +x pi
    title('phase')
    xlim(setaxlim(xplot, 1))
    
    % Relative amplitude difference
    nexttile(5)
    hold on
    box on
    yplot = (amps1 - amps2)./amps1;
    % plot(xplot, yplot, 'ko-')
    title('Relative amp diff')
    xlim(setaxlim(xplot, 1))
    
    % Phase difference from expected
    nexttile(6)
    hold on
    box on
    yplot = (phases1 - phasesX - expectPhases(1)); % - (phases2 - expectPhases(2));
    plot(xplot, yplot, 'ko-')
    title('Phase difference from X pulse')
    xlim(setaxlim(xplot, 1))
    labelaxesfig(tL, "TM amp pulse 1 / V", "")
    ylabel("Diff / deg")
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
%{
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
%}
