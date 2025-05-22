%%
clearvars
addpath(genpath('util/'))

%% PULSE 2
pul2 = load('../data/processed/ZePSI-E-014-007-Nut2.mat');
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
f = figure(17);
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
ending = ["plus", "minus"];

f = figure(18);
f.Name = "Rabi nutations first pulse";
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for iphase = 1:2
    loadPath = "../data/processed/ZePSI-E-014-007-Nut1-" + ...
        ending(iphase) + ".mat";
    pul1 = load(loadPath);
    nMeas = numel(pul1.x);
    
    % Get variables
    x = pul1.x;
    y = pul1.y;
    Param = pul1.Param;
    
    % Variables relative to fitting
    expcos = @(xx, p) [exp(-xx/p(1)).*cos(2*pi*p(2)*xx), ones(numel(xx), 1)];
    fitOpt = optimoptions('lsqnonlin','Display','off');
    w0s = linspace(1, sqrt(30), nMeas).^2*1e-3;
    w0s = flip(w0s);
    p01 = 100;
    
    %------------------------------ FIT --------------------------------------
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
    
    % -------------------------- PLOT ----------------------------------------
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
    errorbar(pAmp, freqRabi, dfreqRabi, 'o')
    xlim(setaxlim(pAmp, 1.05))
end

% Taken from misc_E_002_003_mpfuBrX.m for the ones where the fit is bad
freqExpected = [2.7219, 3.5798, 4.4428, 5.2997, 6.1490, 7.0064, 7.8794];
freqRabiFin = [freqExpected, freqRabi(8:end)];
piPerfFin = [1/2./freqExpected*1e3, piPerf(8:end)];

%% ------------------------ SAVE INTO ESEEM STRUCTURES ---------------------
loadEseemPath = append(loadPath(1:end-9), '-ESEEM-light.mat');
loadEseem = load(loadEseemPath);
turningAngle = zeros(nMeas, 1);
pulseLength = getparampulsespel(loadEseem.Param{1}, 'p1 ');
for ii = 1:nMeas
    turningAngle(ii) = pulseLength/piPerfFin(ii)*pi*23.11/24; % *(pulseLength/RABI_LENGTH_PI_PULSE);
    loadEseem.Param{ii}.turningAngle = turningAngle(ii);
end
save(loadEseemPath, '-struct', 'loadEseem')

