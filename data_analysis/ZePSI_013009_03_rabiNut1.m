clearvars
addpath(genpath('util/'))

%% IMPORT

loadPath = '../data/processed/ZePSI-E-013009-Nut1.mat';
load(loadPath)

%% RABI NUTATIONS EXPERIMENT PULSE 1

disp('Integrate echo...')

iMax = 196;
integWidth = 90;

nTau = size(y{1}, 1);
nMeas = numel(y);
figure(99)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    [integWindow{ii}, integWindowPlot{ii}] = deal(zeros(size(y{ii})));
    for itau = 1:nTau
        valMax = max(real(y{ii}(itau, :)));
        iInteg = (iMax - integWidth/2):(iMax + integWidth/2);
        integWindow{ii}(itau, iInteg) = ones(1, integWidth + 1);
    end
    y2{ii} = integWindow{ii}.*y{ii};  % Signal in the window

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x{1}, x{2}, real(y{ii}));
    hold on
    plot(sax, x{1}, x{2}, imag(y{ii}));
    plot(sax, x{1}, x{2}, real(integWindow{ii})*0.2e4);
    ylim([-0.5e4, 0.5e4])
end

% Integrate to get Rabi nutations
xrabi = x{2}*2 + 2;  % Each step is an increment of 2 ns of the pulse length
rabii = zeros(nMeas, nTau);
% Create window function
winham = windowhamming(nTau, 1);
for ii = 1:nMeas
    rabii(ii, :) = sum(y2{ii}(1:nTau, :), 2);
    % rabii2(ii, :) = sum(y22{ii}(1:nTau, :), 2);

    nexttile
    plot(xrabi, real(rabii(ii, :)), '-')
    hold on
    plot(xrabi, imag(rabii(ii, :)), '-')
    plot(xrabi, winham*max(real(rabii(ii, :))))
    plot(xrabi, winham.*real(rabii(ii, :)))
end

%% FIT cosine

rabii = rabii(:, 1:nTau);  % Recover rabii before zero filling
expcos = @(xx, A, tau, w) A*exp(-xx/tau).*cos(2*pi*w*xx);
fitmodel = @(p) expcos(xrabi', 1, p(1), p(2));  % A = 1

fitOpt = optimoptions('lsqnonlin','Display','off');

clf
w0 = (linspace(1, sqrt(40), nMeas)).^2*1e-3;
for ii = 1:nMeas
    p0 = [100, w0(ii)];
    ydata = real(rabii(ii, :))/max(real(rabii(ii, :)));

    [pfit{ii}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pci(ii, :, :) = nlparci(pfit{ii}, residual, 'jacobian', jacobian);

    [yfit{ii}, pfit{ii}(3), pfit{ii}(4)] = mldividefun(...
        fitmodel, ydata, pfit{ii});
    nexttile
    plot(xrabi, ydata, 'o')
    hold on
    plot(xrabi, yfit{ii})
    title(sprintf('%d: %.4f MHz, %.2f ns', ...
          ii, pfit{ii}(2)*1e3, 1/2/pfit{ii}(2)))
    % plotText = strsplit(Param{ii}.TITL, '-');
    % text(gca, 0.8, 0.8, plotText{end}, 'Units', 'normalized')
end

% ------------------------------------------------------------------------
FREQ_PI_PULSE = 43.4;  % MHz, taken from the rabiNut2.m data analysis file
% ------------------------------------------------------------------------
figure()
clf
for ii = 1:nMeas
    freq(ii) = pfit{ii}(2);
end
xAmp = [78:0.2:81, 82.4:0.2:83];
errorbar(xAmp, freq, pci(:, 2, 1) - freq', ...
    pci(:, 2, 2) - freq', 'o')
ylim(setaxlim(freq, 1.05));
ylabel('Rabi freq / MHz')
yyaxis right
plot(xAmp, freq/FREQ_PI_PULSE*1e3, 'Visible', 'off')
xlim(setaxlim(xAmp, 1.05))
ylim(setaxlim(freq/FREQ_PI_PULSE*1e3, 1.05))
ylabel('Fractions of pi')

loadEseemPath = '../data/processed/ZePSI-E-013009-ESEEM.mat';
LoadEseem = load(loadEseemPath);
turningAngle = zeros(nMeas);
for ii = 1:nMeas
    turningAngle(ii) = pfit{ii}(2)*1e3*pi/FREQ_PI_PULSE;
    LoadEseem.Param{ii}.turningAngle = turningAngle(ii);
end
% save(loadEseemPath, '-struct', 'LoadEseem')

%% FIT SIGMOIDAL FUNCTION

sigfun = @(xx, p) p(1)./(1 + p(2)*exp(-p(3)*(xx - p(4)))) + p(5);
xx = linspace(78, 83, 1000);
p0 = [5, 1.8, 80.5];
fitmodel = @(p) sigfun(xAmp, [1, p, 0]);
yover = sigfun(xx, [0.04, p0, 0]);

% [pfit{ii}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
%     @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
ydata = freq;
[pfitsig, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
    @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
pcisig = nlparci(pfitsig, residual, 'jacobian', jacobian);
[yfitsig, pfitsigA, pfitsig(end + 1)] = mldividefun(...
    fitmodel, ydata, pfitsig);
pfitsig = [pfitsigA, pfitsig];

clf
errorbar(xAmp, freq', pci(:, 2, 1) - freq', ...
    pci(:, 2, 2) - freq', 'o')
hold on
% plot(xx, yover)
plot(xAmp, yfitsig)
plot(xx, sigfun(xx, pfitsig))

% save('../data/processed/ZePSI_013009_mpfuYch.mat', 'sigfun', 'pfitsig')

%% FFT

APPLY_WINDOW = 0;
tStep = xrabi(2) - xrabi(1);
fSampl = 1/tStep;
nzf = 1024;  % Zero filling
if APPLY_WINDOW
    winrabii = winham.*rabii;  % Apply window function
else
    winrabii = rabii;
end
if nzf ~= 0
    fxrabi = fSampl/nzf*(-nzf/2:nzf/2 - 1);  % Freq axis
    frabii = zeros(nMeas, nzf);  % Initialize fft arrays
    winrabii(:, nzf) = zeros(nMeas, 1);  % Zero filling
else
    fxrabi = fSampl/nTau*(-nTau/2:nTau/2 - 1);  % Freq axis
    frabii = zeros(nMeas, nTau);  % Initialize fft arrays
end

figure()
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    frabii(ii, :) = fft(winrabii(ii, :));

    nexttile
    plot(fxrabi, abs(fftshift(frabii(ii, :))), 'o-')
    xlim([-0.1, 0.1])
end
