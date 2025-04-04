clearvars
addpath(genpath('util/'))

%% IMPORT

loadPath = '../data/processed/ZePSI-E-013012-Nut1.mat';
load(loadPath)
xAmp = table2array(readtable(...
    "../data/raw/ZePSI-E-013012/ZePSI-E-013012-param.txt", "Range", "B1"));

%% RABI NUTATIONS EXPERIMENT PULSE 1

disp('Integrate echo...')

iMax = 196;
integWidth = 90;

nTau = size(y{1}, 1);
nMeas = numel(y);
figure(1)
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
    % sax = ScrollableAxes('Index', 1);
    % plot(sax, x{1}, x{2}, real(y{ii}));
    hold on
    % plot(sax, x{1}, x{2}, imag(y{ii}));
    % plot(sax, x{1}, x{2}, real(integWindow{ii})*0.2e4);
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
    % plot(xrabi, real(rabii(ii, :)), 'o-')
    hold on
    % plot(xrabi, imag(rabii(ii, :)), 'o-')
    % plot(xrabi, winham*max(real(rabii(ii, :))))
    % plot(xrabi, winham.*real(rabii(ii, :)))
end

%% FIT cosine

rabii = rabii(:, 1:nTau);  % Recover rabii before zero filling
expcos = @(xx, A, tau, w) A*exp(-xx/tau).*cos(2*pi*w*xx);
fitmodel = @(p) expcos(xrabi', 1, p(1), p(2));  % A = 1

fitOpt = optimoptions('lsqnonlin','Display','off');

figure(100)
clf
w0 = (linspace(sqrt(3), sqrt(48), nMeas)).^2*1e-3;
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
    plotText = sprintf("%.2f", xAmp(ii));
    text(gca, 0.8, 0.8, plotText{end}, 'Units', 'normalized')
    xline(32)
end

%%
% ------------------------------------------------------------------------
FREQ_PI_PULSE = 46.74;  % MHz, taken from the rabiNut2.m data analysis file
% ------------------------------------------------------------------------

for ii = 1:nMeas
    freq(ii) = pfit{ii}(2);
end

diff = zeros(1, nMeas - 1);
for ii = 1:nMeas - 1
    diff(ii) = freq(ii + 1) - freq(ii);
end
figure()
plot(1:nMeas - 1, diff)

% xAmp = [78.00, 79.39, 79.79, 80.04, 80.23, 80.39, 80.53, 80.65, 80.76, ...
%     80.86, 80.96, 81.05, 81.15, 81.24, 81.34, 81.43, 81.53, 81.63, 81.74, ...
%     81.86, 82.00, 82.16, 82.35, 82.60, 83.00];



figure(74)
clf
errorbar(xAmp, freq, pci(:, 2, 1) - freq', ...
    pci(:, 2, 2) - freq', 'o')
ylim(setaxlim([1, FREQ_PI_PULSE]*1e-3, 1.05))
ylabel('Rabi freq / GHz')
% yyaxis right
% plot(xAmp, freq/FREQ_PI_PULSE*1e3, 'Visible', 'off')
xlim(setaxlim(xAmp, 1.05))
% ylim(setaxlim(freq/FREQ_PI_PULSE*1e3, 1.05))
% ylabel('Fractions of pi')
yline(FREQ_PI_PULSE*1e-3)

loadEseemPath = '../data/processed/ZePSI-E-013012-ESEEM.mat';  % MODIFY
LoadEseem = load(loadEseemPath);
turningAngle = zeros(nMeas);
for ii = 1:nMeas
    turningAngle(ii) = pfit{ii}(2)*1e3*pi/FREQ_PI_PULSE;
    LoadEseem.Param{ii}.turningAngle = turningAngle(ii);
end
save(loadEseemPath, '-struct', 'LoadEseem')

%% FIT SIGMOIDAL FUNCTION

sigfun = @(xx, p) p(1)./(1 + p(2)*exp(-p(3)*(xx - p(4)))) + p(5);
xx = linspace(78, 83, 1000);
p0 = [5, 1.8, 80.5];
fitmodel = @(p) sigfun(xAmp, [1, p, 0]);
yover = sigfun(xx, [0.04, p0, 0]);

% [pfit{ii}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
%     @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
ydata = freq';
[pfitsig, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
    @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
pcisig = nlparci(pfitsig, residual, 'jacobian', jacobian);
[yfitsig, pfitsigA, pfitsig(end + 1)] = mldividefun(...
    fitmodel, ydata, pfitsig);
pfitsig = [pfitsigA, pfitsig];

clf
errorbar(xAmp, freq, pci(:, 2, 1) - freq', ...
    pci(:, 2, 2) - freq', 'o')
hold on
% plot(xx, yover)
plot(xAmp, yfitsig)
plot(xx, sigfun(xx, pfitsig))

% save('../data/processed/ZePSI_013012_mpfuYch.mat', 'sigfun', 'pfitsig')
