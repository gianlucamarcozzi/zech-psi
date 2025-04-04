clearvars
addpath(genpath('util/'))

%% IMPORT

loadPath = '../data/processed/ZePSI-E-013023-Nut1.mat';
load(loadPath)
nMeas = numel(x);

%% FIT cosine

expcos = @(xx, A, tau, w) A*exp(-xx/tau).*cos(2*pi*w*xx);
fitmodel0 = @(xx, p) expcos(xx, 1, p(1), p(2));  % A = 1

fitOpt = optimoptions('lsqnonlin','Display','off');

clf
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
w0 = (linspace(sqrt(6), sqrt(37), nMeas)).^2*1e-3;
for ii = 1:nMeas
    p0 = [40, w0(ii)];
    fitmodel = @(p) fitmodel0(x{ii}, p);
    ydata = real(y{ii});
    ydata = ydata/max(ydata);

    [pfit{ii}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p), p0, [], [], fitOpt);
    pci(ii, :, :) = nlparci(pfit{ii}, residual, 'jacobian', jacobian);

    [yfit{ii}, pfit{ii}(3), pfit{ii}(4)] = mldividefun(...
        fitmodel, ydata, pfit{ii});

    nexttile
    plot(x{ii}, ydata, 'o')
    hold on
    plot(x{ii}, yfit{ii})
    title(sprintf('%d: %.4f MHz, %.2f ns', ...
          ii, pfit{ii}(2)*1e3, 1/2/pfit{ii}(2)))
end

nexttile
for ii = 1:nMeas
    freq(ii) = pfit{ii}(2);
end
% errorbar(1:nMeas, freq, pci(:, 2, 1) - freq', pci(:, 2, 2) - freq', 'o')
errorbar(1:nMeas, freq, pci(:, 2, 1) - freq', ...
    pci(:, 2, 2) - freq', 'o')
xlim([0.9, nMeas + 0.1])

%%
% ------------------------------------------------------------------------
FREQ_PI_PULSE = 39.97;  % MHz, taken from the rabiNut2.m data analysis file
% ------------------------------------------------------------------------

for ii = 1:nMeas
    freq(ii) = pfit{ii}(2);
end

diff = zeros(1, nMeas - 1);
for ii = 1:nMeas - 1
    diff(ii) = freq(ii + 1) - freq(ii);
end

for ii = 1:nMeas
    xAmp(ii) = Param{ii}.mpfuXAmp;
end

figure(74)
clf
errorbar(xAmp, freq, pci(:, 2, 1) - freq', ...
    pci(:, 2, 2) - freq', 'o')
ylim(setaxlim([1, FREQ_PI_PULSE]*1e-3, 1.05))
ylim(setaxlim([1, FREQ_PI_PULSE + 7]*1e-3, 1.05))
ylabel('Rabi freq / GHz')
% yyaxis right
% plot(xAmp, freq/FREQ_PI_PULSE*1e3, 'Visible', 'off')
xlim(setaxlim(xAmp, 1.05))
% ylim(setaxlim(freq/FREQ_PI_PULSE*1e3, 1.05))
% ylabel('Fractions of pi')
yline(FREQ_PI_PULSE*1e-3)

loadEseemPath = append(loadPath(1:end-9), '-ESEEM.mat');  % MODIFY
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
fitmodel = @(p) sigfun(xAmp', [1, p, 0]);
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

% save('../data/processed/ZePSI_013022_mpfuXch.mat', 'sigfun', 'pfitsig')
