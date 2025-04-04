clearvars
addpath(genpath('util/'))

%% IMPORT

load('../data/processed/ZePSI-E-013031-Nut2.mat')
nMeas = numel(x);

%% FIT cosine

expcos = @(xx, A, tau, w) A*exp(-xx/tau).*cos(2*pi*w*xx);
fitmodel0 = @(xx, p) expcos(xx, 1, p(1), p(2));  % A = 1

fitOpt = optimoptions('lsqnonlin','Display','off');

clf
w0 = 47*1e-3;
for ii = 1:nMeas
    p0 = [40, w0];
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
errorbar(1:nMeas, freq, pci(:, 2, 1) - freq', pci(:, 2, 2) - freq', 'o')
xlim([0.9, nMeas + 0.1])
deltapci = pci(:, 2, 2) - pci(:, 2, 1);
fprintf("Average Rabi frequency: %.2f\n", sum(freq*1e3./deltapci')/sum(1./deltapci))
