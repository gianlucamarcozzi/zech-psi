clearvars
addpath(genpath('util/'))

%% IMPORT

loadPath = '../data/processed/ZePSI-E-013002-ESEEM.mat';
load(loadPath)

%% ESEEM INTEGRATION

% Find and plot integration window
integWidth = 150;
iMax = 213;

nMeas = numel(y);
nTau = 200;
figure(1)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
% for ii = nMeas - 3:nMeas
for ii = 15:nMeas
    [integWindow{ii}, integWindowPlot{ii}] = deal(zeros(size(y{ii})));
    for itau = 1:nTau
        valMax = max(real(y{ii}(itau, :)));
        iInteg = (iMax - integWidth/2 + itau*2 - 2):(iMax + integWidth/2 + ...
            itau*2 - 2);
        integWindow{ii}(itau, iInteg) = ones(1, integWidth + 1);
        integWindowPlot{ii}(itau, iInteg) = ...
            valMax*integWindow{ii}(itau, iInteg);
    end
    y2{ii} = integWindow{ii}.*y{ii};  % Signal in the window

    nexttile
    sax = ScrollableAxes('Index', 1);
    plot(sax, x{1}, x{2}, real(y{ii}));
    hold on
    plot(sax, x{1}, x{2}, imag(y{ii}));
    plot(sax, x{1}, x{2}, real(integWindowPlot{ii}));
    % plot(sax, x{1}, x{2}, imag(y2{ii}));
    ylim(setaxlim([-3e4, 6e4], 1.1))
    title(Param{ii}.turningAngle)
    % xline([-1, 1]/2*integWidth + x{1}(iMax), '--')
    % fprintf("%d max at:\t%d\n", [ii, iMax])
end

% Integrate to get Eseem
eseem = zeros(nMeas, nTau);
figure(22)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
xeseem = 150 + (0:2:(nTau - 1)*2);
for ii = 15:nMeas
    eseem(ii, :) = sum(y2{ii}(1:nTau, :), 2);

    nexttile
    plot(xeseem, real(eseem(ii, :)))
    hold on
    plot(xeseem, imag(eseem(ii, :)))

    title(Param{ii}.turningAngle/pi*100)
end

figure(23)
clf
for ii = [15, 24]
    plot(xeseem, imag(eseem(ii, :)))
    hold on
end

%% PLOTS

iy = 14;
nSpectrum = 5;
i1 = round(linspace(1, size(y{iy}, 1), nSpectrum));
i2 = 30:size(y{iy}, 2);
% cmap = viridis(nSpectrum);
% figure()
clf
tL = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
yplot = real(y{iy})*1e-4;
for ii = 1:nSpectrum
    plot(x{1}(i2), yplot(i1(ii), i2),  'DisplayName', append(...
        char(hex2dec('03C4')), sprintf(' = %d ns', 150 + i1(ii) - 1)))
    hold on
end
xlim(setaxlim(x{1}(i2), 1))
ylim(setaxlim(yplot(i1(1), i2), 1.05))
xticks(100:200:1000)
yticks(0:0.2:0.6)
legend()
nexttile
yplot = imag(y{iy})*1e-4;
for ii = 1:nSpectrum
    plot(x{1}(i2), yplot(i1(ii), i2))  % + (ii - 1)*0.6e4
    hold on
end
xlim(setaxlim(x{1}(i2), 1))
ylim(setaxlim([yplot(i1(1), i2), yplot(i1(2), i2)], 1.05))
xticks(100:200:1000)
yticks(0:2)
labelaxesfig(tL, "Time / ns",...
    ["In-phase / a.u.", "Out-of-phase / a.u."])
% nexttile(1); ylabel("In-phase channel / a.u.")
% nexttile(2); ylabel("Out-of-phase channel / a.u.")
savePath = "../images/" + "ZePSI-013002-04-eseemRingdownSubtracted.png";
% saveas(gcf, savePath)

%%

ii = 14;
figure(2)
clf
xeseem = 150 + (0:2:(nTau - 1)*2);
yplot = eseem(ii, :)*1e-5;
plot(xeseem, real(yplot))
hold on
plot(xeseem, imag(yplot))

xlim(setaxlim(xeseem, 1))
% ylim(setaxlim([yplot(i1(1), i2), yplot(i1(2), i2)], 1.05))
xticks(200:100:600)
yticks(-9:3:9)
labelaxesfig(gca, "Time / ns", "ESEEM / a.u.")
legend("In-phase", "Out-of-phase")
savePath = "../images/" + "ZePSI-013002-04-eseemPiOver20.png";
% saveas(gcf, savePath)

%%

expcos = @(xx, A, tau, w, phi) A*exp(-xx/tau).*cos(2*pi*w*xx + pi);
fitmodel = @(p) expcos(xeseem, 1, p(1), p(2));  % A = 1

fitOpt = optimoptions('lsqnonlin','Display','off');
figure(7)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
w0 = [3.5]*1e-3;
for ii = 14

    p0 = [500, w0];
    ydata = imag(eseem(ii, :));
    ydata = ydata'/max(ydata);
    clf
    plot(xeseem, ydata, 'o') 
    hold on
    plot(xeseem, fitmodel(p0))

    [pfit{ii}, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(...
        @(p) ydata - mldividefun(fitmodel, ydata, p)', ...
        p0, [], [], fitOpt);
    pci(ii, :, :) = nlparci(pfit{ii}, residual, 'jacobian', jacobian);

    [yfit{ii}, pfit{ii}(end + 1), pfit{ii}(end + 2)] = mldividefun(...
        fitmodel, ydata, pfit{ii});
    % nexttile
    % clf
    % plot(xeseem, ydata, 'o')
    hold on
    plot(xeseem, yfit{ii})
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
xlim(setaxlim(1:nMeas, 1.1))

%% FFT

tStep = xeseem(2) - xeseem(1);
fSampl = 1/tStep;
nzf = 1000;  % Zero filling
if nzf ~= 0
    fxeseem = fSampl/nzf*(-nzf/2:nzf/2 - 1);
else
    fxeseem = fSampl/nTau*(-nTau/2:nTau/2 - 1);
end

figure(4)
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 14
    if nzf ~= 0
        eseem(ii, nzf) = 0;  % Zero filling
        % eseem2(ii, nzf) = 0;  % Zero filling
    end
    % feseem = zeros()
    % feseem(ii, :) = fft(eseem(ii, :) - mean(eseem(ii, :)));
    feseem(ii, :) = fft(eseem(ii, :));
    % feseem2(ii, :) = fft(eseem2(ii, :) - mean(eseem(ii, :)));

    [~, ifMax(ii)] = max(abs(fftshift(feseem(ii, :))));
    % [~, ifMax2(ii)] = max(abs(fftshift(feseem2(ii, :))));

    nexttile
    plot(fxeseem, abs(fftshift(feseem(ii, :))), 'o-')
    % yyaxis right
    % plot(fxeseem, abs(fftshift(feseem2(ii, :))), 'o-')
    xlim([-0.1, 0.1])
end


%%
%{
% Fit Hahn-echo decay
expdecay = @(p) p(1)*exp(-(xeseem)/p(2));
iplot = 1;
for ii = 1:nMeas
    p0 = [1e5, 200];
    [pfit{ii}, ~, residual, ~, ~, ~, jacobian] = ...
        lsqnonlin(@(p) hedecay(ii, :) - expdecay(p), p0);
    pci(ii, :, :) = nlparci(pfit{ii}, residual, 'jacobian', jacobian);

    nexttile(iplot)
    hold on
    plot(xeseem, expdecay(pfit{ii}))
    iplot = iplot + 1;
end

% Plot best fit results
iarr = 1;
for ii = ies
    ampls(iarr) = pfit{ii}(1);
    taus(iarr) = pfit{ii}(2);
    iarr = iarr + 1;
end

figure(3)
clf
errorbar(1:numel(ies), ampls, ...
    pci(ies, 1, 1) - ampls', pci(ies, 1, 2) - ampls', 'o-')
yyaxis right
errorbar(1:numel(ies), taus, ...
    pci(ies, 2, 1) - taus', pci(ies, 2, 2) - taus', 'o-')
xlim(setaxlim(1:numel(ies), 1.1))

% Traces imag
figure(2)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:nMeas
    nexttileend
    sax = ScrollableAxes();
    plot(sax, x0{ii}{1}, x0{ii}{2}, imag(y0{ii}));
    hold on
end

%}
