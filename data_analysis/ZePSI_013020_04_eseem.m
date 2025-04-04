clearvars
addpath(genpath('util/'))

%% IMPORT

loadPath = '../data/processed/ZePSI-E-013020-ESEEM.mat';
load(loadPath)
nMeas = numel(y);
for ii = 1:nMeas
    xAmp(ii) = Param{ii}.mpfuXAmp;
end

%% ESEEM INTEGRATION

% Find and plot integration window
integWidth = 100;
iMax = 190;
DELTA_TAU = 8;  % ns
nMeas = numel(y);
figure(40)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
% for ii = nMeas - 3:nMeas
for ii = 1:nMeas
    [integWindow{ii}, integWindowPlot{ii}] = deal(zeros(size(y{ii})));
    nTau = size(y{ii}, 1);
    for itau = 1:nTau
        valMax = max(real(y{ii}(itau, :)));
        iInteg1 = iMax - integWidth/2 + (itau - 1)*DELTA_TAU;
        iInteg2 = iMax + integWidth/2 + (itau - 1)*DELTA_TAU;
        iInteg = iInteg1:iInteg2;
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
    ylim(setaxlim([-1e4, 2e4], 1.1))
    % titleStr = sprintf("%.1f deg, %.2f MPFU", ...
    %     Param{ii}.turningAngle/pi*180, xAmp(ii));
    % title(titleStr)
    % xline([-1, 1]/2*integWidth + x{1}(iMax), '--')
    % fprintf("%d max at:\t%d\n", [ii, iMax])
end

%% INTEGRATE TO GET ESEEM

eseem = zeros(nMeas, nTau);
xeseem = 180 + (0:8:(nTau - 1)*8);
for ii = 1:nMeas
    eseem(ii, :) = sum(y2{ii}(1:nTau, :), 2);
end

%%

I_BEST = [18];
figure(24)
clf
cmap = viridis(nMeas-3);
for ii = 1:nMeas-3
    displayName = sprintf("%d deg", round(Param{ii}.turningAngle/pi*180/1.07));
    if false
        plot(xeseem, imag(eseem(ii, :)), 'x', 'Color', cmap(ii, :), ...
            'DisplayName', displayName)
    else
        plot(xeseem, imag(eseem(ii, :)), '-', 'Color', cmap(ii, :), ...
            'DisplayName', displayName)
    end
    hold on
end
xlim(setaxlim(xeseem, 1))
ylim(setaxlim([min(imag(eseem(:))), max(imag(eseem(:))), 7e5], 1.05))
% xline(xeseem(I_BEST), "HandleVisibility", "off")
yline(0, "HandleVisibility", "off")
% plot(xeseem, imag(eseemcorr), 'Color', 'r', 'DisplayName', 'Corr')
labelaxesfig(gca, "Time / ns", "Intensity / a.u.")
legend('NumColumns', 3, "FontSize", 9, 'Location', 'north')


saveName = strsplit(loadPath, '-');
saveas(gcf, append('../images/', saveName{3}, '-04-01-oop.png'))

%%
ybeta = zeros(1, nMeas);
xbeta = zeros(1, nMeas);
for ii = 1:nMeas
    for jj = 1:numel(I_BEST)
        ybeta(ii, jj) = imag(eseem(ii, I_BEST(jj)));
    end
    xbeta(ii) = Param{ii}.turningAngle;    
end

figure(33)
clf
idxs = 1:size(ybeta, 1)-3;
for ii = 1:numel(I_BEST)
    yplot = ybeta(idxs, ii)/max(abs(ybeta(idxs, ii)));
    yplot = -sign(yplot(1))*yplot;
    xplot = xbeta(idxs)*180/pi./linspace(1, 1.07, numel(idxs));
    displayName = sprintf("%d ns", xeseem(I_BEST(ii)));
    plot(xplot, yplot, 'o-', 'DisplayName', displayName)
    hold on
end
yline(0, 'HandleVisibility', 'off')
legend('Location', 'northwest')

aa = load("/home/gianluca/files/projects/oop-ciss-calculations/data/digitized/zech_p46_oopEseem_expData.csv");
plot(aa(:, 1), aa(:, 2)/max(abs(aa(:, 2))), 'kx-', 'DisplayName', 'Zech')

xlim([0, 180])
ylim(setaxlim(aa(:, 2)/max(abs(aa(:, 2))), 1.05))
labelaxesfig(gca, "Beta / deg", "Intensity / a.u.")
saveName = strsplit(loadPath, '-');
saveas(gcf, append('../images/', saveName{3}, '-04-02-beta.png'))