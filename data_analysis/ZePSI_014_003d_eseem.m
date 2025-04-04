clearvars
addpath(genpath('util/'))

%% IMPORT
loadPath = '../data/processed/ZePSI-E-014003-ESEEM.mat';
load(loadPath)
nMeas = numel(y);

% Find and plot integration window
integWidth = 120;
iMax = 196;
deltaTau = getparampulsespel(param{1}, 'd30');  % ns
nMeas = numel(y);

% for ii = nMeas - 3:nMeas
for jj = 1:2
    figure(999 + jj)
    clf
    tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
    ii1 = round((jj - 1)*nMeas/2 + 1);
    ii2 = round((jj)*nMeas/2);
    for ii = ii1:ii2
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
    
        nexttile
        sax = ScrollableAxes();
        plot(sax, x.x1, x.x2, real(y{ii}));
        hold on
        plot(sax, x.x1, x.x2, imag(y{ii}));
        plot(sax, x.x1, x.x2, imag(y2{ii}));
        ylim(setaxlim([-1e4, 2e4], 1.1))
        plotTitle = strsplit(param{ii}.TITL, '-');
        title(plotTitle{end}(1:3))
    end
end

% INTEGRATE TO GET ESEEM
eseem = zeros(nMeas, nTau);
d1 = getparampulsespel(param{1}, 'd1 ');
xeseem = d1 + (0:deltaTau:(nTau - 1)*deltaTau);
for ii = 1:nMeas
    eseem(ii, :) = sum(y2{ii}(1:nTau, :), 2);
end

%%
I_BEST = [31:33];
figure(23)
clf

cmap = viridis(nMeas);
for ii = 1:nMeas
    displayName = sprintf("%d deg", round(param{ii}.turningAngle/pi*180));
    plot(xeseem, imag(eseem(ii, :)), '-',  'Color', cmap(ii, :), ...
        'DisplayName', displayName)
    hold on
end
xline(xeseem(I_BEST), "HandleVisibility", "off")
xlim(setaxlim(xeseem, 1))
ylim(setaxlim([min(imag(eseem(:))), max(imag(eseem(:)))], 1.05))
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
        % ybeta(ii, jj) = yfin(ii, I_BEST(jj));
    end
    xbeta(ii) = param{ii}.turningAngle;    
end

figure(33)
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