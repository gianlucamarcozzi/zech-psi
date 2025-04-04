clearvars
addpath(genpath('util/'))

%% IMPORT

loadPath = '../data/processed/ZePSI-E-013031-ESEEM.mat';
load(loadPath)
nMeas = numel(y);hold on

for ii = 1:nMeas
    xAmp(ii) = Param{ii}.mpfuXAmp;
end

%% ESEEM INTEGRATION

% Find and plot integration window
integWidth = 100;
iMax = 190;
deltaTau = getparampulsespel(Param{1}, 'd30');  % ns
nMeas = numel(y);

% for ii = nMeas - 3:nMeas
for jj = 1:2
    figure(39 + jj)
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
        plot(sax, x{1}, x{2}, real(y{ii}));
        hold on
        plot(sax, x{1}, x{2}, imag(y{ii}));
        % plot(sax, x{1}, x{2}, real(integWindowPlot{ii}));
        plot(sax, x{1}, x{2}, imag(y2{ii}));
        % area(sax, x, y, 'FaceColor', 'b', 'FaceAlpha', 0.3);
        ylim(setaxlim([-1e4, 2e4], 1.1))
        plotTitle = strsplit(Param{ii}.TITL, '-');
        title(plotTitle{end})
        % titleStr = sprintf("%.1f deg, %.2f MPFU", ...
        %     Param{ii}.turningAngle/pi*180, xAmp(ii));
        % title(titleStr)
        % xline([-1, 1]/2*integWidth + x{1}(iMax), '--')
        % fprintf("%d max at:\t%d\n", [ii, iMax])
    end
end

% INTEGRATE TO GET ESEEM
eseem = zeros(nMeas, nTau);
d1 = getparampulsespel(Param{1}, 'd1 ');
xeseem = d1 + (0:deltaTau:(nTau - 1)*deltaTau);
for ii = 1:nMeas
    eseem(ii, :) = sum(y2{ii}(1:nTau, :), 2);
end

%%
I_BEST = [1, 5, 11];
figure(23)
clf

cmap = viridis(nMeas);
for ii = 1:nMeas
    % if ii > 11 && ii < 17
        % disp('Skipping some spectra')
    if false
        plot(xeseem, imag(eseem(ii, :)), 'x', 'Color', cmap(ii, :), ...
            'DisplayName', string(Param{ii}.turningAngle/pi*180))
    else
        plot(xeseem, imag(eseem(ii, :)), 'o-',  'Color', cmap(ii, :), ...
            'DisplayName', string(Param{ii}.turningAngle/pi*180))
    end
    hold on
end
xline(xeseem(I_BEST), "HandleVisibility", "off")
yline(0, "HandleVisibility", "off")
% plot(xeseem, imag(eseemcorr), 'Color', 'r', 'DisplayName', 'Corr')
legend('NumColumns', 3)

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
for ii = 1:numel(I_BEST)
    yplot = ybeta(:, ii)/max(abs(ybeta(:, ii)));
    yplot = -sign(yplot(1))*yplot;
    plot(xbeta*180/pi/1.075, yplot, 'o-', 'DisplayName', string(I_BEST(ii)))
    hold on
end
yline(0, 'HandleVisibility', 'off')
legend('Location', 'northwest')

aa = load("/home/gianluca/files/projects/oop-ciss-calculations/data/digitized/zech_p46_oopEseem_expData.csv");
plot(aa(:, 1), aa(:, 2)/max(abs(aa(:, 2))), 'kx-', 'DisplayName', 'Zech')

