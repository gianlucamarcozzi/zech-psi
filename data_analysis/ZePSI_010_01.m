% Recording transients of oop-eseem to subtract light minus dark
clearvars

%% Import
expName = "../data/raw/ZePSI-E-010";
measNo = 1:2;
nMeas = numel(measNo);

for jj = 1:nMeas
    filename = expName + "00" + string(measNo(jj)) + ".DTA";
    [xRaw{jj}, yRaw{jj}, Param{jj}] = eprload(filename);
    xRaw{jj}{1} = xRaw{jj}{1}';
end

[nt, ntrans] = size(yRaw{1});
yRingdown = yRaw{1}(:, end);
x = xRaw;
for jj = 1:nMeas
    y1{jj} = yRaw{jj}; % - repmat(yRingdown, [1, ntrans]);
end
yMinus = yRaw{2} - yRaw{1};

%% Scrollable traces

iPlot1 = 1;
iPlot2 = 2;
clf
h = ScrollableAxes();
% plot(h, x{iPlot1}{1}', x{iPlot1}{2}', real(y1{iPlot1})');
% hold on
% plot(h, x{iPlot2}{1}', x{iPlot2}{2}', real(y1{iPlot2})');
plot(h, x{iPlot1}{1}', x{iPlot1}{2}', imag(yMinus)');

%% Phase correction

model = @(y, phi) y*exp(1i*phi*pi/180);
% for ii = 1:nMeas
%     [y{ii}, phase{ii}] = correctphase(yRaw{2}, 'Maximum');
% end

phis = 0:180;
for iphi = 1:numel(phis)
    yphi(:, iphi) = model(y1{2}(:, 25), phis(iphi));
end

% figure()
clf
h = ScrollableAxes();
plot(h, x{iPlot1}{1}', phis', real(yphi)');
hold on
plot(h, x{iPlot1}{1}', phis', imag(yphi)');


iPlot1 = 1;
iPlot2 = 2;
% clf
h = ScrollableAxes();
% plot(h, x{iPlot1}{1}', x{iPlot1}{2}', imag(model(y1{iPlot1}, 0))');
% hold on
% plot(h, x{iPlot2}{1}', x{iPlot2}{2}', imag(model(y1{iPlot2}, 0))');
% plot(h, x{iPlot1}{1}', x{iPlot1}{2}', imag(yMinus)');


%%

% modelFit = @(p) p(1)* (gaussianresonancebsweep(xa, p(1), p(2));
p0 = [100, 30];

xa = x{1}{1};
ya = yMinus(:, 1);
ya = ya/max(abs(ya));
clf
plot(xa, imag(ya))
hold on
plot(xa, modelFit(p0))

pfit = nlinfit(xa, ya, modelFit, p0);
yFit = modelFit(pfit);

plot(xa, yFit)

%%
function y1 = gaussianresonancebsweep(wReson, mwFreq, c, mode)
    arguments
        wReson
        mwFreq  
        c    (1, 1) double     
        mode string = "lwpp"
    end    

    normTerm = sqrt(2/pi)/c;
    y1 = normTerm*exp(-2 * (wReson - mwFreq).^2 ./ c^2);

end


