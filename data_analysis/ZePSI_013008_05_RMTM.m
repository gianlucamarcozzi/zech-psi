clearvars, clc %, close all

figure(1)
clf
tL = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for ii = 1:13
    filename = sprintf("../data/raw/ZePSI-E-013008/ZePSI-E-013008-%03u001", ii);
    yy = readtable(filename + ".Wfm.csv");
    yrm{ii} = yy{:, 1};
    ytm{ii} = yy{:, 2};
    Param = readtable(filename + ".csv");
    xx = linspace(Param{12, 2}, Param{13, 2}, Param{14, 2})*1e9;
    
    nexttile()
    plot(xx, yrm{ii})
    yyaxis right
    plot(xx, ytm{ii})
    yyaxis left
end

labelaxesfig(tL, 'Time / ns', '');
legendfirsttile(tL, 'RM', 'TM');
% savePath = '../images/RM_TM_pulseLengthDependence.png';

idx1 = 501:2750;
idx2 = (idx1(end) + 1):numel(xx);
if numel(idx1) ~= numel(idx2)
    error("idx1 and idx2 should be same length")
end
nx = numel(idx1);

figure(2)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
% Data correction
for ii = 1:13
    % Subtract constant background
    yrm{ii} = yrm{ii} - mean(yrm{ii}(1:500));  

    % Separate first pulse from second one
    p0{ii} = yrm{ii}(idx1);
    xx0 = xx(idx1);
    p1{ii} = yrm{ii}(idx2);
    xx1 = xx(idx2);
    
    nexttile
    plot(xx, yrm{ii}, xx0, p0{ii}, xx1, p1{ii})

end

%% FFT

tStep = xx(2) - xx(1);
fSampl = 1/tStep;
nzf = 4096;  % Zero filling
if nzf <= nx && nzf ~= 0
    warning("nzf <= nx. Continuing without zero-filling.")
    nzf = 0;
end
if nzf == 0
    ff = fSampl/(nx)*(-(nx)/2:(nx)/2 - 1);
else
    ff = fSampl/(nzf)*(-(nzf)/2:(nzf)/2 - 1);
end

ff = ff*1e3;  % MHz
figure(3)
clf
% tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")
for ii = 1:13

    if nzf ~= 0  % Zero filling
        p0{ii}(nzf) = 0;  % Zero filling
        p1{ii}(nzf) = 0;  % Zero filling
    else  % No zero filling (adjust values in memory in case)
        p0{ii} = p0{ii}(1:nx);
        p1{ii} = p1{ii}(1:nx);
    end

    fp0{ii} = fft(p0{ii});
    fp1{ii} = fft(p1{ii});
    
    ax1(ii) = nexttile(ii);    
    plot(ff, abs(fftshift(fp0{ii})), 'o-')
    hold on
    plot(ff, abs(fftshift(fp1{ii})), 'o-')
    xlim(ax1(ii), setaxlim(ff, 0.05))
    
    % ax1(ii).TickDir = 'out';
    % ax1(ii).XAxis(2).TickLabels = [];
    % ax1(ii).XAxis(2).TickLength = [0 0];
    set(gca,'box','off')
end

% Secondary x-axis
for ii = 1:13
    ax2(ii) = axes('Position', ax1(ii).Position, 'XAxisLocation', 'top', ...
        'YAxisLocation', 'right', 'Color', 'none');
    ffb = mhz2mt(ff)*10;
    xlim(ax2(ii), setaxlim(ffb, 1))
    yticks(ax2(ii) , [])
    xlim(ax2(ii), setaxlim(ffb, 0.05))
    xline([-10, 10])

    % line(ff2, abs(fftshift(fp0{ii})), 'Parent', ax2(ii), 'Color', 'none');
    % 
    % ax2(ii).XColor = 'r';
    % ax2(ii).XLabel.String = 'Secondary X-axis';
    % ax2(ii).XLabel.Color = 'r';
end

%{
% Create secondary x-axis
ax1 = gca;
ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', 'Color', 'none');
x2 = x1 * 2; % Example conversion for secondary x-axis values
line(x2, y1, 'Parent', ax2, 'Color', 'none');

% Link the two x-axes
linkaxes([ax1 ax2], 'x');

% Customize secondary x-axis appearance
ax2.XColor = 'r';
ax2.XLabel.String = 'Secondary X-axis';
ax2.XLabel.Color = 'r';

% Ensure primary axes are active for further modifications
axes(ax1);
%}
