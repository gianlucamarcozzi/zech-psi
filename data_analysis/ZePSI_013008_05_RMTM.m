clearvars, clc %, close all

figure()
clf
tL = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for ii = 1:13
    filename = sprintf("../data/raw/ZePSI-E-013008/ZePSI-E-013008-%03u001", ii);
    yy = readtable(filename + ".Wfm.csv");
    Param = readtable(filename + ".csv");
    xx = linspace(Param{12, 2}, Param{13, 2}, Param{14, 2})*1e9;
    
    nexttile()
    plot(xx, yy{:, 1})
    yyaxis right
    plot(xx, yy{:, 2})
    yyaxis left
end

labelaxesfig(tL, 'Time / ns', '');
legendfirsttile(tL, 'RM', 'TM');
% savePath = '../images/RM_TM_pulseLengthDependence.png';

