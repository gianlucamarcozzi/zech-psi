%%
clearvars

%%
loadFolder = "../data/raw/";
expNames = "ZePSI-E-" + ["014-001"];
nMeas = numel(expNames); 

for ii = 1:nMeas
    filename = loadFolder + expNames(ii);
    [x{ii}, y0, param{ii}] = eprload(filename);
    y{ii} = y0{1} + 1i*y0{2};
end

figure(1)
clf
for ii = 1:nMeas
    nexttile
    plot(x{ii}, real(y{ii}))
    hold on
    plot(x{ii}, imag(y{ii}))
    xlim(setaxlim(x{ii}, 1))
    ylim(setaxlim(real(y{ii}), 1.05))
end
legend("Real", "Imag")
labelaxesfig(gca, "Magnetic Field / G", "Intensity / a.u.")