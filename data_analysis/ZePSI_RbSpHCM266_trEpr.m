%%
clearvars

%%
loadFolder = "../data/raw/";
expNames = "ZePSI-E-" + ["014-002"];
nMeas = numel(expNames); 

for ii = 1:nMeas
    filename = loadFolder + expNames(ii);
    [x0, y{ii}, param{ii}] = eprload(filename);
    x{ii}.x1 = x0{1};
    x{ii}.x2 = x0{2};
end

figure(1)
clf
sax = ScrollableAxes();
for ii = 1:nMeas
    plot(sax, x{1}.x1, x{1}.x2, y{1}')
end
figure(2)
clf
sax = ScrollableAxes();
for ii = 1:nMeas
    plot(sax, x{ii}.x2, x{ii}.x1, y{ii})
end

