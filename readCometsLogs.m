% Script for reading and plotting biomass, media, and flux logs from COMETS
% simulations involving two metabolic models. Calculates cross-correlation
% coefficients between metabolite uptake fluxes and biomass flux for
% cross-feeder model.
%
% Alan R. Pacheco 04.2024

clearvars
close all

%% Define file locations and metabolites of interest

% Directory where simulation results are contained
cometsDirectory = 'Simulation/';

% Minimal medium composition
mediumFile = 'Medium/minMed.mat';

% Format names for metabolites of interest
metDictMets = {'xylan4[e]';'ala__L[e]';'xyl__D[e]';'ac[e]';'pnto__R[e]';'nac[e]';'dtbt[e]';'btn[e]'};
metDictMetNames = {'Xylan';'L-alanine';'D-xylose';'Acetate';'Pantothenate';'Niacin';'Dethiobiotin';'Biotin'};

% Custom colors for strains
colorDictStrains = {'Leaf257','Leaf68'};
colorDictColors = [166,209,121;136,162,226]./255;

%% Load the medium and COMETS layout

load([cometsDirectory '/comets_layout.mat'])
load(mediumFile)

%% Process and plot biomass data

modelNames = cell(length(layout.models),1);
for m = 1:length(layout.models)
    modelNames{m} = layout.models{m}.modelID;
end

biomassLogRaw = parseBiomassLog([cometsDirectory '/' layout.params.biomassLogName]);

%% Read and plot biomass logs

biomassLogTotalPerStrain = zeros(layout.params.maxCycles,length(modelNames));
for t = 1:layout.params.maxCycles+1
    for m = 1:length(modelNames)
        biomassLogTotalPerStrain(t,m) = sum(biomassLogRaw.biomass(intersect(find(biomassLogRaw.t == t-1),find(biomassLogRaw.model == m-1))));
    end
end

figure
for m = 1:length(modelNames)
        plot([1:layout.params.maxCycles+1]*layout.params.timeStep,biomassLogTotalPerStrain(:,m),'LineWidth',3,'Color',colorDictColors(find(ismember(colorDictStrains,modelNames{m})),:))
    hold on
end
set(gca,'FontSize',20)
ylabel('Biomass (gDW)')
xlabel('Time (h)')
lgd = legend(modelNames);
lgd.FontSize = 16;
hold off

figure
relAbuDegrader = biomassLogTotalPerStrain(:,1)./sum(biomassLogTotalPerStrain,2);
plot([1:layout.params.maxCycles+1]*layout.params.timeStep,relAbuDegrader,'LineWidth',3,'Color',[243,86,88]./255)
set(gca,'FontSize',20,'YTick',[0:0.2:1])
ylabel(['Relative abundance, ' modelNames{1}])
xlabel('Time (h)')
xlim([0 round(max([1:layout.params.maxCycles+1]*layout.params.timeStep))])
ylim([-0.025,1.025])

%% Read and plot media logs

allMetsFromModels = layout.mets;
COMETSCycles = layout.params.maxCycles;

mediaLogRaw = parseMediaLog([cometsDirectory '/' layout.params.mediaLogName]);

metsOfInterest = setdiff(mediaLogRaw.metname(unique(mediaLogRaw.met(find(mediaLogRaw.amt > 1e-12)))),[minMed;{'co2[e]';'co[c]';'4hba[e]';'4hbz[e]';'amob[c]'}]);
metsOfInterestNames = metsOfInterest;
for i = 1:length(metsOfInterest)
    idx = find(ismember(metDictMets,metsOfInterest{i}));
    if ~isempty(idx)
        metsOfInterestNames(i) = metDictMetNames(idx);
    end
end

% Plot metabolite abundances per box
cmap = parula(length(metsOfInterest));
figure

for i = 1:length(metsOfInterest)
    mediaLogCurr = mediaLogRaw.amt(find(ismember(mediaLogRaw.metname, metsOfInterest{i})));
    mediaLogCurr(find(mediaLogCurr <= 1e-14)) = 1e-14;
    plot([1:layout.params.maxCycles+1]*layout.params.timeStep,mediaLogCurr,'Color',cmap(i,:),'LineWidth',4)
    hold on
end

set(gca,'FontSize',20,'yscale','log')
ylabel('Metabolite abundance (mmol)')
xlabel('Time (h)')
xlim([0 round(max([1:layout.params.maxCycles+1]*layout.params.timeStep))])
lgd = legend(metsOfInterestNames);
lgd.FontSize = 16;

%% Read and plot flux logs

load([cometsDirectory '/models.mat']) % Load the models

fluxesToTrack = mediaLogRaw.metname(unique(mediaLogRaw.met(find(mediaLogRaw.amt > 1e-12))));
fluxesToPlot = {'xylan4[e]';'ala__L[e]';'xyl__D[e]';'ac[e]';'pnto__R[e]';'nac[e]';'dtbt[e]'};
fluxesToTrack = unique([fluxesToTrack;fluxesToPlot;'Biomass']);
fluxesMat = zeros(length(fluxesToTrack),layout.params.maxCycles,length(modelNames));

%The log is a script that loads a struct named "fluxes"
run([cometsDirectory '/fluxLog.m'])

figure
for m = 1:length(modelNames)
    modelCurr = models.(modelNames{m});
    excRxnsCurr = find(findExcRxns(modelCurr));

    for f = 1:length(fluxesToTrack)
        if strcmp(fluxesToTrack{f},'Biomass')
            rxnIdx = find(ismember(modelCurr.rxns,'Growth'));
        else
            rxnIdx = intersect(excRxnsCurr,find(ismember(modelCurr.rxns,findRxnsFromMets(modelCurr,fluxesToTrack{f}))));
        end
        if ~isempty(rxnIdx)
            fluxVec = zeros(layout.params.maxCycles,1);
            if length(rxnIdx) == 2

                if strcmp(fluxesToTrack{f},'xyl__D[e]') % For xylose, plot only import
                    importersXylose = intersect(findRxnsFromMets(modelCurr,'xyl__D[e]'),findRxnsFromMets(modelCurr,'xyl__D[c]'));
                    importersXylose = [importersXylose;intersect(findRxnsFromMets(modelCurr,'xyl__D[e]'),findRxnsFromMets(modelCurr,'xyl__D[p]'))];
                    fluxVecsXylose = zeros(layout.params.maxCycles,length(importersXylose));
                    for x = 1:length(importersXylose)
                        rxnCoef = modelCurr.S(find(ismember(modelCurr.mets,'xyl__D[e]')),find(ismember(modelCurr.rxns,importersXylose{x})));
                        for t = layout.params.fluxLogRate:layout.params.fluxLogRate:length(fluxVec)
                            fluxVecsXylose(t,x) = fluxes{1,t}{1,1}{1,1}{1,m}(1,find(ismember(modelCurr.rxns,importersXylose{x})));
                        end
                        fluxVecsXylose(:,x) = fluxVecsXylose(:,x)./rxnCoef;
                    end
                    fluxVec = sum(fluxVecsXylose,2);
                else
                    error(['More than one exchange reaction for ' fluxesToTrack{f} ', ' modelNames{m}])
                end

            elseif length(rxnIdx) == 1
                for t = layout.params.fluxLogRate:layout.params.fluxLogRate:length(fluxVec)
                    fluxVec(t) = fluxes{1,t}{1,1}{1,1}{1,m}(1,rxnIdx);
                end
            else
                warning(['Could not plot flux for ' fluxesToTrack{f} ', ' modelNames{m}])
            end
            fluxesMat(f,:,m) = fluxVec;
        end
    end
end
fluxesMat(:,setdiff([1:1:length(fluxVec)],[layout.params.fluxLogRate:layout.params.fluxLogRate:length(fluxVec)]),:) = [];

for m = 1:length(modelNames)
    subplot(1,length(modelNames),m)
    plot([layout.params.fluxLogRate:layout.params.fluxLogRate:length(fluxVec)]*layout.params.timeStep,squeeze(fluxesMat(find(ismember(fluxesToTrack,fluxesToPlot)),:,m)),'LineWidth',2)
    set(gca,'FontSize',20)
    if m == 1; ylabel('Flux (mmol/gDW/hr)'); end
    xlabel('Time (h)')
    title(modelNames{m})
    if m == 1
        fluxesToPlotNames = fluxesToTrack(find(ismember(fluxesToTrack,fluxesToPlot)));
        for i = 1: length(fluxesToPlotNames)
            idx = find(ismember(metDictMets,fluxesToPlotNames{i}));
            if ~isempty(idx)
                fluxesToPlotNames(i) = metDictMetNames(idx);
            end
        end
        legend(fluxesToPlotNames)
    end
end

%% Compute and plot cross correlations between fluxes
maxLag = 500/(layout.params.timeStep*layout.params.fluxLogRate);
[xcorrMatBiomassMedia,lagsBiomassMedia,xcorrMatGRFluxes,lagsGRFluxes] = deal(NaN(length(metDictMets),length(modelNames),maxLag*2+1));
for m = 1:length(modelNames)
    
    biomassCurr = biomassLogTotalPerStrain([layout.params.fluxLogRate:layout.params.fluxLogRate:length(fluxVec)],m);
    growthRateCurr = fluxesMat(find(ismember(fluxesToTrack,'Biomass')),:,m)';
    
    for i = 1:length(metDictMets)
        
        fluxesCurr = fluxesMat(find(ismember(fluxesToTrack,metDictMets{i})),:,m)'*(-1);
        mediaCurrAll = mediaLogRaw.amt(find(ismember(mediaLogRaw.metname,metDictMets{i})));
        mediaCurr = mediaCurrAll([layout.params.fluxLogRate:layout.params.fluxLogRate:length(fluxVec)]);

        [xcorrCurr,lagsCurr] = xcorr(biomassCurr,mediaCurr,maxLag,'normalized');
        xcorrMatBiomassMedia(i,m,:) = xcorrCurr;
        lagsBiomassMedia(i,m,:) = lagsCurr;
        
        if ~isempty(fluxesCurr)
            [xcorrCurr,lagsCurr] = xcorr(growthRateCurr,fluxesCurr,maxLag,'normalized');
            xcorrMatGRFluxes(i,m,:) = xcorrCurr;
            lagsGRFluxes(i,m,:) = lagsCurr;
        end
    end
end

maxXCorrGRFluxesCrossFeeder = zeros(1,length(metDictMets));
for i = 1:length(metDictMets)
    [maxXCorrGRFluxesCrossFeeder(i),locMaxXCorrGRFluxesCrossFeeder(i)] = max(xcorrMatGRFluxes(i,2,:));
end

locMaxXCorrGRFluxesCrossFeeder = (locMaxXCorrGRFluxesCrossFeeder-maxLag-1)*layout.params.timeStep*layout.params.fluxLogRate;