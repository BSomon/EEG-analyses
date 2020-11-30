function IAF =Peak_IAF_fin(EEG,nFreq,sf,label)
%Measures the IAF from any dataset with the Pwelch method on each trial if
%epochs and on the whole dataset if raw signal. Uses a first window of
%7-13Hz to then determine the IAF either on occipital and/or parietal
%electrodes if the electrodes names and locations are provided and
%according to the 10-20 system; or on all the electrodes if there are no
%parietal or occipital electrodes.
% input :	- EEG: a matrix of the shape channels x time points x trials
% 			- nFreq: number of points over which the PSD estimate is
% 			calculated. Default: 2-second window
% 			- sf: EEG data sampling frequency
% 			- label: a cell containing the labels of the electrodes in the
% 			dataset. Default is based on Biosemi 64
% output: - IAF (individual alpha frequency)

% 14/09/20 -- B.S.
[nCh, nSa, nTr]=size(EEG);
if isempty(nFreq)
	nFreq = 1/2; 
end
%Gets the electrodes labels in case it is not provided, according to the
%default Biosemi electrode montage
if isempty(label)
	dChan = [{'Fp1'},{'AF7'},{'AF3'},{'F1'},{'F3'},{'F5'},{'F7'},{'FT7'},{'FC5'},...
		{'FC3'},{'FC1'},{'C1'},{'C3'},{'C5'},{'T7'},{'TP7'},{'CP5'},{'CP3'},...
		{'CP1'},{'P1'},{'P3'},{'P5'},{'P7'},{'P9'},{'PO7'},{'PO3'},{'O1'},{'Iz'},...
		{'Oz'},{'POz'},{'Pz'},{'CPz'},{'Fpz'},{'Fp2'},{'AF8'},{'AF4'},{'Afz'},...
		{'Fz'},{'F2'},{'F4'},{'F6'},{'F8'},{'FT8'},{'FC6'},{'FC4'},{'FC2'},...
		{'FCz'},{'Cz'},{'C2'},{'C4'},{'C6'},{'T8 (T4)'},{'TP8'},{'CP6'},{'CP4'},...
		{'CP2'},{'P2'},{'P4'},{'P6'},{'P8'},{'P10'},{'PO8'},{'PO4'},{'O2'}]; 
	label = dChan(1:nCh);
end

nSection = floor(sf/(nFreq*4.5)); %size of sections (in data points) into which the window is divided so that we have 8 segments with 50% overlap
vWindow = hamming(nSection);
% sWindow = length(vWindow); %size of the windows used to divide the signal into segments (with nOverlap as percent overlap)
nOverlap = floor(nSection/2);
nFFT = max(256, 2^nextpow2(nSection)); %50% overlap of the data for pwelch estimation

IAF = [];

IAFband = [7 14]; %Initiation of a large alpha frequency band

k=1;
vChan=[];

%Selects the occipital and parietal electrodes (for which the label
%contains either an O or a P)
for h=1:size(label,2)
        if strncmp(label{h},'P',1)||strncmp(label{h},'O',1)
            vChan(k)=h;
            k=k+1;
        end
end

%In case the system contains no parietal or occipital electrodes: use all
%channels
if length(vChan)<=1
    vChan = 1:size(label, 2);
end


chName=label(vChan);

eegspec = cell(length(vChan),1);
x=EEG;

for iChan = 1:length(vChan)% scan channels or component 
    for iTr = 1:nTr %Scan trials to perform the pwelch
        tmpdata = x(vChan(iChan),:,:); % channel activity
        [tmpspec,f] = pwelch(matsel(tmpdata,nSa,0,1,iTr),vWindow,nOverlap,nFFT,sf);
        eegspec{iChan} = cat(1, eegspec{iChan}, tmpspec');
    end
end

%Average PSD value over all trials
if nTr>1 
    eegspec = cellfun(@mean, eegspec, 'UniformOutput', false);
end

eegspec = cell2mat(eegspec);
SpecPower = eegspec(1:length(vChan),:)/nTr; %normalize by the number of sections
SpecPower=SpecPower';

%Finds the indeces of the frequency decomposition closest to the alpha band
%limits
[~, iFreq] = closest(f, IAFband);
MeanSpecPower = SpecPower(iFreq(1):iFreq(2),:,:);
vF = f(iFreq(1):iFreq(2),:,:);

% plot(f, SpecPower)
threshold=mean(MeanSpecPower);

chIAF=[];
for iChan=1:size(MeanSpecPower,2)
    [ampPeak, freqPeak] = Peak_localisation(MeanSpecPower(:,iChan), vF);
    %[m i]= max(PfmeanTr(:,ch));
    chIAF = cat(1,chIAF, [ampPeak, freqPeak]);
end

%Average IAF over all selected occipital and parietal channels (or all
%channels if this is the case)
IAFmean=nanmean(chIAF(:,2));
IAF=cat(2, IAF, round(IAFmean,1));
