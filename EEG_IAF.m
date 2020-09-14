function [EEG, vFreq] = EEG_IAF(EEG, tIAF);
IAFt = [];
sf = EEG.srate;
EEGdata = EEG.data;
NFFT = 2*sf; %2-second window for computing pwelch

if isempty(EEG.chanlocs)
	channels = [];
else
	channels = {EEG.chanlocs(:).labels};
end

if size(EEGdata,2)<(tIAF*sf)
	disp(['You are not providing enough data. Your dataset must be at least '...
		num2str(tIAF*sf), ' points long.']);
	return;
else
		IAFt = cat(1, IAFt, Peak_IAF_fin(EEGdata(:, tIAF*sf,:),sf/NFFT,sf,channels));
		iIAF = round(nanmean(IAFt));
end
EEG.vTheta = [iIAF-6, iIAF-2]; % Individual theta frequency band
EEG.vAlpha = [iIAF-2, iIAF+2]; % Individual alpha frequency band
EEG.vBeta = [iIAF+3, iIAF+20]; % Individual beta frequency band
vFreq = cat(1,vTheta, vAlpha, vBeta);
end