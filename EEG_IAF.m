%Main function associated to the individual alpha frequency definition
%(Peak_IAF_fin)
%Inputs:	- EEG: an EEGLAB EEG structure containing at least the srate, data
%				and chanlocs (even if empty) fields.
%			- tIAF: the time period over which you want to compute the
%				individual alpha frequency measure - at least 60 seconds is
%				ideal
%Outputs:	- EEG: an EEGLAB EEG structure with 3 new fields containing the
%				frequency bands for theta alpha and beta based on the IAF
%				definition
%			- vFreq: a vector containing the frequency bands limits for
%			thetas, alpha and beta concatenated over lines

% 14/09/20 -- B.S.

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
		IAFt = cat(1, IAFt, Peak_IAF_fin(EEGdata(:, 1:tIAF*sf,:),sf/NFFT,sf,channels));
		iIAF = round(nanmean(IAFt));
end
EEG.vTheta = [iIAF-6, iIAF-2]; % Individual theta frequency band
EEG.vAlpha = [iIAF-2, iIAF+2]; % Individual alpha frequency band
EEG.vBeta = [iIAF+3, iIAF+20]; % Individual beta frequency band
vFreq = cat(1,vTheta, vAlpha, vBeta);
end
