%% This script allows to extract and plot ERPs, time-frequency and frequency data at the study level on EEGLAB 2019.x
% Prerequisites are to load a STUDY on EEGLAB and pre-compute the ERPs, TF
% and Spec data. When launching the script, it will ask you what electrodes
% you want to compute your data at, and whether you wish to equalize the
% number of trials in all conditions for frequecy and time-frequency data
% analyses (If you do not know why this can be relevant, I deeply enourage
% you to read Mike X. Cohen's book "Analyzing Neural Time-Series Data"). 
% You can launch each part (Spectral analysis, ERSP analysis, and ITC analysis)
% separately. Don't forget to launch however the first lines for the frequency
% range, channels and trial equalization.
% If you have any difficulty, contact me: bertille.somon@isae-supaero.fr

% !!! If you wish to obtain the exact p-values (pcond_ex, pinter_ex, pgroup_ex),
% you need to replace the std_stat() function in EEGLAB 
% (eeglab\functions\studyfunc\std_stat.m) by the one provided in this repository.

% 27/08/20 - BS

iChan = input('Which electrode do you wish to analyze? (Provide the EEGLAB electrode number)\n');
equalize = input('Do you wish to perform data equalization? \n 0: No \n 1: Yes \n');
FreqRange = input('What frequency range od you wish to analyze?');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%SPECTRAL DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recovering the same number of trials for participants in the odd and
%standard condition
s = RandStream('mlfg6331_64');
singletrials = 'on';

STUDY = pop_specparams(STUDY, 'freqrange', FreqRange, 'plotgroups', 'together',...
	'plotconditions', 'together', 'subtractsubjectmean', 'on', 'averagechan', 'on');

[~, datainit] = std_readdata(STUDY, ALLEEG, 'datatype', 'spec',...
	'subject', STUDY.design(STUDY.currentdesign).cases.value{1}, 'singletrials',...
	'on', 'channels', {STUDY.changrp(1).name}, 'rmsubjmean', STUDY.etc.specparams.subtractsubjectmean);
specdata=cell(size(datainit));
subj_by_subj_spec = specdata;

%Selects (randomly) the same number of trials in every condition at the level of the participant
for iSubj = 1:length(STUDY.design(STUDY.currentdesign).cases.value)
	[STUDY, specdata_subj, specfreqs_subj] = std_readdata(STUDY, ALLEEG, 'datatype', 'spec',...
		'subject', STUDY.design(STUDY.currentdesign).cases.value{iSubj}, 'singletrials',...
		singletrials, 'channels', {STUDY.changrp(iChan).name}, 'rmsubjmean', STUDY.etc.specparams.subtractsubjectmean);
	if equalize
		nMin = min(cell2mat(cellfun(@(x)size(x,ndims(x)), specdata_subj, 'UniformOutput', false)), [], 'all');
		iTrials = cellfun(@(x)randsample(s, [1:size(x,ndims(x))], nMin), specdata_subj, 'UniformOutput', false);
		specdata = cellfun(@(x,y)cat(ndims(y)+1,x,y), specdata, cellfun(@(x,y)mean(x(:,:,y),ndims(x)),...
			specdata_subj, iTrials, 'UniformOutput', false), 'UniformOutput', false);
		subj_by_subj_spec = cellfun(@(x, y, z)cat(ndims(y), x, y(:,:,z)), subj_by_subj_spec,...
			specdata_subj, iTrials, 'UniformOutput', false);
	else
		specdata = cellfun(@(x,y)cat(ndims(y)+1,x,y), specdata, cellfun(@(x)mean(x,ndims(x)),...
			specdata_subj, 'UniformOutput', false), 'UniformOutput', false);
		subj_by_subj_spec = cellfun(@(x, y)cat(ndims(y), x, y), subj_by_subj_spec,...
			specdata_subj, 'UniformOutput', false);
	end
	if strcmpi(STUDY.etc.specparams.subtractsubjectmean, 'on')
		count    = 0;
		for iSpec = 1:length(specdata(:))
			if ~isempty(specdata{iSpec})
				if count == 0, meanspec = zeros(size(specdata{iSpec},1),1); end
				if strcmpi(singletrials, 'on')
					meanspec = meanspec + mean(specdata{iSpec},2);
				else
					meanspec = meanspec + specdata{iSpec};
				end
				count = count+1;
			end
		end
		meanspec = meanspec/count;
		specdata  = cellfun(@(x)bsxfun(@minus, x, meanspec), specdata, 'uniformoutput', false);
	end
end

specfreqs = specfreqs_subj;

if size(specdata, 2) > 1
	specdata_chan_chan = specdata;
	specdata = cellfun(@(x)squeeze(mean(x,2)), specdata, 'UniformOutput', false);
end

STUDY = pop_specparams(STUDY, 'freqrange', FreqRange, 'plotgroups', 'together',...
	'plotconditions', 'together', 'subtractsubjectmean', 'on', 'averagechan', 'on');

%If you want the exact p-values
[pcond, pgroup, pinter, statscond, statsgroup, statsinter, pcond_ex, pgroup_ex, pinter_ex]...
	= std_stat(specdata, 'groupstats', 'on', 'condstats', 'on', 'method', 'permutation',...
	'alpha', 0.05, 'mcorrect', 'fdr', 'mode', 'eeglab', 'paired', {'on' 'on'});
%Otherwise use the following line:
%[pcond, pgroup, pinter, statscond, statsgroup, statsinter, pcond_ex, pgroup_ex, pinter_ex]...
%	= std_stat(specdata, 'groupstats', 'on', 'condstats', 'on', 'method', 'permutation',...
%	'alpha', 0.05, 'mcorrect', 'fdr', 'mode', 'eeglab', 'paired', {'on' 'on'});


std_plotcurve(specfreqs, specdata, 'groupstats', pgroup, 'legend', 'on',...
	'condstats', pcond, 'interstats', pinter, 'plotsubjects', 'off',...
	'plotstderr', 'on', 'threshold',0.05, 'unitx', 'Hz', 'filter',[], 'plotgroups','together', ...
	'effect', 'main', 'plotconditions', 'together');
std_plotcurve(specfreqs, specdata, 'legend', 'on', 'plotsubjects', 'off',...
	'plotstderr', 'on', 'condnames', STUDY.design(STUDY.currentdesign).variable(1).value,...
	'groupnames', STUDY.design(STUDY.currentdesign).variable(2).value, ...
	'threshold',0.05, 'unitx', 'Hz', 'filter',[], 'plotgroups','together', ...
	'effect', 'main', 'plotconditions', 'together');


%Saves statistical data in a variable
stat_spec = cell(3,3);
stat_spec = {pcond, pgroup, pinter; statscond, statsgroup, statsinter; pcond_ex, pgroup_ex, pinter_ex};

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TIME-FREQUENCY DATA%%%%%%%%%%%%%%%%%%%%%%%
s = RandStream('mlfg6331_64');
singletrials = 'on';
[~, datainit] = std_readdata(STUDY, ALLEEG, 'datatype', 'ersp',...
	'subject', STUDY.design(STUDY.currentdesign).cases.value{1}, 'singletrials',...
	singletrials, 'channels', {STUDY.changrp(iChan).name});
erspdata = cell(size(datainit));
subj_by_subj_ersp = erspdata;

for iSubj = 1:length(STUDY.design(STUDY.currentdesign).cases.value)
	[STUDY, erspdata_subj, ersptimes_subj, erspfreqs_subj, ~, paramsersp] =...
		std_readdata(STUDY, ALLEEG, 'datatype', 'ersp',	'subject',...
		STUDY.design(STUDY.currentdesign).cases.value{iSubj}, 'singletrials',...
		singletrials, 'channels', {STUDY.changrp(iChan).name});
	
	if equalize
		nMin = min(cell2mat(cellfun(@(x)size(x,ndims(x)), erspdata_subj, 'UniformOutput', false)), [], 'all');
		iTrials = cellfun(@(x)randsample(s, [1:size(x,ndims(x))], nMin), erspdata_subj, 'UniformOutput', false);
		erspdata = cellfun(@(x,y)cat(ndims(y)+1,x,y), erspdata, cellfun(@(x,y)mean(x(:,:,:,y),ndims(x)),...
			erspdata_subj, iTrials, 'UniformOutput', false), 'UniformOutput', false);
		subj_by_subj_ersp = cellfun(@(x, y, z)cat(ndims(y), x, y(:,:,:,z)), subj_by_subj_ersp,...
			erspdata_subj, iTrials, 'UniformOutput', false);
	else
		erspdata = cellfun(@(x,y)cat(ndims(y)+1,x,y), erspdata, cellfun(@(x)mean(x,ndims(x)),...
			erspdata_subj, 'UniformOutput', false), 'UniformOutput', false);
		subj_by_subj_ersp = cellfun(@(x, y)cat(ndims(y), x, y), subj_by_subj_ersp,...
			erspdata_subj, 'UniformOutput', false);
	end
	
end

ersptimes = ersptimes_subj;
erspfreqs = erspfreqs_subj;
paramsersp.singletrials = 'off';
paramsersp.commonbase = 'off';

STUDY = pop_erspparams(STUDY, 'timerange', [], 'freqrange', FreqRange, 'ersplim', [],...
	'subbaseline', paramsersp.commonbase, 'singletrials', 'off', 'averagemode', 'ave', 'averagechan', 'on');


erspdatabln = newtimefbaseln(erspdata, ersptimes, paramsersp);
erspdatabln = cellfun(@(x)squeeze(mean(x,3)), erspdatabln, 'UniformOutput', false);

%If you want the exact p-values
[pcond, pgroup, pinter, statscond, statsgroup, statsinter, pcond_ex, pgroup_ex, pinter_ex]...
	= std_stat(erspdatabln, 'groupstats', 'on', 'condstats', 'on', 'method', 'permutation',...
	'alpha', 0.05, 'mcorrect', 'fdr', 'mode', 'eeglab', 'paired', {'on' 'on'});
%Otherwise use the following line:
%[pcond, pgroup, pinter, statscond, statsgroup, statsinter, pcond_ex, pgroup_ex, pinter_ex]...
%	= std_stat(specdata, 'groupstats', 'on', 'condstats', 'on', 'method', 'permutation',...
%	'alpha', 0.05, 'mcorrect', 'fdr', 'mode', 'eeglab', 'paired', {'on' 'on'});

std_plottf(ersptimes, erspfreqs, erspdatabln, 'datatype', 'ersp', ...
	'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'plottopo', 'off', 'plotmode', ...
	'normal', 'unitcolor', 'dB', 'ersplim', [], 'threshold',   .05, ...
	'effect', 'main', 'maskdata', 'off', 'averagemode', 'ave');

%Topographies
locsOri = eeg_mergelocs(ALLEEG.chanlocs);
locsname = {locsOri(iChan).labels};
locs = locsOri(std_chaninds(STUDY, locsname));

Freqs = ('Frequency range to display topographies? (as a vector [min max])');
Times = ('Time range to display topographies? (as a vector [min max in ms])');

[~,idxt1] = closest(ersptimes,min(Times));
[~,idxt2] = closest(ersptimes,max(Times));
[~,idxf1] = closest(erspfreqs,min(Freq));
[~,idxf2] = closest(erspfreqs,max(Freq));
for index = 1:length(erspdatabln(:))
	erspdatabln{index} = mean(mean(erspdatabln{index}(idxf1:idxf2,idxt1:idxt2,:,:),1),2);
	erspdatabln{index} = reshape(erspdatabln{index}, [1 size(erspdatabln{index},3) size(erspdatabln{index},4) ]);
end

std_chantopo(erspdatabln,'chanlocs', locs);

stat_ersp = cell(3,3);
stat_ersp = {pcond, pgroup, pinter; statscond, statsgroup, statsinter; pcond_ex, pgroup_ex, pinter_ex};


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ITC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = RandStream('mlfg6331_64');
singletrials = 'on';
itctype = 'phasecoher';
[~, datainit] = std_readdata(STUDY, ALLEEG, 'datatype', 'itc',...
	'subject', STUDY.design(STUDY.currentdesign).cases.value{1}, 'singletrials',...
	singletrials, 'channels', {STUDY.changrp(iChan).name});
itcdata = cell(size(datainit));
subj_by_subj_itc = itcdata;
for iSubj = 1:length(STUDY.design(STUDY.currentdesign).cases.value)
	[STUDY, itcdata_subj, itctimes_subj, itcfreqs_subj] = std_readdata(STUDY, ALLEEG, 'datatype', 'itc',...
		'subject', STUDY.design(STUDY.currentdesign).cases.value{iSubj}, 'singletrials',...
		singletrials, 'channels', {STUDY.changrp(iChan).name});
	if equalize
		nMin = min(cell2mat(cellfun(@(x)size(x,ndims(x)), itcdata_subj, 'UniformOutput', false)), [], 'all');
		iTrials = cellfun(@(x)randsample(s, [1:size(x,ndims(x))], nMin), itcdata_subj, 'UniformOutput', false);
		if length(iChan) > 1
			itcdata = cellfun(@(x,y)cat(ndims(y)+1,x,y), itcdata, cellfun(@(x,y)abs(permute(newtimefitc(permute(...
				x(:,:,:,y),[3 1 2 4]), itctype),[2 3 1])),...
				itcdata_subj, iTrials, 'UniformOutput', false), 'UniformOutput', false);
			subj_by_subj_itc = cellfun(@(x, y, z)cat(ndims(y), x, abs(y(:,:,:,z))), subj_by_subj_itc,...
				itcdata_subj, iTrials, 'UniformOutput', false);
		else
			itcdata = cellfun(@(x,y)cat(ndims(y)+1,x,y), itcdata, cellfun(@(x,y)abs(newtimefitc(x(:,:,y), itctype)),...
				itcdata_subj, iTrials, 'UniformOutput', false), 'UniformOutput', false);
			subj_by_subj_itc = cellfun(@(x, y, z)cat(ndims(y), x, abs(y(:,:,z))), subj_by_subj_itc,...
				itcdata_subj, iTrials, 'UniformOutput', false);
		end
	else
		if length(iChan) > 1
			itcdata = cellfun(@(x,y)cat(ndims(y)+1,x,y), itcdata, cellfun(@(x)abs(permute(newtimefitc(permute(...
				x,[3 1 2 4]), itctype),[2 3 1])),...
				itcdata_subj, 'UniformOutput', false), 'UniformOutput', false);
			subj_by_subj_itc = cellfun(@(x, y)cat(ndims(y), x, abs(y)), subj_by_subj_itc,...
				itcdata_subj, 'UniformOutput', false);
		else
			itcdata = cellfun(@(x,y)cat(ndims(y)+1,x,y), itcdata, cellfun(@(x)abs(newtimefitc(x, itctype)),...
				itcdata_subj, 'UniformOutput', false), 'UniformOutput', false);
			subj_by_subj_itc = cellfun(@(x, y)cat(ndims(y), x, abs(y)), subj_by_subj_itc,...
				itcdata_subj, 'UniformOutput', false);
		end
	end
end

itctimes = itctimes_subj;
itcfreqs = itcfreqs_subj;
paramsersp.singletrials = 'off';
paramsersp.commonbase = 'on';
STUDY = pop_erspparams(STUDY, 'timerange', [], 'freqrange', FreqRange, 'ersplim', [],...
	'subbaseline', paramsersp.commonbase, 'singletrials', 'off', 'averagemode', 'ave', 'averagechan', 'on');

%If you want the exact p-values
[pcond, pgroup, pinter, statscond, statsgroup, statsinter, pcond_ex, pgroup_ex, pinter_ex]...
	= std_stat(itcdata, 'groupstats', 'on', 'condstats', 'on', 'method', 'permutation',...
	'alpha', 0.05, 'mcorrect', 'fdr', 'mode', 'eeglab', 'paired', {'on' 'on'});
%Otherwise use the following line:
%[pcond, pgroup, pinter, statscond, statsgroup, statsinter, pcond_ex, pgroup_ex, pinter_ex]...
%	= std_stat(specdata, 'groupstats', 'on', 'condstats', 'on', 'method', 'permutation',...
%	'alpha', 0.05, 'mcorrect', 'fdr', 'mode', 'eeglab', 'paired', {'on' 'on'});

std_plottf(itctimes, itcfreqs, itcdata, 'datatype', 'itc', ...
	'groupstats', pgroup, 'condstats', pcond, 'interstats', pinter, 'plottopo', 'off', 'plotmode', ...
	'normal', 'unitcolor', 'dB', 'ersplim', [], 'threshold',   .05, ...
	'effect', 'main', 'maskdata', 'off', 'averagemode', 'ave');

stat_itc = cell(3,3);
stat_itc = {pcond, pgroup, pinter; statscond, statsgroup, statsinter; pcond_ex, pgroup_ex, pinter_ex};


save('Data_equal.mat', 'subj_by_subj_spec', 'specdata', 'spectimes', 'stat_spec',...
	'subj_by_subj_ersp', 'erspdata','ersptimes', 'erspfreqs', 'stat_ersp',...
	'subj_by_subj_itc', 'itcdata', 'itctimes', 'itcfreqs', 'stat_itc');
