
function[tempWF]=getAllSpikesWaveforms(st,mmf,clu)

%% this code extracts all spike waveforms 
%% size(tempWF)= number of spikes across all channels/clusters X 61 samples of 2 ms
channel_ID=tdfread('cluster_info.tsv');
wfWin = [-30:30];
tempWF=zeros(length(st),length(wfWin));
for num_spikes=1:length(st)-1
    clus_num=clu(num_spikes);
    channel_num=channel_ID.ch(find(channel_ID.cluster_id==clus_num));
    theseST = st(num_spikes); % spike times for cluster 19
    extractST = theseST;
  
    tempWF(num_spikes,:) = mmf.Data.x(channel_num+1,extractST+wfWin(1):extractST+wfWin(end));
end