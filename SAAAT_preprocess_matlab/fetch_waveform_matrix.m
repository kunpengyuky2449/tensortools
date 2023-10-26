clear
%%
listing = dir('Y:\Users\hsrivastava\Ephys\AE_trained_Behaving');
listing = listing(~ismember({listing.name},{'.','..'}));
folder_path_EphysRaw = cell(0);
for i = 1:numel(listing)
    sub_listing = dir(fullfile(listing(i).folder,listing(i).name));
    sub_listing = sub_listing(~ismember({sub_listing.name},{'.','..'}));
    for j = 1:numel(sub_listing)
        subsub_listing = dir(fullfile(sub_listing(j).folder,sub_listing(j).name));
        subsub_listing = subsub_listing(~ismember({subsub_listing.name},{'.','..','__pycache__'}));
        destination = subsub_listing( contains({subsub_listing.name},"W3") & ...
            [subsub_listing.isdir]).name; 
        folder_path_EphysRaw{end+1} = fullfile(subsub_listing(1).folder, destination);
    end
end
%%
for i = 1:numel(folder_path_EphysRaw)
    target_folder = folder_path_EphysRaw{i};
    save_folder = 'D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\SpikeAnalysis\Spikes_Waveform_Matrix';
    cd(target_folder)
    %%
    clusterID = readNPY('spike_clusters.npy');
    pp='amplifier.dat';
    nChansInFile = 128;  % neuropixels phase3a, from spikeGLX
    d = dir(pp);
    nSamps = d.bytes/2/nChansInFile;
    mmf = memmapfile(pp, 'Format', {'int16', [nChansInFile nSamps], 'x'});
    st = readNPY('spike_times.npy');
    
    channel_ID=tdfread('cluster_info.tsv');
    wfWin = -30:30;
    tempWF=zeros(length(st),length(wfWin));
    wbar = waitbar(0,'Start Processing...');
    for num_spikes=1:length(st)-1
        clus_num=clusterID(num_spikes);
        channel_num=channel_ID.ch(find(channel_ID.cluster_id==clus_num));
        theseST = st(num_spikes); % spike times for cluster 19
        extractST = theseST;
      
        tempWF(num_spikes,:) = mmf.Data.x(channel_num+1,extractST+wfWin(1):extractST+wfWin(end));
        if rem(num_spikes,1000) == 0
            temp_progress = (num_spikes-1)/numel(st);
            waitbar(temp_progress, wbar, strcat('Progress:',num2str(temp_progress*100),'%'));
        end
    end
    close(wbar)
    %%
    temp2 = strfind(target_folder,'\');
    file_name = strcat(target_folder(temp2(end)+1:temp2(end)+8), ...
        'waveform_matrix_', datestr(now,'yymmdd_HH_MM_SS'));
    KSLabel = channel_ID.KSLabel;
    save(fullfile(save_folder, file_name),"tempWF","KSLabel")
    disp(strcat(num2str(i),' Sessions have been done'));
end