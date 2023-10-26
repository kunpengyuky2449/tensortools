clear
%%
load('SAAAT_session_summary_073123.mat')
%wbar = waitbar(0,'Start plotting...');
listing_waveform = dir("Y:\Users\KunpengY\Waveform_matrix2");
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
        folder_path_EphysRaw{end+1} = fullfile(subsub_listing(1).folder, destination); %#ok
    end
end

%%
save_folder = 'D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\neuron_active_window_test_with_drift_and_waveform';
if ~exist(save_folder, 'dir')
   mkdir(save_folder)
end
wbar = waitbar(0,'Start Processing...');
for i = 17:numel(session_summary.file_path_PARAMS) %!
    waitbar((i-1)/numel(session_summary.file_path_PARAMS), ...
                wbar,strcat('Processing session:',num2str(i)));
    load(session_summary.file_path_Spikes{i})
    load(session_summary.file_path_PARAMS{i})
    Fs = PARAMS.Fs;
    ds_sp_time = round(sp_times/Fs*10);
    clear Target_Start TC_Start sp_times
    target_folder = folder_path_EphysRaw{i};
    cd(target_folder)
    raw_clusterID = readNPY('spike_clusters.npy');
    raw_st = readNPY('spike_times.npy')/Fs;
    load("rez.mat","rez");
    batch_samples = rez.ops.NT;
    dshift = rez.dshift;
    dshift_time = (1:height(dshift))'* batch_samples/Fs/60;
    spike_clusters=spike_clusters(1:numel(ds_sp_time));
    uniq_neurons = unique(spike_clusters);
    temp = strfind(session_summary.file_path_Spikes{i},'Spikes');
    session_name = session_summary.file_path_Spikes{i}(temp-8:temp-1);
    temp = contains({listing_waveform.name},session_name);
    load(fullfile(listing_waveform(temp).folder,listing_waveform(temp).name))
    clear rez PARAMS
    for n = 1:numel(uniq_neurons)
        close
        % Filter session activity and fetch good sig window
        ds_sp_time_iit = ds_sp_time(spike_clusters == uniq_neurons(n));
        firing_rate_array = zeros(91*60*10,1);
        [N_spike,time_bin] = groupcounts(ds_sp_time_iit);
        firing_rate_array(time_bin,1)=...
            firing_rate_array(time_bin,1) + N_spike;
        % L_gaus = 600; %unit: 0.1s, decide number of sample per filter to covers
        % a_gaus = 4; %unit: standard deviation, decide how wide the distribution spread, better keep this same while change L_gaus
        % w_gaus = gausswin(L_gaus, a_gaus);
        %wvtool(w_gaus);
        windowSize = 150; 
        w_ma = (1/windowSize)*ones(1,windowSize);
        firing_rate_array = filter(w_ma, 1, firing_rate_array);
        FR_sorted = sort(firing_rate_array);
        FR_CSum = cumsum(FR_sorted);
        [a1, a2] = min(abs(FR_CSum-0.05*max(FR_CSum)));
        threshold = FR_sorted(a2)/2;
        %disp(threshold)
        index_temp =  (firing_rate_array>=threshold);
        index_goodsig = zeros(size(index_temp));
        for window = 1:numel(index_temp)-600
            if sum(index_temp(window:window+600)) > 480
                index_goodsig(window:window+600) = 1;
            end
        end

        % Get concatenated waveform over time
        this_WF = tempWF(raw_clusterID == uniq_neurons(n), :);
        this_st = raw_st(raw_clusterID == uniq_neurons(n), :);
        sampling_windows = [5,15,25,35,45,55,65,75,85]'*60;
        concat_WF = zeros(151*8,1);
        sudo_time_WF = linspace(5,85,151*8)+3;
        sudo_spike_center = linspace(5,75,8)+5;
        for k = 1:numel(sampling_windows)-1
            temp = (this_st >= sampling_windows(k) & this_st <= sampling_windows(k+1));
            if sum(temp) > 100
                temp_WF = this_WF(temp,:);
                temp_WF = datasample(temp_WF,100);
                temp_WF = temp_WF-mean(temp_WF(:,1:20),2);
                concat_WF((k-1)*151+1:(k-1)*151+121) = mean(temp_WF,1);
                concat_WF((k-1)*151+121:(k-1)*151+151) = NaN;
            else
                concat_WF((k-1)*151+1:k*151) = NaN;
            end
        end
        
        axis_temp = (1:numel(firing_rate_array))'/600;
        fig1 = figure("color","white",'Position', [50 50 650 750]);
        set(fig1, 'Visible', 'off');
        subplot(3,1,1);
        
        hold on
        % histogram(firing_rate_array)
        % set(gca, 'YScale', 'log')
        firing_rate_array_good = firing_rate_array;
        firing_rate_array_good(index_goodsig==0) = NaN;
        firing_rate_array_bad = firing_rate_array;
        firing_rate_array_bad(index_goodsig==1) = NaN;
        hold on
        plot(axis_temp, ...
            firing_rate_array_good,'r','LineWidth',1)
        plot(axis_temp, ...
            firing_rate_array_bad, 'b','LineWidth',1)
        set(gca, 'XGrid', 'on', 'YGrid', 'off')
        yline(threshold,'k--')
        xlim([0,90])
        xlabel('Session time(min)')
        title('Spiking activity filtered by moving aerage of 15s')
        subplot(3,1,2);
        plot(dshift_time,dshift)
        set(gca, 'XGrid', 'on', 'YGrid', 'off')
        xlim([0,90])
        xlabel('Session time(min)')
        ylim([-35,35])
        ylabel('um')
        title('Estimated Drifting (by KS) of 9 sampled channels')
        subplot(3,1,3);
        plot(sudo_time_WF,concat_WF,'LineWidth',1)
        hold on
        for k = sudo_spike_center
            xline(k,'r--')
        end
        xlim([0,90])
        title("Mean of (100)sampled spikes' waveform in each 10 min window")
        sgtitle(strcat(session_name,': n=',num2str(n), ...
            '; KSlabel = ',KSLabel(n,:), ...
            '; good-sig window = ', num2str(sum((index_goodsig==1))/600), ' min'), ...
            'Interpreter', 'none')
        %disp('next')
        save_path = fullfile(save_folder,strcat(session_name,'n',num2str(n)));
        %saveas(fig1,save_path)
        exportgraphics(fig1,strcat(save_path,'.jpg'),'Resolution',200)
    end
end
close(wbar)
%%
