clear
load('D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\good_behavior_good_neuron_project\gbgn_summary.mat','gbgn_summary')
%%
target_folder = "D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\good_behavior_good_neuron_project";
selected_folder = "D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\good_behavior_good_neuron_project\Selected";
excluded_folder = "D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\good_behavior_good_neuron_project\Excluded";
gbgn_summary.selected_cluster = cell(0);
for i = 1:numel(gbgn_summary.session_list)
    load(gbgn_summary.spikes_path{i})
    load(gbgn_summary.params_path{i})
    load(gbgn_summary.waveform_path{i})
    Fs = PARAMS.Fs;
    ds_sp_time = round(sp_times/Fs*10);
    clear Target_Start TC_Start sp_times
    cd(gbgn_summary.RawEphys_folder{i})
    raw_clusterID = readNPY('spike_clusters.npy');
    raw_st = readNPY('spike_times.npy')/Fs;
    load("rez.mat","rez");
    batch_samples = rez.ops.NT;
    dshift = rez.dshift;
    dshift_time = (1:height(dshift))'* batch_samples/Fs/60;
    spike_clusters=spike_clusters(1:numel(ds_sp_time));
    uniq_neurons = unique(spike_clusters);
    clear rez PARAMS
    selected_cluster = [];
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
        % histogram(firing_rate_array)
        % set(gca, 'YScale', 'log')
        firing_rate_array_good = firing_rate_array;
        firing_rate_array_good(index_goodsig==0) = NaN;
        firing_rate_array_bad = firing_rate_array;
        firing_rate_array_bad(index_goodsig==1) = NaN;
        hold on
        window_start = gbgn_summary.good_behave_time{i}(1)/60;
        window_end = gbgn_summary.good_behave_time{i}(2)/60;
        plot(axis_temp, ...
            firing_rate_array_good,'r','LineWidth',1)
        plot(axis_temp, ...
            firing_rate_array_bad, 'b','LineWidth',1)
        patch([window_start, window_end, window_end,window_start], ...
            [max(ylim),max(ylim),0,0],'k','EdgeColor','none',...
            'FaceAlpha',.1)
        set(gca, 'XGrid', 'on', 'YGrid', 'off')
        yline(threshold,'k--')
        xlim([0,90])
        xlabel('Session time(min)')
        title('Spiking activity filtered by moving aerage of 15s')
        subplot(3,1,2);
        plot(dshift_time,dshift)
        patch([window_start, window_end, window_end,window_start], ...
            [max(ylim),max(ylim),0,0],'k','EdgeColor','none',...
            'FaceAlpha',.1)
        set(gca, 'XGrid', 'on', 'YGrid', 'off')
        xlim([0,90])
        xlabel('Session time(min)')
        ylim([-35,35])
        ylabel('um')
        title('Estimated Drifting (by KS) of 9 sampled channels')
        subplot(3,1,3);
        plot(sudo_time_WF,concat_WF,'LineWidth',1)
        patch([window_start, window_end, window_end,window_start], ...
            [max(ylim),max(ylim),0,0],'k','EdgeColor','none',...
            'FaceAlpha',.1)
        hold on
        for k = sudo_spike_center
            xline(k,'r--')
        end
        xlim([0,90])
        title("Mean of (100)sampled spikes' waveform in each 10 min window")
        sgtitle(strcat(gbgn_summary.session_list{i},': n=',num2str(n), ...
            '; KSlabel = ',KSLabel(n,:), ...
            '; good-sig window = ', num2str(sum((index_goodsig==1))/600), ' min'), ...
            'Interpreter', 'none')
        window_range = round(gbgn_summary.good_behave_time{i}*10);
        good_sig_ratio = sum(index_goodsig(window_range(1):window_range(2)))/...
            (window_range(2)-window_range(1));
        if good_sig_ratio < 0.7
            save_path = fullfile(excluded_folder, ...
                strcat(gbgn_summary.session_list{i},'n',num2str(n)));
            exportgraphics(fig1,strcat(save_path,'.jpg'),'Resolution',200)
            disp(strcat("Cluster ",num2str(n)," is inactive in target window, thus excluded"))
        else
            set(fig1, 'Visible', 'on');
            prompt = "Select this cluster? 1/0 [Y]: ";
            txt = input(prompt,"s");
            if isempty(txt)
                txt = "1";
            end
            if txt == "1"
                selected_cluster(end+1) = n; %#ok
                save_path = fullfile(selected_folder, ...
                    strcat(gbgn_summary.session_list{i},'n',num2str(n)));
                exportgraphics(fig1,strcat(save_path,'.jpg'),'Resolution',200)
                disp(strcat("Cluster ",num2str(n)," is SELECTED"))
            else
                save_path = fullfile(excluded_folder, ...
                    strcat(gbgn_summary.session_list{i},'n',num2str(n)));
                exportgraphics(fig1,strcat(save_path,'.jpg'),'Resolution',200)
                disp(strcat("Cluster ",num2str(n)," is excluded"))
            end
        end
    end
    gbgn_summary.selected_cluster{i} = selected_cluster;
end
%%
save_name = fullfile(target_folder,'gbgn_summary_with_cluster_selection');
save(save_name,"gbgn_summary")