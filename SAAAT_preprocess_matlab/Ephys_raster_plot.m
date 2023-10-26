clear
%%
load('SAAAT_session_summary_073123.mat')
wbar = waitbar(0,'Start plotting...');
for i = 8:numel(session_summary.file_path_PARAMS) %!
    load(session_summary.file_path_PARAMS{i})
    load(session_summary.file_path_Spikes{i})
    tbd1 = (PARAMS.CT == 1) | (PARAMS.LickRTwrtTC' < PARAMS.TCdur_sec);
    Target_Start = PARAMS.Target_Start(tbd1 == 0);
    TC_Start = PARAMS.TcStart_LV(tbd1 == 0)';
    ds_sig_start = round(Target_Start/PARAMS.Fs*1000);
    ds_tc_start = round(TC_Start /PARAMS.Fs*1000);
    ds_sp_time = round(sp_times/PARAMS.Fs*1000);
    clear Target_Start TC_Start sp_times
    spike_clusters=spike_clusters(1:numel(ds_sp_time));
    uniq_neurons = unique(spike_clusters);
    subplot_row = 0;
    plot_id = 0; 
    for n = 1:numel(uniq_neurons)
        % extract spiking data from the specific neuron cluster
        % Down sample spikes sampling rate to ms
        ds_sp_time_iit = ds_sp_time(spike_clusters == uniq_neurons(n));
        % Map spiking time to a binary time series
        if subplot_row == 0
            fig = figure('color','white','Position',[50,50,1250,850]); 
            set(fig, 'Visible', 'off');
            for p = 1:8
                subplot(4,2,p)
                hold on
                xline(0,'r--')
                xlim([-2000,+2000])
            end
            subplot(4,2,1)
            title('Pink Noise To Tone Cloud')
            subplot(4,2,2)
            title('Tone Cloud To Signal')
            title_new = session_summary.file_path_Spikes{i}( ...
                strfind(session_summary.file_path_Spikes{i}, ...
                'preprocess\')+11:end);
            sgtitle(title_new, 'Interpreter', 'none')
            waitbar((i-1)/numel(session_summary.file_path_PARAMS) + ...
                n/numel(uniq_neurons)/ numel(session_summary.file_path_PARAMS), ...
                wbar,strcat('Processing session:',num2str(i)));
        end
        trial_temp_tc = [];
        trial_temp_sig  = [];
        spikes_temp_tc = [];
        spikes_temp_sig = [];
        for t = 1:numel(ds_sig_start)
            temp_1 = ds_sp_time_iit(ds_sp_time_iit< ds_tc_start(t)+2000 & ...
                ds_sp_time_iit > ds_tc_start(t)-2000) - ds_tc_start(t);
            temp_2 = ds_sp_time_iit(ds_sp_time_iit< ds_sig_start(t)+2000 & ...
                ds_sp_time_iit > ds_sig_start(t)-2000) - ds_sig_start(t);
            spikes_temp_tc = cat(1, spikes_temp_tc, temp_1);
            spikes_temp_sig = cat(1, spikes_temp_sig, temp_2);
            trial_temp_tc = cat(1, trial_temp_tc,t*ones(size(temp_1)));
            trial_temp_sig = cat(1, trial_temp_sig, t*ones(size(temp_2)));
        end
        subplot(4,2,subplot_row*2+1)
        scatter(spikes_temp_tc,trial_temp_tc,15,'k|')
        ylabel(strcat('Neuron',num2str(n),'-Trial:'))
        subplot(4,2,subplot_row*2+2)
        scatter(spikes_temp_sig,trial_temp_sig,15,'k|')
        subplot_row = subplot_row+1;
        
        if subplot_row == 4 || n == numel(uniq_neurons)
            title_new(strfind(title_new,'\'))='_'; %#ok
            title_new = strcat(title_new(1:end-4),num2str(plot_id));
            plot_id = plot_id + 1;
            file_title = fullfile('D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\Neuron_raster_plot',...
                title_new);
            %saveas(fig,file_title)
            saveas(fig,file_title,'jpeg')
            subplot_row = 0;
        end
    end
    % Rlease ram
    clear PARAMS spike_clusters ds_sp_time_iit sp_time_iit
end
close(wbar)