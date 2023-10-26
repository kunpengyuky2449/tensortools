clear
%% Find the path of all PARAMS data and spiking data per session per animal
project_path = "D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\good_behavior_good_neuron_project";
load("D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\good_behavior_good_neuron_project\gbgn_summary_with_cluster_selection.mat")
%% hyperparameters
L_gaus = 800; %unit: ms, decide number of sample per filter to cover
a_gaus = 4; %unit: standard deviation, decide how wide the distribution spread, better keep this same while change L_gaus
%w_gaus = gausswin(L_gaus,a_gaus);
%wvtool(w_gaus);
%% Loop through session params to matched reward blocks, find qualified trials, reward/response time
% Qualified trial standard: 
% 1) tc after first recorded neuron spike; 
% 2) signal 3s before last recorded neuron spike
% 3) Not a catch trial (which has only tc but now signal tone)
% 4) tc length > 1000 ms; 
% 5) As requested, now low reward block 1 is always excluded
% sig:signal; tc: tone cloud; tbd: to be deleted
% First hit trial in each block always included regardless qualified standard (1)
%PARAMS.TrialResult: 0 Hit,1 Miss,2 FA,3 CR

%Parameters for saving as half product
wbar = waitbar(0,'Start Processing...');
for i = 1:numel(gbgn_summary.params_path)
    load(gbgn_summary.params_path{i})
    load(gbgn_summary.spikes_path{i})
    time_first_sp = sp_times(1)/PARAMS.Fs*1000;
    time_last_sp = sp_times(end)/PARAMS.Fs*1000-3000;
    ds_sp_time = round(sp_times/PARAMS.Fs*1000);
    spike_clusters=spike_clusters(1:numel(ds_sp_time));
    clear sp_times 
    %remove catch trials and false alarm trials since no signal played
    tbd1 = (PARAMS.CT == 1) | (PARAMS.LickRTwrtTC' < PARAMS.TCdur_sec);
    Target_Start = PARAMS.Target_Start(tbd1 == 0);
    TCdur = PARAMS.TCdur_sec(tbd1 == 0);
    Block_type = PARAMS.Block_type(tbd1 == 0);
    Trial_result = PARAMS.TrialResult(tbd1 == 0);
    Lick_fromTC = PARAMS.LickRTwrtTC(tbd1 == 0)';
    TC_Start = PARAMS.TcStart_LV(tbd1 == 0)';
    hit_trials = Lick_fromTC > TCdur;
    reaction_time =  Lick_fromTC - TCdur;
    % Find first evidednce of changing reward rate, 1: decrease to low reward, 
    % 2: increase to high reward, those trials are important and wont be
    % excluded
    reward_evidednce = zeros(size(Target_Start));
    reward_expectation = 1;
    for j = 1:numel(hit_trials)
        reward_sensed = Block_type(j);
        if  hit_trials(j) == 1 && reward_sensed ~= reward_expectation
            reward_expectation = reward_sensed;
        end
        reward_evidednce(j) = reward_expectation;
    end
    HR_next5 = double(hit_trials);
    flt_temp = cat(1,ones(5,1)*0.2,zeros(6,1));
    HR_next5 = filter(flt_temp, 1, HR_next5);
    % 
    ds_sig_start = round(Target_Start/PARAMS.Fs*1000);
    ds_tc_start = round(TC_Start /PARAMS.Fs*1000);
    Tensor_PerSession_Summary.walktimes = get_walktimes( ...
        PARAMS.WheelSpeed, ...
        PARAMS.WheelTimeSec);
    if height(Tensor_PerSession_Summary.walktimes) > 0
        walk_checks_temp = [];
        for j = 1:height(Tensor_PerSession_Summary.walktimes)
            walk_checks_temp = cat(1, walk_checks_temp, ...
                (Tensor_PerSession_Summary.walktimes(j,1):0.1: ...
                Tensor_PerSession_Summary.walktimes(j,2))');
        end
    end
    tbd2 = [];
    % Qualified trial standard: 
    temp = [-1;diff(reward_evidednce)];
    trials_window = gbgn_summary.good_behave_time{i};
    first_high_context_evidence = find(reward_evidednce == 2,1);
    for j = 1:numel(ds_sig_start)
        % 5) As requested, now low reward block 1 is always excluded
        if j < first_high_context_evidence
            tbd2(end+1) = j;%#ok
            continue
        end
        
        % 1) remove signal displayed 2s before target window start; 
        if ds_tc_start(j) < (trials_window(1)-2)*1000
            tbd2(end+1) = j; %#ok
            continue
        end
        % 2) remove signal displayed 2s after target window end; 
        if ds_sig_start(j) > (trials_window(2)+2)*1000
            tbd2(end+1) = j; %#ok
            continue
        end
        % First hit trial in each block always included regardless
        % qualifying standard (3) and (4)
        
        if temp(j) ~= 0
            continue
        end
        
        % 3) remove trials that mice has locomotion: any walks in sample
        % window -1s to +3s around signal display
        % To prevent dependency that locomotuion led by rewarding,
        % locomotion after mouse hit-lick wont exlcude the trial
        if reaction_time(j) > 0  
            move_check = sum(...
                walk_checks_temp > ds_sig_start(j)/1000-1  & ...
                walk_checks_temp < ds_sig_start(j)/1000+reaction_time(j), ...
                [1,2]);
        else
            move_check = sum(...
                walk_checks_temp > ds_sig_start(j)/1000-1  & ...
                walk_checks_temp < ds_sig_start(j)/1000+3, ...
                [1,2]);
        end

        if move_check > 0
            tbd2(end+1) = j;%#ok
            continue
        end

        % (removed standard)4) tc length < 1000 ms; 
        if TCdur(j)*1000 < 1000
            tbd2(end+1) = j; %#ok
            %disp(TCdur(j))
            continue
        end
    end

    ds_tc_start(tbd2) = [];
    ds_sig_start(tbd2) = [];
    reward_evidednce(tbd2) = [];
    
    reaction_time(tbd2) = [];
    hit_trials(tbd2) = [];
    HR_next5(tbd2) = [];
    % Each sensed-reward block includes trials when reward_expectation 
    % holds same until next reward_evidence different from reward_expectation
    % Initial expectation to be low reward, initial block is low reward
    % block.
    % Calculate how many trials are included in each block
    % save preprocessed data for this session 
    Tensor_PerSession_Summary.ds_sig_start = ds_sig_start; 
    Tensor_PerSession_Summary.ds_tc_start = ds_tc_start;
    Tensor_PerSession_Summary.reward_evidednce = reward_evidednce;
    Tensor_PerSession_Summary.reaction_time = reaction_time; 
    Tensor_PerSession_Summary.hit_trials = hit_trials; 
    Tensor_PerSession_Summary.HR_next5 = HR_next5;
    Tensor_PerSession_Summary.L_gaus = L_gaus;
    Tensor_PerSession_Summary.a_gaus = a_gaus;
    Tensor_PerSession_Summary.trials_PerBlock =  diff([1; ...
        find(diff(reward_evidednce)~=0)+1;numel(reward_evidednce)+1]);
    % Remove neurons that only fire at some portion of the session
    % 2023-10-23: For this subproject good neuron with good behavior window
    % has already been selected in previous steps, no need of further
    % filtering.
    % Remove decide good signal window based on heavily filtered firing rate
    % Remove those with good sig window length < 10 min
    % replace cpf difference to uniform distribution cdf method
    % replace pairwise eucidiean disdance method
    % Now decide the shape of the tensor
    T_timestamps = numel(2000-999:2000+3000);  % Trial Window around signal -1000ms to 3000ms
    K_trials = numel(ds_sig_start);
    N_neurons = numel(gbgn_summary.selected_cluster{i});
    Tensor_PerSession_Summary.neuron_depth = depths_with_clus_id(gbgn_summary.selected_cluster{i});
    
    Tensor_PerSession_Summary.N_neurons = N_neurons;
    Tensor_PerSession_Summary.T_timestamps = T_timestamps;
    Tensor_PerSession_Summary.K_trials = K_trials;
    % Create Tensor for this project, no need of mask tensor for now
    Session_tensor = zeros(Tensor_PerSession_Summary.N_neurons, ...
        Tensor_PerSession_Summary.T_timestamps, ...
        Tensor_PerSession_Summary.K_trials,'single');
    uniq_neurons = unique(spike_clusters); 
    uniq_neurons = uniq_neurons(gbgn_summary.selected_cluster{i});
    Tensor_PerSession_Summary.selected_uniq_neurons = uniq_neurons;
    for n = 1:numel(uniq_neurons)
        %update waitbar
        waitbar((i-1)/numel(gbgn_summary.params_path) + ...
                n/numel(uniq_neurons)/ numel(gbgn_summary.params_path), ...
                wbar,strcat('Processing session:',num2str(i)));
        % extract spiking data from the specific neuron cluster
        firing_rate_array = zeros(91*60*1000,1);
        ds_sp_time_iit = ds_sp_time(spike_clusters == uniq_neurons(n));
        [N_spike,time_bin] = groupcounts(ds_sp_time_iit);
        firing_rate_array(time_bin,1)=...
            firing_rate_array(time_bin,1) + N_spike;
        w_gaus = gausswin(Tensor_PerSession_Summary.L_gaus, ...
            Tensor_PerSession_Summary.a_gaus);
        firing_rate_array = filter(w_gaus, 1, firing_rate_array);
        for k = 1:numel(Tensor_PerSession_Summary.ds_sig_start)
            sig_time = Tensor_PerSession_Summary.ds_sig_start(k);
            tobeinsert = single(firing_rate_array(sig_time-999:sig_time+3000));
            Session_tensor(n,:,k) = tobeinsert;
        end
    end
    Tensor_PerSession_Summary.file_path_Spikes = gbgn_summary.spikes_path{i};
    Tensor_PerSession_Summary.file_path_PARAMS = gbgn_summary.params_path{i};
    % Save tenor and tensor summary
    file_name = gbgn_summary.params_path{i}( ...
                strfind(gbgn_summary.params_path{i}, ...
                'preprocess\')+11:end);
    file_name(strfind(file_name,'\'))='_';
    file_name = file_name(1:strfind(file_name,'PARAMS')-2);
    Tensor_PerSession_Summary.session_name = file_name;
    name_temp = strcat(file_name,'_signal_tensor_GBGN_summary_', datestr(now,'yymmdd_HH_MM_SS')); %#ok
    save(fullfile(project_path, name_temp),"Tensor_PerSession_Summary")
    name_temp = strcat(file_name,'_signal_tensor_GBGN_NTK_', datestr(now,'yymmdd_HH_MM_SS')); %#ok
    save(fullfile(project_path, name_temp),"Session_tensor",'-v7.3')
    clear Session_tensor Tensor_PerSession_Summary spike_clusters ds_sp_time tobeinsert
end
close all
close(wbar)
clear
%%
% test_fr = firing_rate_matrix(:,1);
% %hyperparameters
% L_gaus = 400; %unit: ms, decide number of sample per filter to cover
% a_gaus = 4; %unit: standard deviation, decide how wide the distribution spread
% w_gaus = gausswin(L_gaus,a_gaus);
% %wvtool(w_gaus);
% y = filter(w_gaus, 1, test_fr);
% plot(y)
%%


