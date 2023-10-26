clear
%% Find the path of all PARAMS data and spiking data per session per animal
%load('C:\Users\78184\OneDrive\Desktop\TCA_project\SAAAT_preprocess_matlab\test_data_step12\W3333_1_Spikes.mat')
%load('C:\Users\78184\OneDrive\Desktop\TCA_project\SAAAT_preprocess_matlab\test_data_step12\PARAMS.mat')
listing = dir('D:\KunpengY\AE_trained_behaving\preprocess');
listing = listing(~ismember({listing.name},{'.','..'}));
file_path_PARAMS = cell(0);
file_path_Spikes = cell(0);
for i = 1:numel(listing)
    sub_listing = dir(fullfile(listing(i).folder,listing(i).name));
    sub_listing = sub_listing(~ismember({sub_listing.name},{'.','..'}));
    for j = 1:numel(sub_listing)
        file_path_PARAMS{end+1} = fullfile(sub_listing(j).folder,sub_listing(j).name,'PARAMS');%#ok<SAGROW> 
        %find the file name of spiking data
        file_listing = dir(fullfile(sub_listing(j).folder,sub_listing(j).name));
        sp_file_name = file_listing(contains({file_listing.name},"Spike")).name;
        file_path_Spikes{end+1} = fullfile(sub_listing(j).folder,sub_listing(j).name,sp_file_name);%#ok<SAGROW> 
    end
end
%% hyperparameters
L_gaus = 400; %unit: ms, decide number of sample per filter to cover
a_gaus = 4; %unit: standard deviation, decide how wide the distribution spread, better keep this same while change L_gaus
%w_gaus = gausswin(L_gaus,a_gaus);
%wvtool(w_gaus);
%% Loop through session params to matched reward blocks, find qualified trials, reward/response time
% Qualified trial standard: 1) tc length > 1000 ms; 
% 2) tc after first recorded neuron spike; 
% 3) signal 3s before last recorded neuron spike
% 4) Not a catch trial (which has only tc but now signal tone)
% sig:signal; tc: tone cloud; tbd: to be deleted
% First hit trial in each block always included regardless qualified standard (1)
%PARAMS.TrialResult: 0 Hit,1 Miss,2 FA,3 CR

%Parameters for saving as half product
wbar = waitbar(0,'Start Processing...');
for i = 2:numel(file_path_PARAMS)
    load(file_path_PARAMS{i})
    load(file_path_Spikes{i})
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
    for j = 1:numel(ds_sig_start)
        % 1) tc after first recorded neuron spike; 
        if ds_tc_start(j) < time_first_sp
            tbd2(end+1) = j; %#ok
            continue
        end
        % 2) signal 3s before last recorded neuron spike
        if ds_sig_start(j) > time_last_sp
            tbd2(end+1) = j; %#ok
            continue
        end
        % First hit trial in each block always included regardless
        % qualified standard (3) and (4)
        
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
        end

        % (removed standard)4) tc length > 500ms + 1/2 gaussian filter window; 
        % if TCdur(j)*1000 < 1200
        %     tbd2(end+1) = j; %#ok
        %     %disp(TCdur(j))
        %     continue
        % end
    end
    ds_tc_start(tbd2) = [];
    ds_sig_start(tbd2) = [];
    reward_evidednce(tbd2) = [];
    
    reaction_time(tbd2) = [];
    hit_trials(tbd2) = [];
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
    Tensor_PerSession_Summary.L_gaus = L_gaus;
    Tensor_PerSession_Summary.a_gaus = a_gaus;
    Tensor_PerSession_Summary.trials_PerBlock =  diff([1; ...
        find(diff(reward_evidednce)~=0)+1;numel(reward_evidednce)+1]);
    % Remove neurons that only fire at some portion of the session
    % Now use cpf difference to uniform distribution cdf method
    % replace pairwise eucidiean disdance method
    uniq_neurons = unique(spike_clusters);
    % Pair_EDist = [];
    % for n = 1:numel(uniq_neurons)
    %     % extract spiking data from the specific neuron cluster
    %     ds_sp_time_iit = ds_sp_time(spike_clusters == uniq_neurons(n));
    %     ds_rate = floor(numel(ds_sp_time_iit)/3000);
    %     if  ds_rate > 1
    %         ds_sp_time_iit = downsample(ds_sp_time_iit, ds_rate);
    %     end
    %     Pair_EDist(end+1) = mean(pdist(ds_sp_time_iit)); %#ok
    % end
    cdf_diff = [];
    cdf_uniform = linspace(0,1,3000)';
    time_checks = linspace(time_first_sp,time_last_sp,3000);
    for n = 1:numel(uniq_neurons)
        % extract spiking data from the specific neuron cluster
        ds_sp_time_iit = ds_sp_time(spike_clusters == uniq_neurons(n));
        num_sp = numel(ds_sp_time_iit);
        cdf_empr = zeros(numel(time_checks),1);
        for j = 1: numel(time_checks)
            cdf_empr(j) = sum(ds_sp_time_iit <= time_checks(j))/num_sp;
        end
        cdf_diff(end+1) = sum(abs(cdf_uniform-cdf_empr))/3000; %#ok
    end
    % Exclude the 30% neurons that has least paired euclidean distance (least disperse) 
    uniq_neurons = uniq_neurons(cdf_diff < 0.2);
    % Now decide the shape of the tensor
    T_timestamps = numel(2000-999:2000+3000);  % Trial Window around signal -1000ms to 3000ms
    K_trials = numel(ds_sig_start);
    N_neurons = numel(uniq_neurons);
    Tensor_PerSession_Summary.neuron_depth = depths_with_clus_id( ...
        cdf_diff < 0.2,1);
    Tensor_PerSession_Summary.neuron_cdf_uniformity = cdf_diff(cdf_diff<0.2)';
    Tensor_PerSession_Summary.N_neurons = N_neurons;
    Tensor_PerSession_Summary.T_timestamps = T_timestamps;
    Tensor_PerSession_Summary.K_trials = K_trials;

    % Create Tensor
    Session_tensor = zeros(Tensor_PerSession_Summary.N_neurons, ...
        Tensor_PerSession_Summary.T_timestamps, ...
        Tensor_PerSession_Summary.K_trials,'single');

    for n = 1:numel(uniq_neurons)
        %update waitbar
        waitbar((i-1)/numel(file_path_PARAMS) + ...
                n/numel(uniq_neurons)/ numel(file_path_PARAMS), ...
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
            rwd_time = Tensor_PerSession_Summary.ds_sig_start(k);
            tobeinsert = single(firing_rate_array(rwd_time-999:rwd_time+3000));
            % Insert trial data to tensor
            Session_tensor(n,:,k) = tobeinsert;
        end
    end
    Tensor_PerSession_Summary.file_path_Spikes = file_path_PARAMS{i};
    Tensor_PerSession_Summary.file_path_PARAMS = file_path_Spikes{i};
    % Save tenor and tensor summary
    file_name = file_path_Spikes{i}( ...
                strfind(file_path_Spikes{i}, ...
                'preprocess\')+11:end);
    file_name(strfind(file_name,'\'))='_';
    file_name = file_name(1:strfind(file_name,'.mat')-1);
    Tensor_PerSession_Summary.session_name = file_name;
    file_path = 'D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\tensor_perSession_v2';
    name_temp = strcat(file_name,'_signal_tensor_V2_summary_', datestr(now,'yymmdd_HH_MM_SS')); %#ok
    save(fullfile(file_path, name_temp),"Tensor_PerSession_Summary")
    name_temp = strcat(file_name,'_signal_tensor_V2_NTK_', datestr(now,'yymmdd_HH_MM_SS')); %#ok
    save(fullfile(file_path, name_temp),"Session_tensor",'-v7.3')
    clear Session_tensor Tensor_PerSession_Summary spike_clusters ds_sp_time tobeinsert
end
close(wbar)
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


