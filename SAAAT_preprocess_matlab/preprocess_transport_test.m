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
% Qualified trial standard: 1) tc length > 500ms + 1/2 gaussian filter window; 
% 2) tc after first recorded neuron spike; 
% 3) signal 3s before last recorded neuron spike
% 4) Not a catch trial (which has only tc but now signal tone)
% sig:signal; tc: tone cloud; tbd: to be deleted
% First hit trial in each block always included regardless qualified standard (1)
%PARAMS.TrialResult: 0 Hit,1 Miss,2 FA,3 CR

%Parameters for saving as half product
reward_block_qualified_matrix = zeros(12,numel(file_path_PARAMS));
ds_sig_start_all = cell(0);
reward_evidednce_all = cell(0);
reaction_time_all = cell(0);
hit_trials_all = cell(0);
Block_type_all = cell(0);

for i = 1:numel(file_path_PARAMS)
    load(file_path_PARAMS{i})
    load(file_path_Spikes{i})
    time_first_sp = sp_times(1)/PARAMS.Fs*1000;
    time_last_sp = sp_times(end)/PARAMS.Fs*1000-3000;
    clear sp_times spike_clusters depths_with_clus_id
    %remove catch trials and false alarm trials since no signal played
    tbd1 = (PARAMS.CT == 1) | (PARAMS.LickRTwrtTC' < PARAMS.TCdur_sec);
    Target_Start = PARAMS.Target_Start(tbd1 == 0);
    TCdur = PARAMS.TCdur_sec(tbd1 == 0);
    Block_type = PARAMS.Block_type(tbd1 == 0);
    Trial_result = PARAMS.TrialResult(tbd1 == 0);
    Lick_fromTC = PARAMS.LickRTwrtTC(tbd1 == 0)';
    TC_Start = PARAMS.TcStart_LV(tbd1 == 0)';
    hit_trials = Lick_fromTC > TCdur;
    % Find first evidednce of changing reward rate, 1: decrease to low reward, 
    % 2: increase to high reward, those trials are important and wont be
    % excluded
    reward_evidednce = zeros(size(Target_Start));
    reward_expectation = 1;
    for j = 1:numel(hit_trials)
        reward_sensed = Block_type(j);
        if  hit_trials(j) == 1 && reward_sensed ~= reward_expectation
            reward_evidednce(j) = reward_sensed;
            reward_expectation = reward_sensed;
        end
    end
    % 
    ds_sig_start = round(Target_Start/PARAMS.Fs*1000);
    ds_tc_start = round(TC_Start /PARAMS.Fs*1000);
    tbd2 = [];
    % Qualified trial standard: 

    for j = 1:numel(ds_sig_start)
        % 2) tc after first recorded neuron spike; 
        if ds_tc_start(j) < time_first_sp
            tbd2(end+1) = j; %#ok
            continue
        end
        % 3) signal 3s before last recorded neuron spike
        if ds_sig_start(j) > time_last_sp
            tbd2(end+1) = j; %#ok
            continue
        end
        % First hit trial in each block always included regardless qualified standard (1)
        if reward_evidednce(j) > 0
            continue
        end
        % 1) tc length > 500ms + 1/2 gaussian filter window; 
        if TCdur(j)*1000 < 500+L_gaus/2
            tbd2(end+1) = j; %#ok
            %disp(TCdur(j))
            continue
        end
    end
    ds_sig_start(tbd2) = [];
    reward_evidednce(tbd2) = [];
    reaction_time =  Lick_fromTC - TCdur;
    reaction_time(tbd2) = [];
    hit_trials(tbd2) = [];
    Block_type(tbd2) = [];
    % Each sensed-reward block includes trials when reward_expectation 
    % holds same until next reward_evidence different from reward_expectation
    % Initial expectation to be low reward, initial block is low reward
    % block.
    % Calculate how many trials are included in each block
    reward_block_qualified_session = diff([1;find(reward_evidednce>0);numel(reward_evidednce)])+1;
    % save preprocessed data for this session 
    reward_block_qualified_matrix(1:numel(reward_block_qualified_session),i) = ...
        reward_block_qualified_matrix(1:numel(reward_block_qualified_session),i) + ...
        reward_block_qualified_session;
    ds_sig_start_all{end+1} = ds_sig_start; %#ok
    reward_evidednce_all{end+1} = reward_evidednce; %#ok
    reaction_time_all{end+1} = reaction_time; %#ok
    hit_trials_all{end+1} = hit_trials; %#ok
    Block_type_all{end+1} = Block_type; %#ok
    % Set a simple progress indicator
    disp(strcat(num2str(i/numel(file_path_PARAMS)*100),'%'))
end
%% save temporary data
name_temp = strcat('SAAAT_session_summary_', datestr(now,'mmddyy')); %#ok
session_summary.reward_block_qualified_matrix = reward_block_qualified_matrix;
session_summary.ds_sig_start_all = ds_sig_start_all;
session_summary.reward_evidednce_all = reward_evidednce_all;
session_summary.reaction_time_all = reaction_time_all;
session_summary.hit_trials_all = hit_trials_all;
session_summary.Block_type_all = Block_type_all;
session_summary.L_gaus = L_gaus;
session_summary.a_gaus = a_gaus;
session_summary.file_path_PARAMS = file_path_PARAMS;
session_summary.file_path_Spikes = file_path_Spikes;
save(name_temp,"session_summary")

%% Rerun previous part if L_gaus changed, otherwise start from here
clear
load("SAAAT_session_summary_073123")
L_gaus = session_summary.L_gaus; %unit: ms, decide number of sample per filter to cover
a_gaus = session_summary.a_gaus; %unit: standard deviation, decide how wide the distribution spread, better keep this same while change L_gaus
%w_gaus = gausswin(L_gaus,a_gaus);
%% Decide how many session to be included and how many trial per block to be included
reward_block_matrix = session_summary.reward_block_qualified_matrix;
session_index = 1:width(reward_block_matrix);
% Sicne all session has to be alligned by reward-amount block for allow same num of trials in each block
% Sessions having less than 8 block or having trial less than 15 in any block would be discarded
tbd3 = [];
for i = session_index
    if min(reward_block_matrix(1:8,i)) < 15
        tbd3(end+1) = i; %#ok
    end
end
% Sessions to be selected for analysis
session_index(tbd3) = [];
reward_block_matrix(:,tbd3) = [];
reward_block_matrix(9:end,:) = [];
% Decide trials per block to be selected for analysis, which is number of trials
% for the session having minimum trials in that block
trials_per_block = min(reward_block_matrix,[],2);
% Now decide the shape of the tensor
T_timestamps = numel(1000-499:1000+3000);  % Trial Window around signal -500ms to 3000ms
K_trials = sum(trials_per_block);
N_neurons = 0;
neuron_depth_per_session = cell(0);
neuron_number_per_session = cell(0);
for i = session_index
    load(session_summary.file_path_Spikes{i})
    N_neurons = N_neurons + numel(unique(spike_clusters));
    % Save neuron information here to save ram for later
    neuron_depth_per_session{end+1} = depths_with_clus_id(:,1);%#ok<SAGROW>
    neuron_number_per_session{end+1} = numel(unique(spike_clusters));%#ok<SAGROW>
end
tensor_summaray.a_gaus = a_gaus;
tensor_summaray.L_gaus = L_gaus;
tensor_summaray.N_neurons = N_neurons;
tensor_summaray.T_timestamps = T_timestamps;
tensor_summaray.K_trials = K_trials;
tensor_summaray.session_index = session_index;
tensor_summaray.trials_per_block = trials_per_block;
tensor_summaray.neuron_number_per_session = neuron_number_per_session;
tensor_summaray.neuron_depth_per_session = neuron_depth_per_session;
name_temp = strcat('SAAAT_tensor_summary_', datestr(now,'mmddyy')); %#ok
% Save halp product here
save(name_temp,"tensor_summaray")
clear depths_with_clus_id sp_times spike_clusters neuron_depth_per_session tensor_summaray
%% Create the tensor
clear
load('SAAAT_tensor_summary_080123')
load("SAAAT_session_summary_073123")
load(session_summary.file_path_PARAMS{1})
Fs = PARAMS.Fs;
clear PARAMS
%%
session_index = tensor_summaray.session_index;
L_gaus = session_summary.L_gaus; %unit: ms, decide number of sample per filter to cover
a_gaus = session_summary.a_gaus; %unit: standard deviation, decide how wide the distribution spread, better keep this same while change L_gaus
w_gaus = gausswin(L_gaus,a_gaus);
%%
SAAAT_tensor = zeros(tensor_summaray.N_neurons, ...
    tensor_summaray.T_timestamps, ...
    tensor_summaray.K_trials,'single');
N_insert = 1;
K_insert = 1;
k_track = [];
n_track = [];
for i = session_index
    load(session_summary.file_path_Spikes{i})
    % Find trial indice selected for analysis for this session
    new_block = [1;find(session_summary.reward_evidednce_all{i} > 0)];
    trials_selected = [];
    for j = 1:numel(tensor_summaray.trials_per_block)
        trials_selected = cat(2,trials_selected, ...
            new_block(j):(new_block(j)+tensor_summaray.trials_per_block(j)-1));
    end
    % get signal time for those trials as coordinates of spiking data
    ds_sig_selected = session_summary.ds_sig_start_all{i}(trials_selected);% Unit ms
    spike_clusters=spike_clusters(1:numel(sp_times));
    uniq_neurons = unique(spike_clusters);
    for n = 1:numel(uniq_neurons)
        % extract spiking data from the specific neuron cluster
        firing_rate_array = zeros(91*60*1000,1);
        sp_time_iit = sp_times(spike_clusters == uniq_neurons(n));
        % Down sample spikes sampling rate to ms
        ds_sp_time_iit = round(sp_time_iit/Fs*1000);
        % Map spiking time to a binary time series
        [N_spike,time_bin] = groupcounts(ds_sp_time_iit);
        firing_rate_array(time_bin,1)=...
            firing_rate_array(time_bin,1) + N_spike;
        % Apply a gaussian filter to  convert time series of spikes
        % to pseudo-firing rate time series
        firing_rate_array = filter(w_gaus, 1, firing_rate_array);
        % Cut the trial data we need based on signal coordinates
        for k = ds_sig_selected'
            tobeinsert = single(firing_rate_array(k-499:k+3000));
            % Insert trial data to tensor
            SAAAT_tensor(N_insert,:,K_insert) = tobeinsert;
            %k_track(end+1) = K_insert;
            %n_track(end+1) = N_insert;
            K_insert = K_insert + 1;
        end
        K_insert = 1;
        N_insert = N_insert+1;
        % Set a simple progress indicator
        disp(strcat(num2str(N_insert/tensor_summaray.N_neurons*100),'%'))
    end 
end
%%
name_temp = strcat('SAAAT_tensor_NTK_', datestr(now,'yymmdd_HH_MM_SS')); %#ok
save(name_temp,"SAAAT_tensor",'-v7.3')
%%
clear
load('SAAAT_tensor_NTK_230801_13_34_06')
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


