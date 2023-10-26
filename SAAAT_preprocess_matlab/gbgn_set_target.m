clear
%%
target_folder = "D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\good_behavior_good_neuron_project";
gbgn_summary = struct();
gbgn_summary.session_list = {"W3414_", "W3337_6", "W3337_5", "W3337_4", "W3337_3", ...
    "W3337_2", "W3336_4", "W3333_1"};
range_temp ={[1,118],[3,175],[1,150],[1,175],[1,175],[3,260],[3,210],[1,161]};
%% Find the path of all PARAMS data and spiking data per session per animal
%load('C:\Users\78184\OneDrive\Desktop\TCA_project\SAAAT_preprocess_matlab\test_data_step12\W3333_1_Spikes.mat')
%load('C:\Users\78184\OneDrive\Desktop\TCA_project\SAAAT_preprocess_matlab\test_data_step12\PARAMS.mat')
listing = dir('D:\KunpengY\AE_trained_behaving\preprocess');
load('SAAAT_session_summary_073123.mat')
file_path = 'D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\tensor_perSession_v3';
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
clear listing file_listing sp_file_name sub_listing
%%
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
clear listing file_listing sub_listing destination
%%
gbgn_summary.params_path = cell(0);
gbgn_summary.spikes_path = cell(0);
gbgn_summary.waveform_path = cell(0);
gbgn_summary.RawEphys_folder = cell(0);
for i = 1:numel(gbgn_summary.session_list)
    temp = strrep(gbgn_summary.session_list{i},'_','\');
    gbgn_summary.params_path{end+1}  = file_path_PARAMS{contains(file_path_PARAMS, ...
        temp)};
    gbgn_summary.spikes_path{end+1} = file_path_Spikes{contains(file_path_Spikes, ...
        gbgn_summary.session_list{i})};
    gbgn_summary.RawEphys_folder{end+1} = folder_path_EphysRaw{contains(folder_path_EphysRaw, ...
        temp)};
    temp2 = contains({listing_waveform.name}, gbgn_summary.session_list{i});
    gbgn_summary.waveform_path{end+1} = fullfile(listing_waveform(temp2).folder, ...
        listing_waveform(temp2).name);
end
clear temp temp2 subsub_listing listing_waveform
%%
gbgn_summary.trial_NoCT_range = cell(0);
gbgn_summary.good_behave_time = cell(0);
for i = 1:numel(range_temp)
    load(gbgn_summary.params_path{i})
    temp = contains(session_summary.file_path_PARAMS, ...
        gbgn_summary.params_path{i});
    end_time = session_summary.ds_sig_start_all{temp}(range_temp{i}(2))/1000;
    tbd1 = (PARAMS.CT == 1) | (PARAMS.LickRTwrtTC' < PARAMS.TCdur_sec);
    Target_Start = PARAMS.Target_Start(tbd1 == 0);
    TCdur = PARAMS.TCdur_sec(tbd1 == 0);
    Block_type = PARAMS.Block_type(tbd1 == 0);
    Lick_fromTC = PARAMS.LickRTwrtTC(tbd1 == 0)';
    hit_trials = Lick_fromTC > TCdur;
    reaction_time =  Lick_fromTC - TCdur;
    % Find first evidednce of changing reward rate, 1: decrease to low reward, 
    % 2: increase to high reward, those trials are important and wont be
    % excluded
    reward_expectation = 1;
    new_block_index = [];
    for j = 1:numel(hit_trials)
        reward_sensed = Block_type(j);
        if  hit_trials(j) == 1 && reward_sensed ~= reward_expectation
            reward_expectation = reward_sensed;
            new_block_index(end+1) = j; %#ok
        end
    end
    gbgn_summary.trial_NoCT_range{i} = [new_block_index(range_temp{i}(1)), ...
        range_temp{i}(2)];
    gbgn_summary.good_behave_time{i} = [...
        Target_Start(new_block_index(range_temp{i}(1)))/PARAMS.Fs,...
        end_time];
end
%%
save_name = fullfile(target_folder,'gbgn_summary');
save(save_name,"gbgn_summary")