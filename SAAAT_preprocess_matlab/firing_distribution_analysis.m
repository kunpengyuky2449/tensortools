%%
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
%% load each spike data
cdf_diff_sessions = {1,numel(file_path_PARAMS)};
cdf_uniform = linspace(0,1,3000)';
for i = 1:numel(file_path_PARAMS)
    disp(i)
    load(file_path_PARAMS{i})
    load(file_path_Spikes{i})
    time_first_sp = sp_times(1)/PARAMS.Fs*1000;
    time_last_sp = sp_times(end)/PARAMS.Fs*1000;
    ds_sp_time = round(sp_times/PARAMS.Fs*1000);
    spike_clusters=spike_clusters(1:numel(ds_sp_time));
    clear sp_times
    uniq_neurons = unique(spike_clusters);
    cdf_diff = [];
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
    cdf_diff_sessions{i} = cdf_diff;
end
%%
%set(0,'DefaultFigureWindowStyle','docked')
cdf_diff_all_session = [];
save_folder = fullfile( ...
        'D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\spike_distribution_analysis');
if ~exist(save_folder, 'dir')
   mkdir(save_folder)
end
for i = 1:numel(file_path_PARAMS)
    fig = figure('color','white','Position',[50,50,1250,850]); 
    set(fig, 'Visible', 'off');
    histogram(cdf_diff_sessions{1, i},50)
    cdf_diff_all_session = cat(1,cdf_diff_all_session,cdf_diff_sessions{1, i}');
    xlim([0,0.5])
    temp = strfind(file_path_Spikes{1},'\');
    session_name = file_path_Spikes{i}(temp(end)+1:end-4);
    title(strcat(session_name, ...
        ': histgram of CDF diffference from uniform distribution'), 'Interpreter', 'none')
    xlabel('CDF diffference from uniform distribution')
    ylabel('Number of neuron')
    file_name = fullfile(save_folder,...
                strcat(session_name,'_hist_cdf_difference'));
    %saveas(fig,file_title)
    saveas(fig,file_name,'jpeg')
    saveas(fig,file_name,'fig')
end
fig = figure('color','white','Position',[50,50,1250,850]); 
set(fig, 'Visible', 'off');
histogram(cdf_diff_all_session,50)
xlabel('CDF diffference from uniform distribution')
ylabel('Number of neuron')
title('All 17 sessions: histgram of CDF diffference from uniform distribution')
file_name = fullfile(save_folder,...
                strcat('all_17_sessions_hist_cdf_difference'));
%saveas(fig,file_title)
saveas(fig,file_name,'jpeg')
saveas(fig,file_name,'fig')
%set(0,'DefaultFigureWindowStyle','normal')


