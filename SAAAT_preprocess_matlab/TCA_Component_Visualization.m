function fig_output = TCA_Component_Visualization(Session,rank,component,anomaly)
%
% Session = 'W3333_4';
% rank = 7;
% component = 6;
%% Find the path of all pupil data per session per animal
listing = dir('D:\KunpengY\AE_trained_behaving\preprocess');
listing = listing(~ismember({listing.name},{'.','..'}));
file_path_pupil = cell(0);
for i = 1:numel(listing)
    sub_listing = dir(fullfile(listing(i).folder,listing(i).name));
    sub_listing = sub_listing(~ismember({sub_listing.name},{'.','..'}));
    for j = 1:numel(sub_listing)
        %find the file name of spiking data
        file_listing = dir(fullfile(sub_listing(j).folder,sub_listing(j).name));
        sp_file_name = file_listing(contains({file_listing.name},"pupil_preprocessed")).name;
        file_path_pupil{end+1} = fullfile(sub_listing(j).folder,sub_listing(j).name,sp_file_name);%#ok<SAGROW> 
    end
end
pupil_to_load = file_path_pupil{contains(file_path_pupil,Session)};
%Find all TCA results
result_folder = "G:\My Drive\SAAAT\SAAAT_PerSession_Result";
listing = dir(result_folder);
tosearch = strcat("rank*",Session,"*.mat");
IndexC = regexp({listing.name}, regexptranslate('wildcard', tosearch));
Index = not(cellfun('isempty',IndexC));
Result_To_Visualize = fullfile(result_folder,listing(Index).name);
% Find all tensor summary
tensor_summary_folder = "D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\tensor_perSession";
tensor_summary_list = dir(tensor_summary_folder);
tosearch = strcat(Session,"*summary");
IndexC = regexp({tensor_summary_list.name}, regexptranslate('wildcard', tosearch));
Index = find(not(cellfun('isempty',IndexC)));
tensor_summary_To_Load = fullfile(tensor_summary_folder,tensor_summary_list(Index).name);
%%
% Load session sumary and tensor summary store earlier

result_session = load(Result_To_Visualize);
M = readtable( pupil_to_load ) ;
load(tensor_summary_To_Load);
temp = strcat(Session,"_rank",num2str(rank),"_component",num2str(component));
axis_time = -999:3000;
axis_trial = 1:numel(Tensor_PerSession_Summary.ds_sig_start);
axis_pupil = 0.55:0.05:numel(Tensor_PerSession_Summary.ds_sig_start)+0.5;
ds_pupil = [];
% Create the concatenated pupil dynamic of each trial
for t = Tensor_PerSession_Summary.ds_sig_start'/1000
    % Cut out pupil dynamic of target trial then down sample it to 20
    % sample/s
    y = M.pupil(M.time < t+3 & M.time > t-1);
    x = 1:numel(y);
    xq = linspace(1,numel(y),20);
    pupil_to_cat = interp1(x,y,xq);
    pupil_to_cat(1) = NaN;
    ds_pupil = cat(1,ds_pupil,pupil_to_cat');
end
ds_pupil = ds_pupil-min(ds_pupil);
ds_pupil = ds_pupil/max(ds_pupil);
% figure()
% 
% p1 = plot(axis_pupil,ds_pupil,'r','LineWidth',2);
% p1.Color(4) = 0.3;
% grid on
% Create the Gaussian filtered normalized response density curve
reaction_times = Tensor_PerSession_Summary.reaction_time( ...
    Tensor_PerSession_Summary.reaction_time > 0);
% map reaction_times to binary time series
reaction_rate = zeros(4000,1);
reaction_times = round(reaction_times*1000);
[N_react, time_bin] = groupcounts(reaction_times);
reaction_rate(time_bin+1000,1) = reaction_rate(time_bin+1000,1) + N_react;
w_gaus = gausswin(Tensor_PerSession_Summary.L_gaus, ...
    Tensor_PerSession_Summary.a_gaus);
reaction_rate_filtered = filter(w_gaus, 1, reaction_rate);
reaction_rate_filtered = reaction_rate_filtered/max(reaction_rate_filtered);
% If folder of the plots of this session does not exist, create one
save_folder = fullfile( ...
    'D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\grabed_component', ...
    Session);
if ~exist(save_folder, 'dir')
   mkdir(save_folder)
end
name_temp = fullfile(save_folder, ...
    strcat(temp, datestr(now,'_yymmdd_HH_MM_SS'))); %#ok
%% Plot selected component
fig1 = figure("color","white",'Position', [50 50 1150 550]);
%set(fig1, 'Visible', 'off');
a = rank;
j = component;
subplot(2,5,[1,2,3])
scatter(find( ...
    Tensor_PerSession_Summary.reward_evidednce==2), ...
    result_session.factors{a,3}( ...
    Tensor_PerSession_Summary.reward_evidednce==2,j),50,'.')
hold on
scatter(find( ...
    Tensor_PerSession_Summary.reward_evidednce==1), ...
    result_session.factors{a,3}( ...
    Tensor_PerSession_Summary.reward_evidednce==1,j),50,'.')
scatter(find( ...
    Tensor_PerSession_Summary.hit_trials==1), ...
    result_session.factors{a,3}( ...
    Tensor_PerSession_Summary.hit_trials==1,j),100,'r|')
if anomaly == 0
    temp3 = max(result_session.factors{a,3}(:,j));
else
    temp3 = maxk(result_session.factors{a,3}(:,j),anomaly+1);
    temp3 = temp3(anomaly+1);
end
p1 = plot(axis_pupil,ds_pupil*temp3,'k');
p1.Color(4) = 0.1;
ylim([0,temp3*1.1])
legend('high reward','low reward','Hit trial','norm pupil size','location','best')
xlabel('trial')
xlim([0,numel(Tensor_PerSession_Summary.ds_sig_start)])
ylabel('trial factor')
% plot t-test for trial-factors between reward vs no reward and large
% pupil(top  25%) vs base pupil(lower 50%)
trial_mean_pupil = [];
% Create the concatenated mean pupil of each trial
for t = Tensor_PerSession_Summary.ds_sig_start'/1000
    % Cut out pupil dynamic of target trial then down sample it to 20
    % sample/s
    y = M.pupil(M.time < t+1 & M.time > t-1);
    trial_mean_pupil = cat(1,trial_mean_pupil,mean(y));
end
tf_rwd = result_session.factors{a,3}( ...
     Tensor_PerSession_Summary.hit_trials==1,j);
tf_norwd = result_session.factors{a,3}( ...
     Tensor_PerSession_Summary.hit_trials==0,j);
tf_highrwd = result_session.factors{a,3}( ...
     Tensor_PerSession_Summary.reward_evidednce==2,j);
tf_lowrwd = result_session.factors{a,3}( ...
     Tensor_PerSession_Summary.reward_evidednce==1,j);
tf_large_pupil = result_session.factors{a,3}( ...
     trial_mean_pupil > quantile(trial_mean_pupil,0.8),j);
tf_base_pupil = result_session.factors{a,3}( ...
     trial_mean_pupil < quantile(trial_mean_pupil,0.5),j);
subplot(2,5,6)
boxplot([tf_rwd;tf_norwd], ...
    [ones(size(tf_rwd));2*ones(size(tf_norwd))], ...
    'Labels',{'Rewarded','No reward'})
[h,p] = ttest2( tf_rwd , tf_norwd ,'Vartype','unequal');
ylim([0,temp3*1.1])
if h == 1
    color = 'r';
else
    color = 'k';
end
title(strcat('p = ', num2str(p)),'color',color)
xlabel(strcat('rank-',num2str(a),',component-',num2str(j)))
ylabel('Trial Factor')
subplot(2,5,7)
boxplot([tf_large_pupil;tf_base_pupil], ...
    [ones(size(tf_large_pupil));2*ones(size(tf_base_pupil))], ...
    'Labels',{'Large Pupil','Base Pupil'})
ylim([0,temp3*1.1])
[h,p] = ttest2( tf_large_pupil , tf_base_pupil ,'Vartype','unequal');
if h == 1
    color = 'r';
else
    color = 'k';
end
xlabel(strcat('rank-',num2str(a),',component-',num2str(j)))
ylabel('Trial Factor')
title(strcat('p = ', num2str(p)),'color',color)
subplot(2,5,8)
boxplot([tf_highrwd;tf_lowrwd], ...
    [ones(size(tf_highrwd));2*ones(size(tf_lowrwd))], ...
    'Labels',{'High Reward','Low reward'})
ylim([0,temp3*1.1])
[h,p] = ttest2( tf_highrwd , tf_lowrwd ,'Vartype','unequal');
if h == 1
    color = 'r';
else
    color = 'k';
end
xlabel(strcat('rank-',num2str(a),',component-',num2str(j)))
ylabel('Trial Factor')
title(strcat('p = ', num2str(p)),'color',color)
subplot(2,5,[4,5])
hold on
patch([axis_time'; nan],[result_session.factors{a,2}(:,j); nan], ...
    [reaction_rate_filtered; nan],[reaction_rate_filtered; nan], ...
    'edgecolor', 'interp','LineWidth',1.5); 
h = colorbar;colormap(winter);
xline(0,'r--')
xlim([-1000 3000])
xlabel('time(ms)')
ylabel('time factor')
ylabel(h, 'Norm Lick Density')
subplot(2,5,[9,10])
bar(result_session.factors{a,1}(:,j))
xlabel('Neuron ID')
ylabel('Neuron factor')
%figure
sgtitle(temp, 'Interpreter', 'none')
%%
saveas(fig1,name_temp)
exportgraphics(fig1,strcat(name_temp,'.jpg'),'Resolution',300)
fig_output = fig1;
disp("Done")









