clear
%% Find the path of all pupil data per session per animal
listing = dir('D:\KunpengY\AE_trained_behaving\preprocess');
result_folder = "G:\My Drive\SAAAT\SAAAT_session_signal_tensor_V2_Result";
plot_folder = 'D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\session_plot_v2';
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
%Find all TCA results
listing = dir(result_folder);
tosearch = strcat("rank*.mat");
IndexC = regexp({listing.name}, regexptranslate('wildcard', tosearch));
Index = find(not(cellfun('isempty',IndexC)));
Result_To_Visualize = {listing(Index).name}';
% Find all tensor summary
tensor_summary_folder = "D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\tensor_perSession";
tensor_summary_list = dir(tensor_summary_folder);
tensor_summary_list = {tensor_summary_list( ...
    contains({tensor_summary_list.name},'summary')).name};
%%
% Load session sumary and tensor summary store earlier
%set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')
wbar = waitbar(0,'Start Processing...');
for i = 1:numel(Result_To_Visualize)-1
    waitbar((i-1)/numel(Result_To_Visualize), ...
                wbar,strcat('Processing session:',num2str(i)));
    result_session = load(fullfile(result_folder,Result_To_Visualize{i}));
    n = 3; repeats = 3;
    temp_axis = repmat(result_session.rank,width(result_session.error), 1);
    temp_axis = reshape(temp_axis,[numel(temp_axis),1]);
    temp_axis2 = repmat(result_session.rank,width(result_session.simi), 1); 
    temp_axis2 = reshape(temp_axis2,[numel(temp_axis2),1]);
    %fig_err_simi = figure("color","white",'Position', [50 50 1150 650]);
    fig1 = figure("color","white",'Position', [50 50 1150 650]);
    set(fig1, 'Visible', 'off');
    subplot(1,2,1)
    grid on
    hold on
    plot(result_session.rank,min(result_session.error,[],2),'LineWidth',1.5)
    scatter(temp_axis,reshape(result_session.error', ...
        [numel(result_session.error),1]),60,'+')
    legend('minimun error','error of each rep', ...
        'Location','northoutside')
    xlabel('# components')
    ylabel('Error of reconstruction')
    subplot(1,2,2)
    grid on
    hold on
    plot(result_session.rank,mean(result_session.simi,2),'LineWidth',1.5)
    scatter(temp_axis2,reshape(result_session.simi', ...
        [numel(result_session.simi),1]),60,'+')
    legend('mean similarity','similarity of each rep', ...
        'Location','northoutside')
    ylim([0.7,1])
    xlabel('# components')
    ylabel('Similarity comparing to best fit')
    sgtitle(strcat('TCA error and similarity curve (', ...
        num2str(width(result_session.error)),' replicates each rank)'))
    
    result_file_name = char(fullfile(result_folder,Result_To_Visualize{i}));
    temp = strfind(result_file_name,'_');
    temp = temp(temp>strfind(result_file_name,'rank'));
    temp = result_file_name(temp(6)+1:temp(8));
    pupil_to_read = file_path_pupil{contains(file_path_pupil,temp)};
    M = readtable( pupil_to_read ) ;
    load(fullfile(tensor_summary_folder, ...
        tensor_summary_list{contains(tensor_summary_list,temp)}));
    
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
    save_folder = fullfile(plot_folder, temp);
    if ~exist(save_folder, 'dir')
       mkdir(save_folder)
    end
    name_temp = fullfile(save_folder, ...
        strcat('error_similarity', datestr(now,'_yymmdd_HH_MM_SS'))); %#ok
    saveas(fig1,name_temp)
    saveas(fig1,name_temp,'jpeg')
    % Plot each rank
    for a = double(result_session.rank)
        waitbar((i-1)/numel(Result_To_Visualize)+ ...
            a/numel(result_session.rank)/numel(Result_To_Visualize), ...
            wbar,strcat('Processing session:',num2str(i)));

        fig2 = figure("color","white",'Position', [50 50 1250 150*a+50]);
        set(fig2, 'Visible', 'off');
        for j = 1:a
            subplot(a,4,j*4-3)
            %figure
            hold on
            bar(result_session.factors{a,1}(:,j))
            xlabel('Neuron ID')
            ylabel('Neuron factor')
            subplot(a,4,j*4-2)
            hold on
            patch([axis_time'; nan],[result_session.factors{a,2}(:,j); nan], ...
                [reaction_rate_filtered; nan],[reaction_rate_filtered; nan], ...
                'edgecolor', 'interp','LineWidth',1.5); 
            h = colorbar;colormap(winter);
            xline(0,'r--')
            xlim([-1000 3000])
            xlabel('time(ms)')
            ylabel('time factor')
            ylabel(h, 'Lick Density')
            subplot(a,4,[j*4-1,j*4])
            scatter(find( ...
                Tensor_PerSession_Summary.reward_evidednce==2), ...
                result_session.factors{a,3}( ...
                Tensor_PerSession_Summary.reward_evidednce==2,j),'.')
            hold on
            scatter(find( ...
                Tensor_PerSession_Summary.reward_evidednce==1), ...
                result_session.factors{a,3}( ...
                Tensor_PerSession_Summary.reward_evidednce==1,j),'.')
            scatter(find( ...
                Tensor_PerSession_Summary.hit_trials==1), ...
                result_session.factors{a,3}( ...
                 Tensor_PerSession_Summary.hit_trials==1,j),'r|')
            p1 = plot(axis_pupil,ds_pupil*max(result_session.factors{a,3}(:,j)),'k');
            p1.Color(4) = 0.1;
            legend('high reward','low reward','Hit trial','norm pupil size','location','eastoutside')
            xlabel('trial')
            ylabel('trial factor')
            
        end
        name_temp = fullfile(save_folder, ...
            strcat('factor_analysis_rank',num2str(a), datestr(now,'_yymmdd_HH_MM_SS'))); %#ok
        saveas(fig2,name_temp)
        exportgraphics(fig2,strcat(name_temp,'.jpg'),'Resolution',300)

        fig3 = figure("color","white",'Position', [50 50 200*a+50 750]);
        set(fig3, 'Visible', 'off');
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
        for j = 1:a
            % tf: trial factor
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
            subplot(3,a,j)
            boxplot([tf_rwd;tf_norwd], ...
                [ones(size(tf_rwd));2*ones(size(tf_norwd))], ...
                'Labels',{'Rewarded','No reward'})
            [h,p] = ttest2( tf_rwd , tf_norwd ,'Vartype','unequal');
            if h == 1
                color = 'r';
            else
                color = 'k';
            end
            title(strcat('p = ', num2str(p)),'color',color)
            xlabel(strcat('rank-',num2str(a),',component-',num2str(j)))
            ylabel('Trial Factor')
            subplot(3,a,a+j)
            boxplot([tf_large_pupil;tf_base_pupil], ...
                [ones(size(tf_large_pupil));2*ones(size(tf_base_pupil))], ...
                'Labels',{'Large Pupil','Base Pupil'})
            [h,p] = ttest2( tf_large_pupil , tf_base_pupil ,'Vartype','unequal');
            if h == 1
                color = 'r';
            else
                color = 'k';
            end
            xlabel(strcat('rank-',num2str(a),',component-',num2str(j)))
            ylabel('Trial Factor')
            title(strcat('p = ', num2str(p)),'color',color)
            subplot(3,a,2*a+j)
            boxplot([tf_highrwd;tf_lowrwd], ...
                [ones(size(tf_highrwd));2*ones(size(tf_lowrwd))], ...
                'Labels',{'High Reward','Low reward'})
            [h,p] = ttest2( tf_highrwd , tf_lowrwd ,'Vartype','unequal');
            if h == 1
                color = 'r';
            else
                color = 'k';
            end
            xlabel(strcat('rank-',num2str(a),',component-',num2str(j)))
            ylabel('Trial Factor')
            title(strcat('p = ', num2str(p)),'color',color)

        end
        sgtitle(strcat('Box plot and t-test for trial factors of different groups for rank-',num2str(a)))
        name_temp = fullfile(save_folder, ...
            strcat('trial_factor_Ttest_rank',num2str(a), datestr(now,'_yymmdd_HH_MM_SS'))); %#ok
        saveas(fig3,name_temp)
        exportgraphics(fig3,strcat(name_temp,'.jpg'),'Resolution',300)
    end
    
    

end
close(wbar)

