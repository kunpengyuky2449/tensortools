clear
%% Load all stored fitted params and related data into a cell array in order of rank
result_folder = "G:\My Drive\SAAAT\SAAAT_PerSession_Result";
listing = dir(result_folder);
tosearch = strcat("rank*.mat");
IndexC = regexp({listing.name}, regexptranslate('wildcard', tosearch));
Index = find(not(cellfun('isempty',IndexC)));
Result_To_Visualize = {listing(Index).name}';
%%
% Load session sumary and tensor summary store earlier
set(0,'DefaultFigureWindowStyle','docked')
for i = 1:numel(Result_To_Visualize)
    result_session = load(fullfile(result_folder,Result_To_Visualize{i}));
    n = 3; repeats = 3;
    temp_axis = repmat(result_session.rank,width(result_session.error), 1);
    temp_axis = reshape(temp_axis,[numel(temp_axis),1]);
    temp_axis2 = repmat(result_session.rank,width(result_session.simi), 1);
    temp_axis2 = reshape(temp_axis2,[numel(temp_axis2),1]);
    %fig_err_simi = figure("color","white",'Position', [50 50 1150 650]);
    figure("color","white",'Position', [50 50 1150 650]);
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
    

end

%% Plot error drop curve and similarity change curve
plot_summary.ranks_errors = [];
plot_summary.ranks_simis = [];
plot_summary.ranks_fits = [];
plot_summary.errors = [];
plot_summary.errors_mean = [];
plot_summary.simi = [];
plot_summary.simi_mean  = [];
for i = 1:numel(All_fits)
    plot_summary.ranks_fits = cat(1, ...
        plot_summary.ranks_fits, All_fits{i}.rank);
    plot_summary.errors = cat(1, ...
        plot_summary.errors, All_fits{i}.error');
    plot_summary.simi = cat(1, ...
        plot_summary.simi, All_fits{i}.simi');
    plot_summary.ranks_errors = cat(1, ...
        plot_summary.ranks_errors, repmat(All_fits{i}.rank,size(All_fits{i}.error')));
    plot_summary.ranks_simis = cat(1, ...
        plot_summary.ranks_simis, repmat(All_fits{i}.rank,size(All_fits{i}.simi')));
    plot_summary.errors_mean = cat(1, ...
        plot_summary.errors_mean, mean(All_fits{i}.error));
    plot_summary.simi_mean = cat(1, ...
        plot_summary.simi_mean, mean(All_fits{i}.simi));
end
%%
fig = figure("color","white",'Position', [50 50 1150 650]);
subplot(1,2,1)
grid on
hold on
plot(plot_summary.ranks_fits,plot_summary.errors_mean,'LineWidth',1.5)
scatter(plot_summary.ranks_errors,plot_summary.errors,60,'+')
legend('mean error','error of each rep', ...
    'Location','northoutside')
xlabel('# components')
ylabel('Error of reconstruction')
subplot(1,2,2)
grid on
hold on
plot(plot_summary.ranks_fits,plot_summary.simi_mean,'LineWidth',1.5)
scatter(plot_summary.ranks_simis,plot_summary.simi,60,'+')
legend('mean similarity','similarity of each rep', ...
    'Location','northoutside')
xlabel('# components')
ylabel('Similarity compaaring to best fit')
sgtitle('TCA error and similarity curve (4 replicates each rank)')
name_temp = strcat('mat_plot\plot_error_simi_rank',num2str(min(plot_summary.ranks_fits)),'to', ...
    num2str(max(plot_summary.ranks_fits)), datestr(now,'_yymmdd_HH_MM_SS')); %#ok
saveas(fig,name_temp)
saveas(fig,name_temp,'jpeg')
close
%% Get axis for later mapping
plot_summary.neuron_depth = cat(1,tensor_summaray.neuron_depth_per_session{:});
reaction_times = [];
hit_binary = [];
for i = tensor_summaray.session_index
    new_block = [1;find(session_summary.reward_evidednce_all{i} > 0)];
    trials_selected = [];
    for j = 1:numel(tensor_summaray.trials_per_block)
        trials_selected = cat(2,trials_selected, ...
            new_block(j):(new_block(j)+tensor_summaray.trials_per_block(j)-1));
    end
    reaction_times = cat(1, reaction_times, session_summary.reaction_time_all{i}(trials_selected));
    hit_binary = cat(2,hit_binary, session_summary.hit_trials_all{i}(trials_selected));
end
reaction_times = reaction_times(reaction_times > 0);
% map reaction_times to binary time series
reaction_rate = zeros(3500,1);
reaction_times = round(reaction_times*1000);
[N_react, time_bin] = groupcounts(reaction_times);
reaction_rate(time_bin,1) = reaction_rate(time_bin,1) + N_react;
w_gaus = gausswin(session_summary.L_gaus, session_summary.a_gaus);
reaction_rate_filtered = filter(w_gaus, 1, reaction_rate);
plot_summary.reaction_rate_filtered = reaction_rate_filtered/max(reaction_rate_filtered);
plot_summary.trial_time_axis = (-499:3000)';
hit_rate = mean(hit_binary,2);
tbd = cumsum(tensor_summaray.trials_per_block);
hit_rate(tbd(1:end-1)+1) = hit_rate(tbd(1:end-1)+2);
w_gaus2 = gausswin(13, 4);
hit_rate_filtered = filter(w_gaus2, 1, hit_rate);
plot_summary.hit_rate_filtered = hit_rate_filtered/max(hit_rate_filtered);
temp = [0; tbd];
plot_summary.block_low_index = [1:38,58:84,107:127,145:164];
plot_summary.block_high_index = [39:57,85:106,128:144,165:179];
%%
length = numel(plot_summary.neuron_depth);
color_p = repmat([1,0,0],[length,1]);
color_p(:,1) = color_p(:,1) - round( ...
    plot_summary.neuron_depth/max(plot_summary.neuron_depth),2);
color_p(:,3) = color_p(:,3) + round( ...
    plot_summary.neuron_depth/max(plot_summary.neuron_depth),2);
%% Component analysis
a = 8;
fig2 = figure("color","white",'Position', [50 50 1150 150*a+50]);
for i = 1:a
    subplot(a,3,i*3-3+1)
    %figure
    hold on
    bar(All_fits{a}.factors{1}(:,i))
    %scatter(plot_summary.neuron_depth, All_fits{a}.factors{1}(:,i),'.')
    xlabel('Neuron depth')
    ylabel('Neuron factor')
    subplot(a,3,i*3-3+2)
    plot(plot_summary.trial_time_axis,All_fits{a}.factors{2}(:,i))
    xline(0,'r--')
    xlim([-500 3500])
    xlabel('time(ms)')
    ylabel('time factor')
    subplot(a,3,i*3-3+3)
    hold on
    temp = 1:179;
    scatter(temp(plot_summary.block_low_index), ...
        All_fits{a}.factors{3}(plot_summary.block_low_index,i),'.')
    scatter(temp(plot_summary.block_high_index), ...
        All_fits{a}.factors{3}(plot_summary.block_high_index,i),'.')
    plot(temp,plot_summary.hit_rate_filtered*max(All_fits{a}.factors{3}(:,i)))
    legend('low reward','high reward','norm hit rate','location','eastoutside')
    xlabel('trial')
    ylabel('trial factor')
end
%%
name_temp = strcat('mat_plot\rank',num2str(a),'_component_analysis', datestr(now,'_yymmdd_HH_MM_SS')); %#ok
saveas(fig2,name_temp)
saveas(fig2,name_temp,'jpeg')
%% Focus on analyzing neural factor
fig3 = figure("color","white",'Position', [50 50 1150 300]);
a = 5;
temp2 = 1;
for i = [2,3,5,1,4]
    subplot(1,5,temp2)
    temp2 = temp2+1;
    hold on
    for j = 1:length
    bar( j, All_fits{a}.factors{1}(j,i), ...
        'FaceColor', color_p(j,:), ...
        'EdgeColor',color_p(j,:))
    end
    %scatter(plot_summary.neuron_depth, All_fits{a}.factors{1}(:,i),'.')
    xlabel('Neuron depth')
    ylabel('Neuron factor')
end
%%
name_temp = strcat('mat_plot\rank',num2str(a),'_trial_factors_depth', datestr(now,'_yymmdd_HH_MM_SS')); %#ok
saveas(fig3,name_temp)
saveas(fig3,name_temp,'jpeg')
%%
fig3 = figure("color","white",'Position', [50 50 1150 150+50]);
a = 5;
temp3 = [2,3,5,1];
for i = 1:3
    subplot(1,3,i)
    last_component = All_fits{a}.factors{1}(:,temp3(i));
    last_component = last_component/mean(last_component);
    new_component = All_fits{a}.factors{1}(:,temp3(i+1));
    new_component = new_component/mean(new_component);
    diff = (new_component./last_component);
    hold on
    scatter( last_component , new_component ,[],[.5 .5 .5], ...
        '.' )
    scatter( last_component(diff<1/2) , ...
        new_component(diff<1/2) ,[],[1 0 0], ...
        '.' )
    scatter( last_component(diff < 1.3 & diff > 1/1.3) , ...
        new_component(diff < 1.3 & diff > 1/1.3) ,[],[0.2 0.5 0.2], ...
        '.' )
    scatter( last_component(diff>2) , ...
        new_component(diff>2) ,[],[0 0 1], ...
        '.' )
    perc = [sum(diff<1/2)/numel(diff), ...
        sum((diff < 1.3 & diff > 1/1.3))/numel(diff), ...
        sum(diff>2)/numel(diff)]*100;
    legend('',strcat(num2str(perc(1)),'%'),...
        strcat(num2str(perc(2)),'%'),...
        strcat(num2str(perc(3)),'%'))
end
%%
