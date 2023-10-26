clear
%%
load('SAAAT_session_summary_073123.mat')
save_folder = 'D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\session_performance_analysis';
if ~exist(save_folder, 'dir')
   mkdir(save_folder)
end
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
%%
wbar = waitbar(0,'Start Processing...');
for i = 1:numel(session_summary.file_path_PARAMS) %!
    close
    waitbar((i-1)/numel(session_summary.file_path_PARAMS), ...
                wbar,strcat('Processing session:',num2str(i)));
    fig1 = figure("color","white",'Position', [50 50 750 750]);
    set(fig1, 'Visible', 'off');
    subplot(2,1,1);
    hold on
    temp3 = session_summary.reward_evidednce_all{i};
    reward_evidednce_block = zeros(size(session_summary.reward_evidednce_all{i}));
    for j  = 1:numel(temp3)
        if temp3(j) == 1
            reward_evidednce_block(j:end) = 0;
        elseif temp3(j) == 2
            reward_evidednce_block(j:end) = 1;
        end
    end
    s1 = area(reward_evidednce_block,"LineStyle","none");
    alpha(s1,.4)
    temp = session_summary.hit_trials_all{i};
    MW_average = cat(1,ones(10,1),zeros(9,1));
    Hit_ratio_MW_average = filter(MW_average, sum(MW_average), temp);
    plot(Hit_ratio_MW_average,"r","LineWidth",1.5)
    xlabel("Trial number (Only include miss and hit trials)")
    ylabel("Hit to (Hit+Miss) ratio")
    temp2 = strfind(session_summary.file_path_Spikes{i},'Spikes');
    session_name = session_summary.file_path_Spikes{i}(temp2-8:temp2-1);
    title(strcat(session_name,"Hit to (Hit+Miss) ratio over blocks"),'Interpreter',"None")
    legend(["High Context Block","Hit to Hit+Miss Ratio"])
    temp3 = strrep(session_name,'_','\');
    
    cd(folder_path_EphysRaw{contains(folder_path_EphysRaw, ...
        temp3)})
    load("rez.mat","rez");
    load(session_summary.file_path_PARAMS{i})
    Fs = PARAMS.Fs;
    batch_samples = rez.ops.NT;
    dshift = rez.dshift;
    dshift_time = (1:height(dshift))'* batch_samples/Fs;
    dshift_sudo_time = [0];
    for j = 1:numel(session_summary.ds_sig_start_all{i})
        index_dshift = round(session_summary.ds_sig_start_all{i}(j)/(1000*batch_samples/Fs));
        to_concat = linspace(dshift_sudo_time(end),j,index_dshift+1-numel(dshift_sudo_time));
        dshift_sudo_time = cat(1,dshift_sudo_time,to_concat(2:end)');
    end
    subplot(2,1,2)
    min_dshift = 1.1*min(dshift,[],"all");
    max_dshift = 1.1*max(dshift,[],"all");
    hold on
    s2 = area(reward_evidednce_block*(max_dshift-min_dshift)+min_dshift,"LineStyle","none");
    alpha(s2,.2)
    plot(dshift_sudo_time,dshift(1:numel(dshift_sudo_time),:))
    % Create a sudo time axis of trial
    ylim([min_dshift,max_dshift])
    ylabel('um')
    xlabel("Trial number (Only include miss and hit trials)")
    title("Estimated drifting through trials")
    save_path = fullfile(save_folder,strcat(session_name,'_session_hit_ratio'));
    saveas(fig1,save_path)
    exportgraphics(fig1,strcat(save_path,'.jpg'),'Resolution',200)
end
close all
close(wbar)