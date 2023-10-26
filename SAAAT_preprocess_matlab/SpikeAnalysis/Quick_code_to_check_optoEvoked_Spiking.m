%% quick code for looking at laser induced spiking
%% it will plot unit wise rastor plot and save it as  png



%% read the digital input to extract intan syncs
main_path=pwd;
cd 'Z:\Users\hsrivastava\EfferentsTest\Efferents'
load('positions.mat')
cd(main_path)
clearvars -except p
pre=.02;post=.02;
[ConvertedData,cc]=convertTDMS2mat_pipeline;

AllParams=ConvertedData.Data.MeasuredData;
for kk=1:length(AllParams)
    tf=strcmp(AllParams(kk).Name,'Untitled/LaserTP');
    if tf==1
        PARAMS.LaserTime=AllParams(kk).Data;
    end
    tf1=strcmp(AllParams(kk).Name,'Untitled/IntanSync');
    if tf1==1
        PARAMS.LV_Sync=AllParams(kk).Data;

    end

    tf6=strcmp(AllParams(kk).Name,'Settings');
    if tf6==1
        StimProps=AllParams(kk).Property;
    end
end



for kk=1:length(StimProps)
    tf2=strcmp(StimProps(kk).Name,'Sound Configuration.Stimulus Config.Laser N');
    if tf2==1
        PARAMS.NumPulses=str2num(StimProps(kk).Value);
    end

    tf2=strcmp(StimProps(kk).Name,'Sound Configuration.Stimulus Config.Laser D(ms) ');
    if tf2==1
        PARAMS.PulseDur=str2num(StimProps(kk).Value);
    end

    tf2=strcmp(StimProps(kk).Name,'Sound Configuration.Stimulus Config.Laser IPI(s)');
    if tf2==1
        PARAMS.InterPulseDur=str2num(StimProps(kk).Value);
    end
end

move_one_folder_forward

fileinfo = dir('digitalin.dat');
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen('digitalin.dat', 'r');
digital_word = fread(fid, num_samples, 'uint16');
fclose(fid);
ch=1;
digital_input_ch = (bitand(digital_word, 2^ch) > 0);

diff_signal=diff(digital_input_ch);
Intan_Sync_onIntan=find(diff_signal==-1);


LVSync=PARAMS.LV_Sync;
LVSync(find(LVSync>=2.5))=5;
LVSync(find(LVSync<2.5))=0;

diff_signalLV=diff(LVSync);
IntanSync_onLV=find(diff_signalLV==-5);

% remove excess pulse
if length(IntanSync_onLV)<length(Intan_Sync_onIntan)
    temp=Intan_Sync_onIntan(1:length(IntanSync_onLV));
    Intan_Sync_onIntan=[];
    Intan_Sync_onIntan=temp;
end


IntanFs=30000;
LV_Fs=30003;

%% read spike times
format long
spikeTimes_SampleNum = readNPY(['spike_times.npy']);
clusterID = readNPY(['spike_clusters.npy']);
clusterGroups =tdfread('cluster_KSLabel.tsv');
cl_id=clusterID;
try
    clusterinfo=tdfread('cluster_info.tsv');
    channel_ID=tdfread('cluster_info.tsv');
    depths_with_clus_id(:,1)=channel_ID.depth;
    depths_with_clus_id(:,2)=channel_ID.cluster_id;
    depths_with_clus_id(:,3)=channel_ID.ch;
end
spikeTimes_SampleNum=double(spikeTimes_SampleNum);

spikeTimes_Seconds=spikeTimes_SampleNum./IntanFs;
Intan_Sync_onIntan_Seconds=Intan_Sync_onIntan./IntanFs;



pp='amplifier.dat';
nChansInFile = 64;  % neuropixels phase3a, from spikeGLX
d = dir(pp);
nSamps = d.bytes/2/nChansInFile;
mmf = memmapfile(pp, 'Format', {'int16', [nChansInFile nSamps], 'x'});
st = readNPY('spike_times.npy');

[tempWF]=getAllSpikesWaveforms(st,mmf,cl_id);




%% convert spike times on Intan board into LV board
sp_times=[];  %% this takes time to run depending upon number of pulses
relevant_waveforms=[];
disp('converting spiketimes on LV time')
for num_sync=1:length(Intan_Sync_onIntan)-1


    curr_sp=find(spikeTimes_Seconds>=Intan_Sync_onIntan_Seconds(num_sync) & spikeTimes_Seconds<Intan_Sync_onIntan_Seconds(num_sync+1));

    if ~isempty(curr_sp)

        sp_times = [sp_times; IntanSync_onLV(num_sync) + (spikeTimes_Seconds(curr_sp)-Intan_Sync_onIntan_Seconds(num_sync)).*LV_Fs;];

        relevant_waveforms=[relevant_waveforms; tempWF(curr_sp,:)];
    end
    clear curr_sp
end



%% find laser onset times

LaserTemp=PARAMS.LaserTime;
LaserTemp(find(LaserTemp>=2.5))=5;
LaserTemp(find(LaserTemp<2.5))=0;

diff_Laser=diff(LaserTemp);
Laser_Onset=find(diff_Laser==5);

if length(Laser_Onset)~=PARAMS.NumPulses;
    disp('Incorrect number of pulses')
    req_gap=median(diff(Laser_Onset));
    Laser_Onset_copy=Laser_Onset;
    takeDiff=diff(Laser_Onset);
    takeDiff_2(2:length(Laser_Onset))=takeDiff;
    takeDiff_2(1)=Laser_Onset(1);

    Laser_Onset_copy(find(takeDiff_2<(.9.*req_gap)))=[];
    Laser_Onset=[];
    Laser_Onset=Laser_Onset_copy;
    if length(Laser_Onset)==PARAMS.NumPulses;
        disp('Correction done')
    else error('Incorrect number of pulses')
    end
end




%% make raster plots


stim_ind=Laser_Onset;
n_neurons = length(unique(cl_id));
ui = unique(cl_id);
for uu = 1:n_neurons
    this_ui = find(cl_id(1:length(sp_times)) == ui(uu));
    ui_spike_sample = sp_times(this_ui);
    ui_spike_time{uu} = (ui_spike_sample/30003);
end



%% read raw data for spike waveforms
for kk=1:length(ui_spike_time) % run this from unit number where you expect response, in this session it was around 500
    kk
     java.lang.System.gc()
    sptimes=ui_spike_time{kk};
  
    find_figure('Raster');clf
   
    axes('position',p(:,1))
    box on
    hold on
    cnty=0;
    spt1=[];YY=[];spikewaveF=[];new_waveforms=[];
    savename=strcat('Unit_',num2str(kk));
    for jj=1:length(stim_ind)
        this_stim=stim_ind(jj)./30003;
        all=find(sptimes>(this_stim-pre));
        new_sp=sptimes(all);
        new_waveforms=relevant_waveforms(all,:);
        new_waveforms(find(new_sp>(this_stim+post)),:)=[];
        new_sp(find(new_sp>(this_stim+post)))=[];
        spt=new_sp-this_stim;

        if length(spt)>0

            yy=jj.*ones(1,length(spt));
            spt1=[spt1;spt];
            YY=[YY yy];
            spikewaveF=[spikewaveF;new_waveforms];



            clear spt new_sp all


        end
    end

    plot(spt1,YY,'.k')
    xlim([-pre post])
    ylim([0 length(Laser_Onset)])

    x_points = [0, 0, (PARAMS.PulseDur./1000), (PARAMS.PulseDur./1000)];
    y_points = [0, length(Laser_Onset), length(Laser_Onset), 0];
    color = [0, 0, 1];

    hold on;
    a = fill(x_points, y_points, color);
    a.FaceAlpha = 0.1;
   
    hold on
    xline(0)
    title(strcat('Cluster# = ',num2str(depths_with_clus_id(kk,2))))
    ylabel('trial #')
    xlabel('Time (sec)')



     binsize=.001; % binsize in seconds
        totTime=600;
       num_bins=0:binsize:totTime;
        SpikeTrain=zeros(1,length(num_bins)-1);
      
       
            spdata=sptimes;
            for timebins=1:length(spdata)
                SpikeTrain(1,fix(spdata(timebins)/binsize)+1)=SpikeTrain(1,fix(spdata(timebins)/binsize)+1) +1;
            end
            clear spdata

              pre_bins=.02/binsize; % spont bins
        post_bins=.02/binsize;

             for numLaser=1:length(Laser_Onset)
            
            curr_stim_time=Laser_Onset(numLaser)/30003;
            StimLocked_Sp(:,numLaser)=SpikeTrain(:,find(num_bins>curr_stim_time,1,'first')-pre_bins:find(num_bins>curr_stim_time,1,'first')+post_bins);
        end

        axes('position',p(:,3))
%         shadedErrorBar(-pre:binsize:post,mean(StimLocked_Sp,2),(std(StimLocked_Sp,[],2))./sqrt(size(StimLocked_Sp,2)))
            plot(-pre:binsize:post,(mean(StimLocked_Sp,2))./binsize,'LineWidth',2,'Color','k')
        hold on
          x_points = [0-.001, 0-.001, (PARAMS.PulseDur./1000)-.001, (PARAMS.PulseDur./1000)-.001];
    y_points = [0, (max(mean(StimLocked_Sp,2)))./binsize, (max(mean(StimLocked_Sp,2)))./binsize, 0];
    color = [0, 0, 1];

    hold on;
    a = fill(x_points, y_points, color);
    a.FaceAlpha = 0.1;
    ylabel('Firing Rate (Hz)')
    xlabel('Time (sec)')
    %saveas(gcf,savename,'png')

    axes('position',p(:,2))
    hold on


  
    baselinecorr0=mean(spikewaveF(:,1:25),2);
     baselinecorr1=repmat(baselinecorr0,1,61);
     basecorrected=spikewaveF-baselinecorr1;

    plot(basecorrected','color',[.5 .5 .5])

   hold on
   plot(mean(basecorrected),'LineWidth',2,'Color','k')

   clear tempWF extractST theseST nWFsamps

    pause

end




























