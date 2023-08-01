%%
clear
session = 'W3333 session3';
%%
sp_times_converted = sp_times/PARAMS.Fs;
spikes_2plt = 50000;
StartTime_2plt = sp_times(numel(sp_times)-spikes_2plt);
tc_2plt = PARAMS.TcStart_LV( ...
    PARAMS.TcStart_LV>StartTime_2plt)'...
    /PARAMS.Fs;
sig_2plt = PARAMS.Target_Start( ...
    PARAMS.Target_Start>StartTime_2plt)...
    /PARAMS.Fs;
%%
a=12;
figure
subplot(1,2,1)
scatter(sp_times_converted(numel(sp_times)-spikes_2plt:end-30000), ...
    spike_clusters(numel(sp_times)-spikes_2plt:numel(sp_times)-30000), ...
    "|")
title(strcat(session,': Matched by start of each vector'))
xline(tc_2plt(1:a),"k--","LineWidth",1)
xline(sig_2plt(1:a),"r--","LineWidth",1)
legend("Spikes","Tone CLoud Start","Signal Start")
xlabel('time(s)')
ylabel('Neuron')
subplot(1,2,2)
hold on
scatter(sp_times_converted(numel(sp_times)-spikes_2plt:end-30000), ...
    spike_clusters(end-spikes_2plt:end-30000), ...
    "|")
title(strcat(session,': Matched by end of each vector'))
xline(tc_2plt(1:a),"k--","LineWidth",1)
xline(sig_2plt(1:a),"r--","LineWidth",1)
legend("Spikes","Tone CLoud Start",'','','',"Signal Start")
xlabel('time(s)')
ylabel('Neuron')
%%
clear
spikeTimes = readNPY(['spike_times.npy']);
clusterID = readNPY(['spike_clusters.npy']);