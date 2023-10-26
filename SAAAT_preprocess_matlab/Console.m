%%
all = [7 6 9 9 0];
Session = strcat('W333',num2str(all(1)),'_',num2str(all(2)));
rank = all(3);
component = all(4);
anomaly = all(5);
fig_output = TCA_Component_Visualization_V2( ...
    Session, ...
    rank, ...
    component, ...
    anomaly);
%%

load(['G:\My Drive\SAAAT\SAAAT_session_signal_tensor_V2_Result\' ...
    'rank1to13_6rep_ncp_hals__5_W3336_5_Spikes_signal_tensor_V2_NTK_230912_14_39_06.mat'])

neuron_factor_example = factors{5, 1}(:,2);
save("neuron_factor_example.mat","neuron_factor_example")