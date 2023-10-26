clear
%% Find the path of all Tensor in path folder to concat 
target_folder = 'D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\session_2concat_v3';
listing = dir(target_folder);
listing = listing(~ismember({listing.name},{'.','..'}));
tensor_file_name = {listing(contains({listing.name},"NTK")).name};
save(fullfile(target_folder, 'concated_tensor_order'),"tensor_file_name")
%%
new_N = 0;
new_K = 0;
old_T = 0;
for i = 1: numel(tensor_file_name)
    load(fullfile(target_folder,tensor_file_name{i}),"Session_mask")
    new_N = new_N + size(Session_mask,1);
    new_K = new_K + size(Session_mask,3);
    old_T = size(Session_mask,2);
    clear Session_mask
end
%%
pooled_mask = zeros(new_N, old_T, new_K,'logical');
current_N = 0;
current_K = 0;
for i = 1: numel(tensor_file_name)
    load(fullfile(target_folder,tensor_file_name{i}),"Session_mask")
    Session_mask = logical( Session_mask );
    Next_N = current_N+size(Session_mask,1);
    Next_K = current_K+size(Session_mask,3);
    disp(strcat('N = ', num2str(Next_N),', K = ', num2str(Next_K)))
    pooled_mask(current_N+1:Next_N, :, current_K+1:Next_K) = Session_mask;
    current_N = Next_N;
    current_K = Next_K;
    clear Session_mask
end
disp("Finished concating mask")
save_path1 = fullfile(target_folder, ...
    strcat("pooled_mask_",datestr(now,'yymmdd_HH_MM_SS')));
save(save_path1,"pooled_mask","-v7.3")
disp("concat mask saved")
clear pooled_mask
%%
pooled_tensor_NTK = zeros(new_N, old_T, new_K,'single');
current_N = 0;
current_K = 0;
for i = 1: numel(tensor_file_name)
    load(fullfile(target_folder,tensor_file_name{i}),"Session_tensor")
    Next_N = current_N+size(Session_tensor,1);
    Next_K = current_K+size(Session_tensor,3);
    disp(strcat('N = ', num2str(Next_N),', K = ', num2str(Next_K)))
    pooled_tensor_NTK(current_N+1:Next_N, :, current_K+1:Next_K) = Session_tensor;
    current_N = Next_N;
    current_K = Next_K;
    clear Session_tensor
end
disp("Finished concating tensor")
save_path1 = fullfile(target_folder, ...
    strcat("pooled_tensor_NTK_", datestr(now,'yymmdd_HH_MM_SS')));
save(save_path1,"pooled_tensor_NTK","-v7.3")
disp("concat tensor saved")
clear pooled_tensor_NTK