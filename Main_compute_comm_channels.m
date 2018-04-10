% find orthogonal Quorum sensing channels based on user specified
% requirements
% ZAT 2018 Imperial College London

clear variables
close all
clc

%% user variables
filename    = 'AHL receiver devices parameters 6 dev 6 ind.xlsx';
sheet       = 2;
% specify number of communication channels needed
channel_num = 2;
% specify device fold activation threshold needed
act_th      = 2;
% specify crosstalk device[s] fold activation threshold needed
cross_th    = 2;
% specify minimum [AHL] in Molar / 10^X
min_AHL_concentration = -15;
% specify maximum [AHL] in Molar / 10^X
max_AHL_concentration = -4;

% first or all solutions (is) are needed
config.solution    = 'all'; % first or all
% the results are shown on the screen or written into a file
config.show_result = 'screen'; % screen or file
%screen output verbose or minimal?
config.output      = 'verbose'; % verbose or minimal

%% read excel
[num,txt] = xlsread(filename,sheet);
devices   = {txt{2:6:end,1}};
AHLs      = txt(2:7,2)';

num_of_devices = size(devices,2);
num_of_AHLs    = size(AHLs,2);
data_mtx_size  = num_of_devices*num_of_AHLs;


if strcmpi(config.show_result,'screen')
    fileID = 1; % ID for the screen
else
    fileID = fopen(sprintf('results_d_%d_ahl_%d_ch_%d.txt',num_of_devices,num_of_AHLs,channel_num),'w');
end

fprintf('Looking for %d channels with %d devices and %d AHLs\n',channel_num,num_of_devices,num_of_devices);
tmp_cell = cell(data_mtx_size,1);
%% build raw data cell array, each cell contains the four measurement data
for k=1:data_mtx_size
    tmp_cell{k} = num(k,:);
end
data_mtx = reshape(tmp_cell,num_of_AHLs,num_of_devices)'; % order of argument: because of transpose

%% convert the four measurement data in each cells to dose response
%AHL inducer concentration range space
conc_range    = linspace(min_AHL_concentration,max_AHL_concentration,max_AHL_concentration-min_AHL_concentration+1);

dose_curve    = @(x,conc_range) x(1) + ((x(2) - x(1)) ./ ( 1+10.^((log10(x(4))-conc_range).* x(3))));

dose_data_mtx = cellfun( dose_curve,data_mtx,repmat({conc_range},size(data_mtx)),'UniformOutput',false);

%% timing of the algorithm
% trial_num = 100;
% combination_gen_time = zeros(trial_num,1);
% solution_found_time  = zeros(trial_num,1);
% 
% for t=1:trial_num
%% generate possible pairs based on the secifications
tic;
all_possible_channels = get_possible_channels(channel_num,num_of_devices,num_of_AHLs);
fprintf('%d possible channel combinations were generated in %f sec\n',size(all_possible_channels,3),toc);
fprintf(fileID,'---------------------------------------\n');
% combination_gen_time(t) = toc;
num_of_possible_combinations = size(all_possible_channels,3);

% randomize the combinations list
new_idx = randperm(size(all_possible_channels,3));
all_possible_channels = all_possible_channels(:,:,new_idx);

%% calculate the ortoganilty region for the current candidates
score = zeros(num_of_possible_combinations,channel_num);
solution_counter = 0;
tic
% loop through all possible channel combinations
for k = 1:num_of_possible_combinations
    % check all communication pairs in a combination
    signal_foldchange_bin = zeros(channel_num,size(dose_data_mtx{1},2));
    for channel=1:channel_num
        % step 1: thresholding for the signal
        % calculate fold change: compared to the lowest response
        dose_mtx_coord    = all_possible_channels(:,channel,k);
        device_candidate  = dose_mtx_coord(1);
        ahl_candidate     = dose_mtx_coord(2);
        signal_foldchange = dose_data_mtx{device_candidate,ahl_candidate};
        signal_foldchange_bin(channel,:)  = (signal_foldchange./signal_foldchange(1)) > act_th;
        
        % step 2: thresholding for crosstalk
        cross_idx = all_possible_channels(1,:,k);
        % remove the current signal from the list of potential cross talks
        cross_idx((cross_idx == device_candidate)) = [];
        % calculate fold change: compared to the lowest response
        bin_data  = cellfun(@(x) (x/x(1))<cross_th,dose_data_mtx(cross_idx,ahl_candidate),'UniformOutput',false);
        local_mtx = cell2mat(bin_data);
        cross_talk_bin = all(local_mtx,1);
        % compute orthogonality range
        score(k,channel) = sum(signal_foldchange_bin(channel,:) & cross_talk_bin);
    end
    if all(score(k,:)>0)
        solution_counter = solution_counter + 1;
        if strcmpi(config.output ,'verbose') || strcmpi(config.show_result,'file')
            fprintf(fileID,['solution: %d \t \t' repmat('%4i', 1, size(signal_foldchange_bin, 2)) '\n'],solution_counter,[min_AHL_concentration:max_AHL_concentration]);
            for channel = 1:channel_num
                d   = devices{all_possible_channels(1,channel,k)};
                ahl = AHLs{all_possible_channels(2,channel,k)};
                fprintf(fileID,['%s - %s \t' repmat('%4i', 1, size(signal_foldchange_bin, 2)) '\n'], d,ahl,signal_foldchange_bin(channel,:)');
            end
            fprintf(fileID,'---------------------------------------\n');
        end
        if strcmpi(config.solution,'first')
            fprintf('First Solution Found!\n')
            break
        end
    end
end
fprintf('%d solution(s) was/were found in %f sec\n',solution_counter,toc)
% solution_found_time(t) = toc;

% end

%% reporting on the timing
% fprintf('------\n channel num: %d\n Combinations were generated in avg %f sec (min: %f max: %f)\n Solution was found in avg %f sec (min: %f max: %f) \n ------\n',...
%     channel_num,mean(combination_gen_time),min(combination_gen_time),max(combination_gen_time),mean(solution_found_time),min(solution_found_time),max(solution_found_time));

%% reporting
show_zero = 1;
idx_num = 1;
if strcmpi(config.output,'minimal')
    for z=1:size(score,1)
        if show_zero && all(score(z,:) ~= 0)
            report_str = '';
            for channel = 1:channel_num
                d   = devices{all_possible_channels(1,channel,z)};
                ahl = AHLs{all_possible_channels(2,channel,z)};
                report_str = [report_str sprintf('%s - %s: %d',d,ahl,score(z,channel)) '\t| ' ];
            end
            fprintf(['%d: ' report_str '\n'],idx_num);
            idx_num = idx_num + 1;
        end
    end
end