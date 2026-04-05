clear; clc;
addpath('D:\language_template\code');
%% 1. 数据路径和标签文件读取
fprintf('正在读取数据文件列表...\n');

% 控制组数据文件
control_data_list_pre = dir('D:\youyi_fucha\baseline_data\Control\xcp_abcd\sub-*\func\sub-*desc-residual_den-91k_bold.dtseries.nii');
control_data_list_now = dir('D:\youyi_fucha\baseline_data\Control\sub-*-91k_desc-denoised_bold.dtseries.nii');
control_alldata = [control_data_list_pre; control_data_list_now];

% 患者组数据文件
PA_data_list_pre = dir('D:\youyi_fucha\baseline_data\Patient\xcp_abcd\sub-*\func\sub-*desc-residual_den-91k_bold.dtseries.nii');
PA_data_list_now = dir('D:\youyi_fucha\baseline_data\Patient\sub-*-91k_desc-denoised_bold.dtseries.nii');
PA_alldata = [PA_data_list_pre; PA_data_list_now];

fprintf('找到控制组文件: %d个\n', length(control_alldata));
fprintf('找到患者组文件: %d个\n', length(PA_alldata));

% 读取标签文件
dlabel_tmp = cifti_read('D:\language_template\code\CBIG\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\HCP\fslr32k\cifti\Schaefer2018_1000Parcels_Kong2022_17Networks_order.dlabel.nii');
label_unique = unique(dlabel_tmp.cdata);
lack_node_id = setdiff(1:1000, label_unique);
network_label = readtable('D:\youyi_fucha\code\Schaefer_1000_netlabel.csv');
colour_map = dlabel_tmp.diminfo{2}.maps.table;
for n=1:17
    subnetid_pos = find(network_label.x17_network_label==n)+1;
    colorset(n,:) = colour_map(subnetid_pos(1)).rgba;
    names{n,1} = colour_map(subnetid_pos(1)).name;
end
fprintf('标签文件读取完成，有效节点数: %d\n', length(label_unique)-1); % 减1是因为包含0

%% 2. 处理患者组数据
fprintf('\n开始处理患者组数据...\n');
PA_FC_matrices = [];
PA_subject_names = {};

for sub = 1:numel(PA_alldata)
    try
        filedata = [PA_alldata(sub).folder filesep PA_alldata(sub).name];
        fprintf('处理患者 %d/%d: %s\n', sub, numel(PA_alldata), PA_alldata(sub).name);
        
        % 读取和处理数据
        datasub = extract_cortex(cifti_read(filedata));
        schaefer_1000 = compute_mat_base_label(datasub, label_unique, dlabel_tmp.cdata);
        schaefer_1000 = schaefer_1000(2:end, :); % 去除第一行（背景）
        
        % 计算功能连接矩阵
        data_FC = FisherTransform(corr(schaefer_1000'));
        
        % 存储FC矩阵和被试名称
        if sub == 1
            PA_FC_matrices = zeros(size(data_FC, 1), size(data_FC, 2), numel(PA_alldata));
        end
        PA_FC_matrices(:, :, sub) = data_FC;
        PA_subject_names{sub} = PA_alldata(sub).name;
        
    catch ME
        fprintf('处理患者 %d 时出错: %s\n', sub, ME.message);
        continue;
    end
end

%% 3. 处理控制组数据
fprintf('\n开始处理控制组数据...\n');
Control_FC_matrices = [];
Control_subject_names = {};

for sub = 1:numel(control_alldata)
    try
        filedata = [control_alldata(sub).folder filesep control_alldata(sub).name];
        fprintf('处理控制组 %d/%d: %s\n', sub, numel(control_alldata), control_alldata(sub).name);
        
        % 读取和处理数据
        datasub = extract_cortex(cifti_read(filedata));
        schaefer_1000 = compute_mat_base_label(datasub, label_unique, dlabel_tmp.cdata);
        schaefer_1000 = schaefer_1000(2:end, :); % 去除第一行（背景）
        
        % 计算功能连接矩阵
        data_FC = FisherTransform(corr(schaefer_1000'));
        
        % 存储FC矩阵和被试名称
        if sub == 1
            Control_FC_matrices = zeros(size(data_FC, 1), size(data_FC, 2), numel(control_alldata));
        end
        Control_FC_matrices(:, :, sub) = data_FC;
        Control_subject_names{sub} = control_alldata(sub).name;
        
    catch ME
        fprintf('处理控制组 %d 时出错: %s\n', sub, ME.message);
        continue;
    end
end

%% 4. 计算组水平平均功能连接矩阵
fprintf('\n计算组水平平均功能连接矩阵...\n');
% 计算平均FC矩阵
PA_mean_FC = mean(PA_FC_matrices, 3);
Control_mean_FC = mean(Control_FC_matrices, 3);
sublist_CONNFC.PA_mean_FC = PA_mean_FC;
sublist_CONNFC.Control_mean_FC = Control_mean_FC;
save('D:\youyi_fucha\network_segeragation_first_test\sublist_CONNFC.mat','sublist_CONNFC');% 计算组间差异矩阵
Group_diff_FC = PA_mean_FC - Control_mean_FC;
fprintf('患者组平均FC矩阵大小: %dx%d\n', size(PA_mean_FC));
fprintf('控制组平均FC矩阵大小: %dx%d\n', size(Control_mean_FC));
%% 5. 准备网络标签和分隔线位置
fprintf('\n准备网络标签和分隔线位置...\n');

% 获取网络标签（去除背景标签）
network_17_labels = network_label.x17_network_label(1:end); % 去除背景标签
network_7_labels = network_label.x7_network_label(1:end);   % 去除背景标签

% 检查数据结构，确定左右脑分界点
left_brain_end = 1000 / 2; % 假设前一半是左脑，后一半是右脑

fprintf('总节点数: %d, 左脑节点: 1-%d, 右脑节点: %d-%d\n', ...
        total_nodes, left_brain_end, left_brain_end+1, total_nodes);

% 分别处理左脑和右脑的网络边界
% 左脑17网络边界
left_17_labels = network_17_labels(1:left_brain_end);
unique_left_17 = unique(left_17_labels);
left_17_boundaries = [];
for i = 1:length(unique_left_17)
    network_positions = find(left_17_labels == unique_left_17(i));
    if i < length(unique_left_17)
        left_17_boundaries = [left_17_boundaries, max(network_positions) + 0.5];
    end
end

% 右脑17网络边界
right_17_labels = network_17_labels(left_brain_end+1:end);
unique_right_17 = unique(right_17_labels);
right_17_boundaries = [];
for i = 1:length(unique_right_17)
    network_positions = find(right_17_labels == unique_right_17(i));
    if i < length(unique_right_17)
        right_17_boundaries = [right_17_boundaries, max(network_positions) + left_brain_end + 0.5];
    end
end

% 合并左右脑17网络边界
network_17_boundaries = [left_17_boundaries, right_17_boundaries];

% 左脑7网络边界
left_7_labels = network_7_labels(1:left_brain_end);
unique_left_7 = unique(left_7_labels);
left_7_boundaries = [];
for i = 1:length(unique_left_7)
    network_positions = find(left_7_labels == unique_left_7(i));
    if i <= length(unique_left_7)
        left_7_boundaries = [left_7_boundaries, max(network_positions) + 0.5];
    end
end

% 右脑7网络边界
right_7_labels = network_7_labels(left_brain_end+1:end);
unique_right_7 = unique(right_7_labels);
right_7_boundaries = [];
for i = 1:length(unique_right_7)
    network_positions = find(right_7_labels == unique_right_7(i));
    if i <= length(unique_right_7)
        right_7_boundaries = [right_7_boundaries, max(network_positions) + left_brain_end + 0.5];
    end
end

% 合并左右脑7网络边界
network_7_boundaries = [left_7_boundaries, right_7_boundaries];

% 添加左右脑分界线
brain_boundary = left_brain_end + 0.5;

fprintf('左脑17网络边界: %d个\n', length(left_17_boundaries));
fprintf('右脑17网络边界: %d个\n', length(right_17_boundaries));
fprintf('左脑7网络边界: %d个\n', length(left_7_boundaries));
fprintf('右脑7网络边界: %d个\n', length(right_7_boundaries));
fprintf('左右脑分界线位置: %.1f\n', brain_boundary);

%% 6. 可视化组水平网络
fprintf('\n开始可视化组水平网络...\n');

% 创建患者组功能连接图


figure('Position', [100, 100, 800, 700]);
imagesc(PA_mean_FC);
colorcode = 103;
colormap(slanCM(colorcode));
colorbar;
title(sprintf('患者组平均功能连接矩阵 (n=%d)', size(PA_FC_matrices, 3)), ...
      'FontSize', 14, 'FontWeight', 'bold');
xlabel('脑区节点', 'FontSize', 12);
ylabel('脑区节点', 'FontSize', 12);
axis equal tight;
caxis([-0.5, 0.5]);
% 添加网络分隔线 - 患者组
hold on;
box on;
network_7_boundaries = [0,network_7_boundaries];
% 7大网络分隔线（黑色粗线）
for i = 1:length(network_7_boundaries)
    line([network_7_boundaries(i), network_7_boundaries(i)], [0.5, size(PA_mean_FC, 1)+0.5], ...
         'Color', 'black', 'LineWidth', 0.4);
    line([0.5, size(PA_mean_FC, 2)+0.5], [network_7_boundaries(i), network_7_boundaries(i)], ...
         'Color', 'black', 'LineWidth', 0.4);
end
hold off;
set(gca, 'XColor', 'none', 'YColor', 'none');

% 保存患者组图片
saveas(gcf, 'Patient_Group_FC_Network.png', 'png');
saveas(gcf, 'Patient_Group_FC_Network.fig', 'fig');

figure('Position', [100, 100, 800, 700]);
imagesc(Control_mean_FC);
colorcode = 103;
colormap(slanCM(colorcode));
colorbar;
title(sprintf('健康组平均功能连接矩阵 (n=%d)', size(PA_FC_matrices, 3)), ...
      'FontSize', 14, 'FontWeight', 'bold');
xlabel('脑区节点', 'FontSize', 12);
ylabel('脑区节点', 'FontSize', 12);
axis equal tight;
caxis([-0.5, 0.5]);
% 添加网络分隔线 - 患者组
hold on;
box on;
network_7_boundaries = [0,network_7_boundaries];
% 7大网络分隔线（黑色粗线）
for i = 1:length(network_7_boundaries)
    line([network_7_boundaries(i), network_7_boundaries(i)], [0.5, size(PA_mean_FC, 1)+0.5], ...
         'Color', 'black', 'LineWidth', 0.4);
    line([0.5, size(PA_mean_FC, 2)+0.5], [network_7_boundaries(i), network_7_boundaries(i)], ...
         'Color', 'black', 'LineWidth', 0.4);
end
hold off;
set(gca, 'XColor', 'none', 'YColor', 'none');

