clear;clc;
control_data_list_pre = dir('D:\youyi_fucha\baseline_data\Control\xcp_abcd\sub-*\func\sub-*desc-residual_den-91k_bold.dtseries.nii');
control_data_list_now = dir('D:\youyi_fucha\baseline_data\Control\sub-*-91k_desc-denoised_bold.dtseries.nii');
control_alldata = [control_data_list_pre;control_data_list_now];

PA_data_list_pre = dir('D:\youyi_fucha\baseline_data\Patient\xcp_abcd\sub-*\func\sub-*desc-residual_den-91k_bold.dtseries.nii');
PA_data_list_now = dir('D:\youyi_fucha\baseline_data\Patient\sub-*-91k_desc-denoised_bold.dtseries.nii');
PA_alldata = [PA_data_list_pre;PA_data_list_now];

dlabel_tmp = cifti_read('D:\language_template\code\CBIG\stable_projects\brain_parcellation\Schaefer2018_LocalGlobal\Parcellations\HCP\fslr32k\cifti\Schaefer2018_1000Parcels_Kong2022_17Networks_order.dlabel.nii');
label_unique = unique(dlabel_tmp.cdata);
lack_node_id = setdiff(1:1000,label_unique);
network_label = readtable('D:\youyi_fucha\code\Schaefer_1000_netlabel.csv');
for sub=1:numel(PA_alldata)
    filedata = [PA_alldata(sub).folder filesep PA_alldata(sub).name];
    datasub = extract_cortex(cifti_read(filedata));
    schaefer_1000 = compute_mat_base_label(datasub,label_unique,dlabel_tmp.cdata);
    schaefer_1000 = schaefer_1000(2:end,:);
    data_FC = FisherTransform(corr(schaefer_1000'));
    [S, W, B]  = segregation(data_FC,network_label.x7_network_label);
    [S_all, S_same, S_other, W_same, B_all, B_same, B_other]=segregation_by_type_eqcont(data_FC,network_label.x17_network_label,network_label.x7_network_label,'diagzero','negzero');%,'negzero'
    sub_PA_first(sub,:,1) = S_all;
    sub_PA_first(sub,:,2) = S_same;
    sub_PA_first(sub,:,3) = S_other;
    sub_PA_first(sub,:,4) = W_same;
    sub_PA_first(sub,:,5) = B_all;
    sub_PA_first(sub,:,6) = B_same;
    sub_PA_first(sub,:,7) = B_other;
    sub_PA_first(sub,:,8) = S;
    sub_PA_first(sub,:,9) = W;
    sub_PA_first(sub,:,10) = B;
    sub_PA_first_list{sub}=PA_alldata(sub).name;
    disp(num2str(sub));
end


for sub=1:numel(control_alldata)
    filedata = [control_alldata(sub).folder filesep control_alldata(sub).name];
    datasub = extract_cortex(cifti_read(filedata));
    schaefer_1000 = compute_mat_base_label(datasub,label_unique,dlabel_tmp.cdata);
    schaefer_1000 = schaefer_1000(2:end,:);
    data_FC = FisherTransform(corr(schaefer_1000'));
    [S, W, B]  = segregation(data_FC,network_label.x17_network_label);
    [S_all, S_same, S_other, W_same, B_all, B_same, B_other]=segregation_by_type_eqcont(data_FC,network_label.x17_network_label,network_label.x7_network_label,'diagzero','negzero');%,'negzero'
    sub_con_first(sub,:,1) = S_all;
    sub_con_first(sub,:,2) = S_same;
    sub_con_first(sub,:,3) = S_other;
    sub_con_first(sub,:,4) = W_same;
    sub_con_first(sub,:,5) = B_all;
    sub_con_first(sub,:,6) = B_same;
    sub_con_first(sub,:,7) = B_other;
    sub_con_first(sub,:,8) = S;
    sub_con_first(sub,:,9) = W;
    sub_con_first(sub,:,10) = B;
    sub_con_first_list{sub}=control_alldata(sub).name;
    disp(num2str(sub));
end
sub_con_first_list = sub_con_first_list';
sub_PA_first_list = sub_PA_first_list';
sub_con_first(:,:,11) = 0;
sub_PA_first(:,:,11) = 1;
%%
allsub_data = [sub_con_first;sub_PA_first];
Network_label = {'Deafualt','Language','Control','VAN','DAN','SMN','VIS'};
Metrics_name ={'S_all', 'S_same', 'S_other', 'W_same', 'B_all', 'B_same', 'B_other','S','W','B'};
alldata_use = [];
alldata_use_name = {};
iter = 1;
for metric=1:10
    if metric<8
        for n=1:7
           group_info = allsub_data(:,n,11);
           data_metric = allsub_data(:,n,metric);
           alldata_use = [alldata_use,data_metric];
           alldata_use_name{1,iter} = [Network_label{n} '_' Metrics_name{metric}];
           iter = iter+1;
           if any(isnan(data_metric))
               continue;
           else
               [h,p,ci,stats] = ttest2(data_metric(1:71),data_metric(72:end));
               all_result(n,metric,1) = p;
               all_result(n,metric,2) = stats.tstat;
               if p < 0.01
                  disp([Network_label{n} ':  ' Metrics_name{metric} ' has sig diff']);  
               end 
           end
        end
    else
       group_info = allsub_data(:,n,11);
       data_metric = allsub_data(:,1,metric);
       alldata_use = [alldata_use,data_metric];
       alldata_use_name{1,iter} = ['WholeBrain_' Metrics_name{metric}];
       iter = iter+1;
       if any(isnan(data_metric))
           continue;
       else
           [h,p,ci,stats] = ttest2(data_metric(1:71),data_metric(72:end));
           all_result(n,metric,1) = p;
           all_result(n,metric,2) = stats.tstat;
           if p < 0.01
              disp([Network_label{n} ':  ' Metrics_name{metric} ' has sig diff']);  
           end 
       end
    end
end

all_metric_data = array2table(alldata_use);
all_metric_data.Properties.VariableNames = alldata_use_name;
all_metric_data.subname = [sub_con_first_list;sub_PA_first_list];
all_metric_data.group = group_info;
writetable(all_metric_data,'all_metric_data_sub.csv')

