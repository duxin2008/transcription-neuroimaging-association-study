clc, clear, close all
%% select the statistical map
[stat_name, stat_path] = uigetfile('*.nii;*.img','please select the statistical map, any resolution');
stat_file = [stat_path stat_name];
%% select the mask
[mask_name, mask_path] = uigetfile('*.nii;*.img','please select a mask to confine samples, same resolution as statistical map');
mask_file = [mask_path mask_name];
%% select the sample coordinate file
[sample_name, sample_path] = uigetfile('*.mat','please select the sample coordinate file, mat format');
sample_file = [sample_path sample_name];
temp = load(sample_file);
sample_XYZ = temp.SampleCoordinates(:,2:4);
%% sphere radius
cc = inputdlg('How many times of the voxel size (sphere radius)');
radius = str2num(cc{1});
%% read mask coordinate
mask_V = spm_vol(mask_file);
[mask_data, mask_XYZ] = spm_read_vols(mask_V);
mask_index = find(mask_data);
mask_data_use = mask_data(mask_index);
mask_XYZ = mask_XYZ(:,mask_index);
%% read the statistical map
stat_V = spm_vol(stat_file);
stat_data = spm_read_vols(stat_V);
%% compute distance
Distance = pdist2(sample_XYZ,mask_XYZ');
M_distance = min(Distance,[],2);
%% distance threshold (1 voxel), find samples in mask
sample_in_mask_index = find(M_distance<=abs(stat_V.mat(1)));
sample_in_mask_coord = sample_XYZ(sample_in_mask_index,:);
sample_info = [sample_in_mask_index,sample_in_mask_coord];
%% extract the sphere signal of statistical map based on sample site
sample_wise_statistic = zeros(size(sample_in_mask_coord,1),1);
for i = 1:size(sample_in_mask_coord,1)
    disp(['Total sample:' num2str(length(sample_in_mask_index))  '; Now: sample ' num2str(i)]);
    [SphereData,Header] = y_Sphere(sample_in_mask_coord(i,:), radius*abs(stat_V.mat(1)), stat_file);
    Sphere = SphereData.*mask_data;
    sample_wise_statistic(i,1) = mean(stat_data(find(Sphere)));
end
SampleGeneExpression = temp.SampleGeneExpression(sample_in_mask_index,:);
SampleGeneExpression(:,1) = [];
[r,p] = corr(sample_wise_statistic,SampleGeneExpression);
FDR = mafdr(p,'BHFDR',true);
xlswrite('F:\DUXIN\KZBY\Express_Allen\results\wholebrain_Rhippo_sample_wise_statistic_test.xlsx',sample_wise_statistic,'sample_wise_statistic');
xlswrite('F:\DUXIN\KZBY\Express_Allen\results\wholebrain_Rhippo_SampleGeneExpression_test.xlsx',SampleGeneExpression,'SampleGeneExpression');

%xlswrite('F:\DUXIN\KZBY\Express_Allen\results\test.xlsx',r','r');
%xlswrite('F:\DUXIN\KZBY\Express_Allen\results\test.xlsx',p','p');
%xlswrite('F:\DUXIN\KZBY\Express_Allen\results\test.xlsx',FDR','fdr');



