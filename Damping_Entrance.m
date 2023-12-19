% % dataset EO24
% dataset_path = 'App4_SD24_37000_down_offen_220_run3_RPMRange-34100_34700_thres_S7_RotorMethod_movmean_10.mat';
% n_blades = 12;
% EO = [24];
% smoothness = 100;

%% dataset EO8/20
dataset_path = 'App4_DRH30_29300_down_offen_180_run1_RPMRange-28200_29200_thres_S8_RotorMethod_movmean_10.mat';
n_blades = 12;
EO = [8,20];
smoothness = 300;


%% start
data = load(dataset_path, 'mean_RPM' , 'P_Magn', 'Fit_Error');
blade = Damping_PreProcess(EO,data,n_blades);
for i = 1:length(blade)
    n_modes = MDOF_FindNumOfModes(blade{i},smoothness);
    % Damping_MutiDegreeOfFreedom(dataset_path,blade{i},n_modes,smoothness,peak_interval_width);
end



