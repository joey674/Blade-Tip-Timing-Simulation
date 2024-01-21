clear
%% dataset EO24
% dataset_path = 'EO24/App4_SD24_37000_down_offen_220_run3_RPMRange-34100_34700_thres_S7_RotorMethod_movmean_10.mat';
% n_blades = 12;
% EO = [24];  

%% dataset EO08/20
dataset_path = 'EO08_20/App4_DRH30_29300_down_offen_180_run1_RPMRange-28200_29200_thres_S8_RotorMethod_movmean_10.mat';
n_blades = 12;
EO = [8,20];

%% dataset EO16
% dataset_path = 'EO16/App4_SD24_37000_down_offen_220_run3_RPMRange-35600_36200_thres_S7_RotorMethod_movmean_10.mat';
% n_blades = 12;
% EO = [16];

%% dataset EO17
% dataset_path = 'EO17/App4_SD24_37000_down_offen_220_run3_RPMRange-33600_34200_thres_S7_RotorMethod_movmean_10.mat';
% n_blades = 12;
% EO = [17];

%% dataset EO30
% dataset_path = 'EO30/App4_DRH30_29300_down_offen_180_run1_RPMRange-27200_27700_thres_S8_RotorMethod_movmean_10.mat';
% n_blades = 12;
% EO = [30];


%% start
data = load(dataset_path, 'mean_RPM' , 'P_Magn', 'Fit_Error');
blade_sets = Damping_PreProcess(EO,data,n_blades);
for i = 1:length(blade_sets)
    Damping_MutiDegreeOfFreedom(dataset_path,blade_sets{i},EO(i));
end



