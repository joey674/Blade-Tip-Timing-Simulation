

clear;
datasets = {
    %% EO24 SD24
    '../Dataset/EO24/App4_SD24_37000_down_offen_220_run3_RPMRange-34100_34700_thres_S7_RotorMethod_movmean_10.mat', 12, [24],"Set1";
    % '../Dataset/EO24/SD24/App4_SD24_37000_down_offen_220_run2_RPMRange-34200_34700_thres_S7_RotorMethod_movmean_10.mat', 12, [24],"Set2";
    % '../Dataset/EO24/SD24/App4_SD24_37000_down_offen_220_run3_RPMRange-34200_34800_thres_S7_RotorMethod_movmean_10.mat', 12, [24],"Set3";
    % '../Dataset/EO24/SD24/App4_SD24_33000_up_offen_220_run2_RPMRange-34000_35000_thres_S6_RotorMethod_movmean_10.mat', 12, [24],"Set4";
    % '../Dataset/EO24/SD24/App4_SD24_33000_up_offen_220_run3_RPMRange-34000_35000_thres_S6_RotorMethod_movmean_10.mat', 12, [24],"Set5";
    %% EO8/20
    % '../Dataset/EO08_20/App4_DRH30_29300_down_offen_180_run1_RPMRange-28200_29200_thres_S8_RotorMethod_movmean_10.mat', 12, [20, 8],"Set1";
    % '../Dataset/EO08_20/App4_DRH30_29300_down_offen_150_run2_RPMRange-28200_29200_thres_S8_RotorMethod_movmean_10.mat', 12, [20, 8],"Set2";
    % '../Dataset/EO08_20/App4_DRH30_29300_down_offen_150_run3_RPMRange-28200_29200_thres_S8_RotorMethod_movmean_10.mat', 12, [20, 8],"Set3";
    % '../Dataset/EO08_20/App4_DRH30_29300_down_offen_150_run4_RPMRange-28200_29200_thres_S8_RotorMethod_movmean_10.mat', 12, [20, 8],"Set4";
    %% 
    % '../Dataset/EO16/App4_SD24_37000_down_offen_220_run3_RPMRange-35600_36200_thres_S7_RotorMethod_movmean_10.mat', 12, [16],"Set1";
    % '../Dataset/EO17/App4_SD24_37000_down_offen_220_run3_RPMRange-33600_34200_thres_S7_RotorMethod_movmean_10.mat', 12, [17],"Set1";
    % '../Dataset/EO30/App4_DRH30_29300_down_offen_180_run1_RPMRange-27200_27700_thres_S8_RotorMethod_movmean_10.mat', 12, [30],"Set1";
};

for k = 1:size(datasets, 1)

    dataset_path = datasets{k, 1};
    n_blades = datasets{k, 2};
    EO = datasets{k, 3};
    Tag = datasets{k, 4};
    data = load(dataset_path, 'mean_RPM', 'P_Magn','P_Phase', 'Fit_Error');
    blade_sets = Preprocess(EO, data, n_blades);    
    for i = 1:length(blade_sets)
        Damping_MutiDegreeOfFreedom(blade_sets{i}, EO(i),Tag);
        % Damping_SingleDegreeOfFreedom(blade_sets{i}, EO(i),Tag);
        % Damping_HalfPowerBandWidth(blade_sets{i}, EO(i),Tag);

        % Research_plot(blade_sets{i}, EO(i));
        % Research_FindPeakAutomatic(blade_sets{i}, EO(i));
    end
end



