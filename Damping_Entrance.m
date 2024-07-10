

clear;
datasets = {
    '../Dataset/EO24/App4_SD24_37000_down_offen_220_run3_RPMRange-34100_34700_thres_S7_RotorMethod_movmean_10.mat', 12, [24];
    % '../Dataset/EO08_20/App4_DRH30_29300_down_offen_180_run1_RPMRange-28200_29200_thres_S8_RotorMethod_movmean_10.mat', 12, [20, 8];
    % '../Dataset/EO16/App4_SD24_37000_down_offen_220_run3_RPMRange-35600_36200_thres_S7_RotorMethod_movmean_10.mat', 12, [16];
    % '../Dataset/EO17/App4_SD24_37000_down_offen_220_run3_RPMRange-33600_34200_thres_S7_RotorMethod_movmean_10.mat', 12, [17];
    % '../Dataset/EO30/App4_DRH30_29300_down_offen_180_run1_RPMRange-27200_27700_thres_S8_RotorMethod_movmean_10.mat', 12, [30];
};
for k = 1:size(datasets, 1)
    dataset_path = datasets{k, 1};
    n_blades = datasets{k, 2};
    EO = datasets{k, 3};
    data = load(dataset_path, 'mean_RPM', 'P_Magn','P_Phase', 'Fit_Error');
    blade_sets = Preprocess(EO, data, n_blades);
    for i = 1:length(blade_sets)
        Damping_MutiDegreeOfFreedom(blade_sets{i}, EO(i));
        % Damping_SingleDegreeOfFreedom(blade_sets{i}, EO(i));
        % Damping_HalfPowerBandWidth(blade_sets{i}, EO(i));

        % Research_plot(blade_sets{i}, EO(i));
        % Research_FindPeakAutomatic(blade_sets{i}, EO(i));
    end
end



