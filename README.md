# Explanations On Simulations

## Main Settings
1. homo_case.R: investigates settings S1 to S5 when homo errors.
2. hetero_case.R: investigates settings S1 to S5 when hetero errors.
3. homo_case_S2.R: investigates varying violation strengths ($\tau$) for S2 when homo errors.
4. hetero_case_S2.R: investigates varying violation strengths ($\tau$) for S2 when hetero errors.

## Tuning Parameters
1. simu_illustrate_LUa.R: illustrates varying initial intervals and grid-size (LUa) effect on Searching and Sampling
for S1 and S2 when homo errors.
2. simu_illustrate_Lua_additional.R: Specifically for Searching method, it adds a few more options on grid-size based on the above
one.
3. simu_illustrate_S0.R: investigates the setting S0 for varying thresholds.
4. simu_illustrate_samplingLambda.R: illustrates the effect of varying proportions on Sampling method for settings S1 to S5
5. simu_illustrate_samplingM.R: illustrates the effect of sampling times M on Sampling method for settings S1 and S2
6. simu_bootstrap.R: compares the bootstrap threshold vs. regular threshold for Searching method for settings S1 to S5
7. simu_filtering.R: compares filtering vs. non-filtering effects on Sampling method for settings S1 to S5

## Others
1. simu_highd_pz100.R: investigates SearchingSampling on high-d S2 setting when $p_z=100$
2. simu_highd_pz200.R: investigates SearchingSampling on high-d S2 setting when $p_z=200$
3. simu_CIIV.R: investigates SearchingSampling on settings CIIV-1 and CIIV-2
4. simu_threshold_S2.R: computes the $max_{T_{j,k}}$ for setting S2
5. homo_S2_hists.R: get point estimators of TSHT and CIIV method for setting S2, prepared for later histograms
