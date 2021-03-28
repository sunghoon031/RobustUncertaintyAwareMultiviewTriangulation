# RobustUncertaintyAwareMultiviewTriangulation
MATLAB implementation of our multiview triangulation method proposed in "Robust Uncertainty-Aware Multiview Triangulation" ([arXiv](https://arxiv.org/abs/2008.01258)).

#### Instruction to run the triangulation code: 

1. Download the 3D uncertainty grid (learned from the simulations): [uncertainty.mat](https://github.com/sunghoon031/RobustUncertaintyAwareMultiviewTriangulation/blob/master/uncertainty.mat).
2. In the same folder, download the [script](https://github.com/sunghoon031/RobustUncertaintyAwareMultiviewTriangulation/blob/master/RobustUncertaintyAwareMultiviewTriangulation_ReleaseCode.m).
3. Run the script on Matlab.

#### Instruction to run the uncertainty generation code:

1. Open GenerateUncertaintyData.m
2. Set `load_precomputed_results = false;`
3. Set `n_simulations` to the number of simulations you want to run for each number of cameras in `n_cameras`.
4. Run. If you use the default `n_simulations`, it will take a very long time to finish (possibly, more than a few days)...
5. This will create two files, `uncertainty.mat` and `results_inliers.mat`. **Caution: this will overwrite the original `uncertainty.mat`! Make sure you rename the filename if you want to avoid this!**
6. Set `load_precomputed_results = true;` and run. Then it will plot some of the results.
