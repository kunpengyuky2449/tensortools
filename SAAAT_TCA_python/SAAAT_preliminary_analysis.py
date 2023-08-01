# %% import modules
import numpy as np
import tensortools as tt
import matplotlib.pyplot as plt
import mat73

# %% asd
sss
# %% load NTK tensor
# ... specify a numpy array holding the tensor you wish to fit
data = mat73.loadmat(
    "D:\BCM_projects\TCA_project2\SAAAT_preprocess_matlab\SAAAT_tensor_NTK_230801_13_34_06.mat")
# Access the tensor by data["SAAAT_tensor"]

# %% Fit with Non-negative TCA with rank 1-9
# Fit an ensemble of models, 4 random replicates / optimization runs per model rank
ensemble = tt.Ensemble(fit_method="ncp_hals")
ensemble.fit(data["SAAAT_tensor"], ranks=range(1, 2), replicates=2)

# %%
a = ensemble.factors(1)[0]
#%%
fig, axes = plt.subplots(1, 2)
tt.plot_objective(ensemble, ax=axes[0])  # plot reconstruction error as a function of num components.
tt.plot_similarity(ensemble, ax=axes[1])  # plot model similarity as a function of num components.
fig.tight_layout()

# Plot the low-d factors for an example model, e.g. rank-2, first optimization run / replicate.
num_components = 1
replicate = 0
tt.plot_factors(ensemble.factors(num_components)[replicate])  # plot the low-d factors

plt.show()

rank = 1
error = ensemble.objectives(1)
simi = ensemble.similarities(1)[1:]
#%% save the ensemble for later access by python
import pickle
from datetime import datetime
now = datetime.now()
results_dir = '../results/'
time_str = now.strftime("%Y-%m-%d-%H%M%S")
with open(results_dir+'rank1_2rep' + time_str, 'wb') as f:
    # indent=2 is not needed but makes the file human-readable
    # if the data is nested
    pickle.dump(ensemble, f)

#%% save important variables for matlab
rank = 1
error = ensemble.objectives(1)
simi = ensemble.similarities(1)[1:]
factors = ensemble.factors(1)[0].factors
results_dir_mat = '../SAAAT_preprocess_matlab/'
import scipy.io
scipy.io.savemat(results_dir_mat+'rank1'+time_str+'.mat',
                 dict(rank=rank,
                      error=error,
                      simi=simi,
                      factors=factors))

