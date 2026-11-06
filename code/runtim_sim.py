import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



import numpy as np
import pandas as pd
from debiasm.torch_functions import rescale

import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
import torch
import os

from scipy.stats import wilcoxon
import time

import os
import numpy as np
import pandas as pd
from debiasm.sklearn_functions import OnlineDebiasMClassifier, \
                                            DebiasMClassifier
from debiasm.torch_functions import rescale
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, cross_val_predict, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score
import subprocess

sns.set(rc={'axes.facecolor':'white', 
            'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
            'grid.color': 'black'
            }, 
        
       font_scale=2)

sns.set_style('ticks')

import torch.nn.functional as F
def c_scrub_transform(self, x, p):
    batch_inds, x = x[:, 0], x[:, 1:]

    x = F.normalize( torch.pow(2, self.batch_weights[batch_inds.long()] ) * x, p=1 )


    x -= ( torch.cat([ self.C_SCRuB_Layers[str(batch)].gamma.unsqueeze(0)
                                            for batch in batch_inds.long() ], 
                                dim=0 ) * ( 1 - p ) )

    x = F.normalize( F.relu(x), p=1 )
    return(x)

import sys
sys.path.append('../../Clean_analysis/NEWER_DM_C_SCRuB/c_scrub/')
from sklearn_functions import C_SCRuB_DebiasMClassifier
sys.path.append('../../Clean_analysis/C-SCRuB/')
from c_scrub import C_SCRuB


def set_up_prediction_eval_run(samples_,
                               phenotype_std,
                               read_depth, 
                               contam_levels,
                               n_studies=4, 
                               n_per_study=88,
                               bcf_stdv = 2, 
                               linear_weights_stdv = 2,
                               seed=0):
    
    np.random.seed(seed)
    torch.manual_seed(seed)
    
    n_features=samples_.shape[1]
    
    
    samples=samples_.copy()
    
    
    ## simulated bcfs
    true_bcfs = np.random.normal(scale=bcf_stdv,
                                 size = (n_studies, n_features)
                                 )
    
    ## simulated linear weights
    true_weights = np.random.normal(scale=linear_weights_stdv, 
                                    size = n_features
                                    )
    
    ## enforce sparsity in relevant weights
    true_weights = true_weights

    ## simulated phenotypes / covariates
    phenotypes = samples @ true_weights[:, np.newaxis]
    if phenotype_std > 0:
        phenotypes += np.random.normal(scale = phenotype_std,
                                       size=phenotypes.shape)
    
    covar_weights = np.random.normal(scale=linear_weights_stdv, 
                                     size = n_features
                                     )
    covariates = samples @ true_weights[:, np.newaxis]
    if phenotype_std > 0:
        covariates += np.random.normal(scale = phenotype_std,
                                       size=phenotypes.shape)
    
    
    studies = pd.Series( [ i // n_per_study for i in range(samples.shape[0]) 
                                     ] , 
                        index=samples.index, 
                        name='studies')
    
    
    contaminants = rescale( np.random.lognormal(size=(n_studies, n_features), sigma = 1.5) )
    contaminants[contaminants<1e-2]=0
    contaminants=rescale(contaminants)
    
    
    
    contaminated_samples = ( samples.values * (1-contam_levels) ) + \
                           ( contaminants[studies] * contam_levels )

    samples_transformed = pd.DataFrame( ( rescale( 
                            contaminated_samples *\
                                    np.power(2, true_bcfs)[studies] 
                                    )*read_depth
                              ).astype(int), 
                                       index=samples.index, 
                                       columns=samples.columns
                                      )
    
    y=phenotypes.values[:, 0] > phenotypes.median().values[0] 
    cov=covariates.values[:, 0] > covariates.median().values[0]
    
    run_description = '_'.join( [a.replace('_', '-') + '-' + str(b)
                                 for a,b in 
                                   {'N_studies':n_studies, 
                                   'N_per_study':n_per_study, 
                                   'N_features':n_features, 
                                   'Pheno_noise':phenotype_std, 
                                   'Read_depth':int(read_depth), 
                                   'BCF_std':bcf_stdv,
                                   'run_num':seed}.items()]
                                )
    
    return(samples_transformed, 
           y, 
           cov,
           studies, 
           samples, 
           contam_levels, 
           contaminants
           )

def main()


    all_times = []
    all_n_per_batch = []

    test_grp=3

    for bcf_stdv in [2.5]:
        for n_per_batch in [24, 48, 72, 96]:
            for seed in range(1,11):
                print(seed)

                np.random.seed(seed)

                samples_true = pd.read_csv('../results/multibatch_simulation_data/round_{}/true_samples.csv'.format(seed),
                                              index_col=0)

                samples_true = samples_true.sample(frac=1)


                conc = pd.read_csv('../results/multibatch_simulation_data/round_{}/observed_concentrations.csv'.format(seed), 
                                   index_col=0
                                   ).head(n_per_batch*4).values

                contam_levels = pd.read_csv('../results/multibatch_simulation_data/round_{}/underlying_contamination_levels.csv'.format(seed), 
                                            index_col=0).values


                np.random.seed(seed)
                torch.manual_seed(seed)
                print(seed)

                samples_transformed,  \
                           y,  \
                           cov, \
                           studies, \
                        samples_true, \
                        contam_levels, \
                        contaminations \
                    = set_up_prediction_eval_run(samples_true.head(n_per_batch*4),
                                     phenotype_std=1,
                                     read_depth=1e4, 
                                     bcf_stdv = bcf_stdv, 
                                     seed=seed, 
                                     contam_levels=contam_levels[:n_per_batch*4], 
                                     n_studies=4, 
                                     n_per_study=n_per_batch
                                     )

                samples_transformed_with_batch = pd.concat([studies, 
                                                            samples_transformed
                                                           ], 
                                                           axis=1 ).values





                test_grp=3


                test_inds = studies.values==test_grp

                start = time.time()

                cs_dm = C_SCRuB_DebiasMClassifier(x_val=samples_transformed_with_batch[test_inds].astype(float), 
                                                  mu_val=conc[test_inds][:, 0],
                                                 )

                cs_dm = cs_dm.fit(samples_transformed_with_batch[~test_inds].astype(float), 
                                  y[~test_inds].astype(float),
                                  mu=conc[~test_inds][:, 0], 
                                  )


                self=cs_dm.model

                p=torch.cat([ self.C_SCRuB_Layers[str(batch)].p for batch in self.unique_bs]).unsqueeze(-1)\
                                        [torch.cat([ self.C_SCRuB_Layers[str(batch)].overall_dataset_inds for batch in self.unique_bs ] 
                                                  )]


                runtime = time.time()-start


                all_times.append(runtime)
                all_n_per_batch.append(n_per_batch)


    plot_df = pd.DataFrame({'Time (s)':all_times,
                            'N':np.array(all_n_per_batch)*4
                            })

    for hide_axes in [False, True]:

        out_path = '../results/plots/DMC_SCRuB_Runtime.pdf'

        plt.figure(figsize=(8,8))
        ax = sns.boxplot(x='N', 
                         y='Time (s)', 
                         data=plot_df, 
                         fliersize=0, 
                         color='#FF6065'
                         )

        ax = sns.swarmplot(x='N', 
                           y='Time (s)', 
                           data=plot_df, 
                           ax=ax, 
                           color='black', 
                           )

        plt.yticks([0, 15, 30, 45])

        if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.title(None)
            ax.set(yticks=[0, 15, 30, 45],
                   yticklabels=[], 
                   xticklabels=[]
                   )
            plt.legend().remove()
            out_path=out_path[:-4] + '-no-axes'+out_path[-4:]


        plt.savefig(out_path, 
                    format='pdf',
                    bbox_inches='tight', 
                    dpi=900)

    
if __name__=='__main__':
    main()
