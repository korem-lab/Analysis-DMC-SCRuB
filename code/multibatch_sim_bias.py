import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


from debiasm.torch_functions import rescale

import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
import torch
import os

from scipy.stats import wilcoxon

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

palette = {'Raw':'#1f77b4',
           'Decontam':'#12A380',
           'Decontam (Low Biomass)':'darkgreen',
           'Decontam Combined':'#078A69',
           'DEBIAS-M':'lightblue',
           'SCRuB':'#d62728',
           'B2':'white',
           'B1':'white',
           'Decontam --> DEBIAS-M':'#12A380',
           'DEBIAS-M --> Decontam':'#12A380',
           'DEBIAS-M --> C-SCRuB':'#e377c2', 
           'DM-C-SCRuB':'#FF6065'}

import sys
sys.path.append('../../Clean_analysis/DM_C_SCRuB/c_scrub/')
from sklearn_functions import C_SCRuB_DebiasMClassifier
sys.path.append('../../Clean_analysis/C-SCRuB/')
from c_scrub import C_SCRuB

import torch.nn.functional as F
def c_scrub_transform(self, x, p):
    batch_inds, x = x[:, 0], x[:, 1:]

    x = F.normalize( torch.pow(2, self.batch_weights[batch_inds.long()] ) * x, p=1 )


    x -= ( torch.cat([ self.C_SCRuB_Layers[str(batch)].gamma.unsqueeze(0)
                                            for batch in batch_inds.long() ], 
                                dim=0 ) * ( 1 - p ) )

    x = F.normalize( F.relu(x), p=1 )
    return(x)

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





def main():
    all_run_descs   = []
    all_debias_jsds = []
    all_csdmcr_jsds   = []
    all_sep_jsds    = []
    all_dec_debias_jsds = []
    all_debias_dec_jsds = []
    all_just_dec_jsds = []
    all_raw_jsds = []

    test_grp=3

    for bcf_stdv in [3, 2.5, 2][::-1]:
        for seed in range(1,11):
            print(seed)

            np.random.seed(seed)

            samples_true = pd.read_csv('../results/multibatch_simulation_data/round_{}/true_samples.csv'.format(seed),
                                          index_col=0)

            samples_true = samples_true.sample(frac=1)

            conc = pd.read_csv('../results/multibatch_simulation_data/round_{}/observed_concentrations.csv'.format(seed),
                               index_col=0).values

            contam_levels = pd.read_csv('../results/multibatch_simulation_data/round_{}/underlying_contamination_levels.csv'.format(seed), 
                                        index_col=0).values

            all_rocs=[]
            all_names=[]
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
                = set_up_prediction_eval_run(samples_true,
                                 phenotype_std=1,
                                 read_depth=1e4, 
                                 bcf_stdv = bcf_stdv, 
                                 seed=seed, 
                                 contam_levels=contam_levels, 
                                 n_studies=4, 
                                 n_per_study=96
                                 )

            samples_transformed_with_batch = pd.concat([studies, 
                                                        samples_transformed
                                                       ], 
                                                       axis=1 ).values


            test_grp=3

            test_inds = studies.values==test_grp

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

            # test_inds = studies.values==test_grp
            dm = DebiasMClassifier(x_val=samples_transformed_with_batch[test_inds].astype(float))

            dm = dm.fit(samples_transformed_with_batch[~test_inds].astype(float), 
                        y[~test_inds].astype(float)
                       )


            dm_jsd=pd.Series( jensenshannon( samples_true.values[test_inds].T, 
                                       dm.transform(samples_transformed_with_batch.astype(float) )[test_inds].T
                                     )**2
                                 )


            csdmcr_jsd = pd.Series( 
                jensenshannon( samples_true.values[test_inds].T, 
                                   c_scrub_transform(cs_dm.model, 
                                      torch.Tensor( samples_transformed_with_batch[test_inds].astype(float) ), 
                                      p[test_inds]
                                      ).detach().numpy().T
                                 )**2\
                            )



            debias_transformed = dm.transform(samples_transformed_with_batch[test_inds].astype(float) )

            cscr = C_SCRuB()
            separate_out = cscr.fit_transform(debias_transformed, 
                                              conc[:, 0][test_inds]
                                              )

            dm_to_cscrub_jsd = pd.Series( jensenshannon( samples_true.values[test_inds].T, 
                                               separate_out.T
                                             )**2
                                         )

            print(dm_jsd.median())
            print(csdmcr_jsd.median())
            print(dm_to_cscrub_jsd.median())



            all_run_descs.append({'bias':bcf_stdv, 
                                  'seed':seed
                                  })

            all_debias_jsds.append(dm_jsd)
            all_csdmcr_jsds.append(csdmcr_jsd)
            all_sep_jsds.append(dm_to_cscrub_jsd)   


            pd.DataFrame({'bt':studies, 
                          'mu':conc[:, 0]}).to_csv('tmp-dataset/md.csv')

            pd.DataFrame({'mu':conc[:, 0][test_inds]}).to_csv('tmp-dataset/md_dm_test.csv')

            pd.DataFrame( rescale( samples_transformed_with_batch[:, 1:] ) ).to_csv('tmp-dataset/data.csv')
            pd.DataFrame( debias_transformed ).to_csv('tmp-dataset/data_dm_test.csv')


            ### run the R script
            os.system('conda run -p /Users/george/opt/anaconda3/envs/r43 Rscript run_simulation_decontam.R')


            dec_deb_out = pd.read_csv('tmp-dataset/decontam_deb_out.csv', index_col=0)


            dm_to_decontam_jsd = pd.Series( jensenshannon( samples_true.values[test_inds].T, 
                                                                               dec_deb_out.T
                                                                             )**2
                                                                         )


            dec_raw_out = samples_transformed_with_batch.copy().astype(float)
            dec_raw_out[:, 1:] = pd.read_csv('tmp-dataset/decontam_raw_out.csv', index_col=0).values


            dm = DebiasMClassifier(x_val=dec_raw_out[test_inds].astype(float))

            dm = dm.fit(dec_raw_out[~test_inds].astype(float), 
                        y[~test_inds].astype(float)
                        )

            dec_to_dm_jsd=pd.Series( jensenshannon( samples_true.values[test_inds].T, 
                                     dm.transform(dec_raw_out.astype(float) )[test_inds].T
                                   )**2
                               )


            print(dm_to_decontam_jsd.median())
            print(dec_to_dm_jsd.median())


            all_debias_dec_jsds.append( dm_to_decontam_jsd )
            all_dec_debias_jsds.append( dec_to_dm_jsd )


            all_just_dec_jsds.append( 
                                    pd.Series( jensenshannon( samples_true.values[test_inds].T, 
                                                    dec_raw_out[:, 1:] [test_inds].T
                                   )**2
                               ) 
                                        )

            print( all_just_dec_jsds[-1].median() )

            all_raw_jsds.append( 
                     pd.Series( jensenshannon( samples_true.values[test_inds].T, 
                                                   samples_transformed.values[test_inds].T
                                   )**2
                               ) )


            print( all_raw_jsds[-1].median() )




    res_full= pd.DataFrame( all_run_descs )\
                            .assign(dm_jsd=all_debias_jsds)\
                            .assign(combined_jsd=all_csdmcr_jsds)\
                            .assign(separate_jsd=all_sep_jsds)\
                            .assign(all_debias_dec_jsds=all_debias_dec_jsds)\
                            .assign(all_dec_debias_jsds=all_dec_debias_jsds)\
                            .assign(raw_jsd=all_raw_jsds)\
                            .assign(dec_jsd=all_just_dec_jsds)


    name_map = {'raw_jsd':'Raw',
                'dec_jsd':'Decontam',
                'dm_jsd':'DEBIAS-M',
                'all_dec_debias_jsds':'Decontam --> DEBIAS-M',
                'all_debias_dec_jsds':'DEBIAS-M --> Decontam',
                'separate_jsd':'DEBIAS-M --> C-SCRuB',
                'combined_jsd': 'DM-C-SCRuB', 
                }

    plot_df = res_full\
            .explode(['dm_jsd', 
                      'combined_jsd', 
                      'separate_jsd', 
                      'all_debias_dec_jsds',
                      'all_dec_debias_jsds',
                      'raw_jsd', 
                      'dec_jsd'
                     ])\
                .melt(id_vars = res_full.columns[:2], 
                      value_vars=['dm_jsd', 
                                  'combined_jsd', 
                                  'separate_jsd',
                                  'all_debias_dec_jsds',
                                  'all_dec_debias_jsds', 
                                  'raw_jsd',
                                  'dec_jsd'
                                  ], 
                      value_name='jsd', 
                      var_name='method')\
            .groupby(['bias',
                      'method', 
                      'seed']
                        )['jsd'].median().reset_index()

    plot_df.method=plot_df.method.map(name_map)


    for hide_axes in [False, True]:
        out_path = '../results/plots/multibatch_simulations_varying_bias.pdf'
        plt.figure(figsize=(8,8))
        ax = sns.boxplot(x='bias', 
                         y='jsd', 
                         hue='method', 
                         data=plot_df, 
                         hue_order=[b for a,b in name_map.items()],#[1:]
                         fliersize=0, 
                         palette=palette
                         )

        ax = sns.swarmplot(x='bias', 
                           y='jsd', 
                           hue='method', 
                           data=plot_df, 
                           ax=ax, 
                           color='black', 
                           dodge=True, 
                           hue_order=[b for a,b in name_map.items()]#[1:]
                           )
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend(handles[:5], labels[:5])
        ax.set(yticks=[0.1, 0.3, .5])
        plt.ylim([0, 0.6])

        if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.title(None)
            ax.set_xticks([])
            ax.set(yticks=[0.1, 0.3, .5],
                   yticklabels=[])
            plt.legend().remove()
            out_path=out_path[:-4] + '-no-axes'+out_path[-4:]


        plt.savefig(out_path, 
                    format='pdf',
                    bbox_inches='tight', 
                    dpi=900)



if __name__=='__main__':
    main()


















