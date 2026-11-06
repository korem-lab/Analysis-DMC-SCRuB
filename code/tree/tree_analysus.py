import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from skbio.diversity.alpha import shannon
import os
from sklearn.linear_model import LogisticRegressionCV, RidgeCV
from sklearn.model_selection import cross_val_predict, cross_val_score, LeaveOneGroupOut
from skbio.stats.composition import clr
from debiasm.torch_functions import rescale
import tqdm
import torch
from debiasm import DebiasMClassifier
from sklearn.metrics import roc_auc_score

from sklearn.metrics import roc_auc_score
import sys
sys.path.append('../Clean_analysis/NEWER_DM_C_SCRuB/c_scrub/')
from sklearn_functions import C_SCRuB_DebiasMClassifier

from scipy.stats import combine_pvalues
import delong

def one_tail_delong(y, p1, p2):
    delong_p = 10**( delong.delong_roc_test( y, p1, p2 ) )/2
    
    if roc_auc_score(y, p1) > roc_auc_score(y, p2):
        delong_p = 1 - delong_p
    
    return(delong_p)

def main():

    otu_tbl = pd.read_csv('tree-microbiome/16S/OTU_table.txt', sep='\t', index_col=0)
    otu_tbl = otu_tbl.loc[(otu_tbl.Genus.isna()==False)&\
                          (otu_tbl.Species.isna())]\
                    .fillna('NaN').groupby(['Kingdom', 'Phylum', 'Class',
                     'Order', 'Family', 
                     'Genus', 'Species'])[otu_tbl.columns[:-7]].sum().T
    md_16s = pd.read_csv('tree-microbiome/29205563/16s_w_metadata.csv')

    name_map = pd.DataFrame( {a:b.values[0] for a,b in
                {x: otu_tbl.index[otu_tbl.index.str.contains(x)] for x in md_16s['Sample ID'].unique() }.items()
                 if b.shape[0]==1 }, index=['nm']).T.reset_index()


    otu_tbl=otu_tbl.loc[name_map.nm]
    md_16s = md_16s.set_index('Sample ID').loc[name_map['index']]


    md_full=pd.read_csv('tree-microbiome/29205563/16S_tree_sample_table_with_meta.csv', index_col=0)\
                .loc[name_map.nm]
    md_full=md_full.loc[md_full['plate_identifier.x'].isna()==False]


    pcr_md=pd.read_csv('tree-microbiome/29205563/ddPCR_meta_all_data.csv', index_col=0)\

    pcr_md=pcr_md.groupby('plate_identifier').head(1)

    pcr_md=pcr_md.set_index('plate_identifier').loc[md_full['plate_identifier.x']]
    pcr_md.index=md_full.index

    inds = ( (md_full['material'] == 'Wood')&(md_full['Tree.Core.Species.Code'].isna()==False)
           ).values

    otu_tbl = otu_tbl.loc[:, ( otu_tbl[inds] > 0 ).mean(axis=0) > 0.01 ]
    inds = otu_tbl.sum(axis=1)>0
    md_full=md_full.loc[inds]
    otu_tbl=otu_tbl.loc[inds]

    X_with_batch = np.hstack((pd.Categorical(md_full['Extraction.Plate.ID.x']).codes[:, np.newaxis], 
                                   otu_tbl.values))

    inds = ( (md_full['material'] == 'Wood')&(md_full['Tree.Core.Species.Code'].isna()==False)
           ).values


    pd.DataFrame({'bt':md_full['Extraction.Plate.ID.x'].values[inds], 
                  'mu':md_full['X16S_per_ul'].values[inds]
                 }).to_csv('tmp-dataset/md.csv')

    pd.DataFrame( rescale( X_with_batch[:, 1:][inds] ) ).to_csv('tmp-dataset/data.csv')

    ### run the R script
    os.system('conda run -p /Users/george/opt/anaconda3/envs/r43 Rscript run_simulation_decontam.R')

    decontam_out = pd.read_csv('tmp-dataset/decontam_raw_out.csv', index_col=0).values

    decontam_out_with_batch = np.hstack((X_with_batch[inds][:, :1], decontam_out))


    seed=42
    disease_types = pd.Series( md_full['Tree.Core.Species.Code'].values[inds] ).value_counts()\
        .head(4).index ## min sample count

    complete_run_descs = []
    all_preds = []

    for disease_type in disease_types:
        print(disease_type)
        for test_plate in np.unique( md_full['Extraction.Plate.ID.x'].values[inds] ):
            print(test_plate)


            test_inds = md_full['Extraction.Plate.ID.x'].values[inds] == test_plate
            if ( (md_full['Tree.Core.Species.Code'].values[inds][test_inds]==disease_type).sum()>5 ):
                if True:
                    pred_summaries = {}


                    np.random.seed(seed)
                    torch.manual_seed(seed)

                    X_with_batch[inds][:, 0] = pd.Categorical( X_with_batch[inds][:, 0] ).codes


                    dmc = DebiasMClassifier(x_val=X_with_batch[inds][test_inds].astype(float))

                    dmc = dmc.fit(X_with_batch[inds][~test_inds].astype(float), 
                                  md_full['Tree.Core.Species.Code'].values[inds][~test_inds]==disease_type,
                                  )

                    dmc_auroc = roc_auc_score(md_full['Tree.Core.Species.Code'].values[inds][test_inds]==disease_type,
                                              dmc.val_preds)

                    pred_summaries['Debias'] = dmc.val_preds



                    np.random.seed(seed)
                    torch.manual_seed(seed)

                    cscr_dm = C_SCRuB_DebiasMClassifier(x_val=X_with_batch[inds][test_inds].astype(float), 
                                              mu_val=md_full['X16S_per_ul'].values[inds][test_inds])

                    cscr_dm=cscr_dm.fit(X_with_batch[inds][~test_inds].astype(float), 
                                        md_full['Tree.Core.Species.Code'].values[inds][~test_inds]==disease_type,
                                        md_full['X16S_per_ul'].values[inds][~test_inds].astype(float)
                                        )

                    cscr_auroc = roc_auc_score(md_full['Tree.Core.Species.Code'].values[inds][test_inds]==disease_type,
                                          cscr_dm.val_preds)

                    pred_summaries['DMC'] = cscr_dm.val_preds


                    ## deonctam
                    dmc = DebiasMClassifier(x_val=decontam_out_with_batch[test_inds].astype(float))

                    dmc = dmc.fit(decontam_out_with_batch[~test_inds].astype(float), 
                                  md_full['Tree.Core.Species.Code'].values[inds][~test_inds]==disease_type,
                                  )

                    dec_deb_auroc = roc_auc_score(md_full['Tree.Core.Species.Code'].values[inds][test_inds]==disease_type,
                                                  dmc.val_preds)

                    pred_summaries['Dec-DEB'] = dmc.val_preds



                    np.random.seed(seed)
                    torch.manual_seed(seed)

                    lr = LogisticRegressionCV()

                    lr=lr.fit(
                         clr(1e-6+
                            rescale(X_with_batch[inds][~test_inds].astype(float)[:, 1:] )
                            ), 
                                        md_full['Tree.Core.Species.Code'].values[inds][~test_inds]==disease_type
                                        )

                    lr_auroc = roc_auc_score(md_full['Tree.Core.Species.Code'].values[inds][test_inds]==disease_type,
                                               lr.predict_proba(
                                                     clr(1e-6+
                                                   rescale( X_with_batch[inds][test_inds].astype(float)[:, 1:]) )
                                               )[:, 1]
                                                )


                    pred_summaries['LR'] = lr.predict_proba(
                                                     clr(1e-6+
                                                   rescale( X_with_batch[inds][test_inds].astype(float)[:, 1:]) )
                                               )[:, 1]


                    np.random.seed(seed)
                    torch.manual_seed(seed)

                    lr = LogisticRegressionCV()

                    lr=lr.fit(
                          clr(1e-6+
                              rescale( decontam_out_with_batch[~test_inds].astype(float)[:, 1:] )
                             ), 
                                        md_full['Tree.Core.Species.Code'].values[inds][~test_inds]==disease_type
                                        )

                    lr_dec_auroc = roc_auc_score(md_full['Tree.Core.Species.Code'].values[inds][test_inds]==disease_type,
                                               lr.predict_proba(
                                  clr(1e-6+
                                      rescale( decontam_out_with_batch[test_inds].astype(float)[:, 1:])
                                     ) )[:, 1]
                                               )

                    pred_summaries['LR_Dec'] = lr.predict_proba( clr(1e-6+
                                      rescale( decontam_out_with_batch[test_inds].astype(float)[:, 1:])
                                     ) )[:, 1]


                    pred_summaries['truth'] = md_full['Tree.Core.Species.Code'].values[inds][test_inds]==disease_type


                    run_desc = {'Species type':disease_type, 
                                'Test plate':test_plate,  
                                'DEBIAS-M auROC':dmc_auroc, 
                                'C-SCRuB-DM auROC':cscr_auroc,
                                'Decontam auROC':dec_deb_auroc, 
                                'Linear auROC':lr_auroc, 
                                'Linear decontam auROC':lr_dec_auroc,
                                }


                    all_preds.append(pred_summaries)


                    print(run_desc)

                    complete_run_descs.append(run_desc)



    sns.set(rc={'axes.facecolor':'white', 
                'figure.facecolor':'white', 
                'axes.edgecolor':'black', 
                'grid.color': 'black'
                }, 

           font_scale=2)

    sns.set_style('ticks')

    palette = {'Raw':'#1f77b4',
               'Linear auROC':'#1f77b4',
               'Decontam':'#55a868',
               'Decontam (Low Biomass)':'darkgreen',
               'Decontam Combined':'#078A69',
               'DEBIAS-M auROC':'lightblue',
               'SCRuB':'#d62728',
               'B2':'white',
               'B1':'white',
               'Decontam --> DEBIAS-M auROC':'#12A380',
               'DEBIAS-M --> Decontam':'#12A380',
               'DEBIAS-M --> C-SCRuB':'#e377c2', 
               'C-SCRuB-DM auROC':'#FF6065', 
               'Linear decontam auROC':'#12A380',
               'DMC --> CLR':'#FF6065', 
               }

    for hide_axes in [False, True]:
        out_path = '../../results/plots/Tree-16S-auROC-clr-linears.pdf'

        comp_res = pd.DataFrame(complete_run_descs)


        plot_df = pd.melt(comp_res, 
                          value_vars=comp_res.columns[2:], 
                          id_vars=comp_res.columns[:2],
                          var_name='Group',
                          value_name='auROC'
                          )

        plot_df['Group'] = plot_df.Group.str.replace('Decontam', 'Decontam --> DEBIAS-M', )
        hue_order = plot_df.Group.unique()[[3,4,0,2,1]]

        plt.figure(figsize=(8,8*1.6))
        ax=sns.boxplot(x='Group',
                       y='auROC', 
                       data=plot_df, 
                       order=hue_order, 
                       palette=palette
                       )

        ax=sns.swarmplot(x='Group',
                         y='auROC', 
                         data=plot_df, 
                         s=10, 
                         order=hue_order,
                         color='black'
                         )


        plt.xlabel(None)
        plt.title('Predicting tree species from\ntheir wood composition')
        plt.xticks(rotation=90)

        plt.ylim(0.28, 1.075)

        ax.set_yticks([0.3, .5, .7, 0.9])

        if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.title(None)
            ax.set_xticks([])
            ax.set(yticks=[0.3, .5, .7, 0.9],
                   yticklabels=[])
            plt.legend().remove()
            out_path=out_path[:-4] + '-no-axes'+out_path[-4:]



        plt.savefig(out_path, 
                    format='pdf',
                    bbox_inches='tight', 
                    dpi=900)


    print({a:
     combine_pvalues(
            np.hstack( [one_tail_delong(all_preds[i]['truth'].astype(int),
                                        all_preds[i][a], 
                                        all_preds[i]['DMC']
                                        )
                         for i in range(len(all_preds)) ] )[0] )
            for a in all_preds[0] if a not in ['DMC', 'truth'] 
            })



if __name__=='__main__':
    main()




