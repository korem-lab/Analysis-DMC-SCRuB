
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import torch
import statsmodels
from scipy import stats
from scipy.spatial.distance import jensenshannon
import os
import glob

from c_scrub import C_SCRuB


def main(rd_path):
    all_res_dirs =  glob.glob(rd_path)


    all_raw_jsds = {}
    all_scrub_jsds = {}
    all_cscrub_jsds = {}

    ii=0
    for sim_dir in all_res_dirs:
        torch.manual_seed(42)
        np.random.seed(42)
        ii+=1
        print(ii)
        df = pd.read_csv(os.path.join(sim_dir, 'observed_samples.csv'), index_col=0)
        conc = pd.read_csv(os.path.join(sim_dir, 'true_concentrations.csv'), index_col=0)


        cscr_out={'decontaminated_samples':C_SCRuB().fit_transform(df, conc.x)}

        true_samples = pd.read_csv(os.path.join(sim_dir, 
                                                'true_samples.csv'), index_col=0)

        scr_out = pd.read_csv(os.path.join(sim_dir, 'scrub_out.csv'), index_col=0)

        run_nm=sim_dir.split('/')[-1]

        all_raw_jsds[run_nm] = jensenshannon(df.T, 
                                             true_samples.T)**2

        all_scrub_jsds[run_nm] = jensenshannon(scr_out.T, 
                                               true_samples.T)**2

        all_cscrub_jsds[run_nm] = jensenshannon(cscr_out['decontaminated_samples'].T, 
                                                true_samples.T)**2


    raw_df = pd.DataFrame({a:np.median(b) for a,b in all_raw_jsds.items()}, index=['JSD']).T\
                    .assign(Group='Raw')
    scr_df = pd.DataFrame({a:np.median(b) for a,b in all_scrub_jsds.items()}, index=['JSD']).T\
                    .assign(Group='SCRuB')
    csc_df = pd.DataFrame({a:np.median(b) for a,b in all_cscrub_jsds.items()}, index=['JSD']).T\
                    .assign(Group='C-SCRuB')



    all_dfs = pd.concat([raw_df, scr_df, csc_df])
    all_dfs['round']=all_dfs.index.str.split('_').str[-1].astype(int)

    full_res = pd.read_csv(os.path.join('/'.join(rd_path.split('/')[:-1]), 
                                        'concentration_simulation_results.csv'), index_col=0)

    all_dfs = all_dfs.merge(full_res.drop(['names', 'jsds', 'contam_species'], axis=1).drop_duplicates(), 
                            on='round')

    ratio = 0.69314718
    full_res['JSD']=full_res.jsds * ratio
    full_res['Group']=full_res['names']

    all_dfs = pd.concat([all_dfs, 
                         full_res.loc[full_res.names.isin(['Input', 'SCRUB'])==False]\
                             .groupby(list(all_dfs.columns.values[1:]))['JSD'].median()\
                             .reset_index()[all_dfs.columns]
                         ] )


    ## the w2w level was just the placehold column i tossed this in
    all_dfs['Concentration Noise'] = all_dfs.well_to_well_level.astype(str)

    all_dfs['Concentration STD'] = all_dfs['contam_level']

    all_dfs.Group = all_dfs.Group.str.replace('Restrictive', 'Decontam Frequency')\
                           .str.replace('microDecon', 'Decontam Combined')\
                           .str.replace('Raw', 'No decontamination')

    sns.set_theme(font_scale=2)

    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
                'axes.edgecolor':'black', 
               'grid.color': 'black'}, 
           font_scale=2)

    order=['No decontamination',
           'Decontam',
           'Decontam (Low Biomass)',
           'Decontam Frequency',
           'Decontam Combined',
           'SCRuB',
           'C-SCRuB'
           ]

    ordd = all_dfs['Concentration Noise'].value_counts().index.sort_values()
    all_dfs = all_dfs.reset_index(drop=True)
    
    return(all_dfs)


if __name__=='__main__':
    for rd_path in [
                    '../results/simulation_data/round_*',
                    '../results/simulation_data_bias/round_*']:
        
        results = main(rd_path)
        results.to_csv('../results/JSDs_{}.csv'.format(rd_path.split('/')[2]))
        



    
    
    
    
    
    
    
    
    
    
    
    