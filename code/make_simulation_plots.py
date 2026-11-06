
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
           'grid.color': 'black'}, 
       font_scale=2)


palette = {'No decontamination':'#1f77b4',
           'Decontam':'#55a868',
           'Decontam (Low Biomass)':'darkgreen',
           'Decontam Combined':'#078A69',
           'SCRuB':'#d62728',
           'B2':'white',
           'B1':'white',
           'Decontam Frequency':'#12A380',
           'C-SCRuB':'#e377c2'}


order=['No decontamination',
       'B1',
       'Decontam',
       'Decontam (Low Biomass)',
       'Decontam Combined',
       'SCRuB',
       'B2',
       'Decontam Frequency',
       'C-SCRuB'
       ]


x_val_map = {'../results/JSDs_simulation_data.csv':'Concentration Noise',
             '../results/JSDs_simulation_data_bias.csv':'log(Processing bias) STD',
             }

def main():
    for input_path in x_val_map:
        for hide_axes in [False, True]:
            all_dfs = pd.read_csv(input_path)
            
            if 'bias' in input_path:
                all_dfs = all_dfs.loc[all_dfs['Concentration Noise']>0]

            out_path = '../results/plots/{}.pdf'.format(input_path.split('/')[-1][:-4])

            ordd = all_dfs['Concentration Noise'].value_counts().index.sort_values()
            all_dfs = all_dfs.reset_index(drop=True)

            if 'bias' in input_path:
                plt.figure(figsize=( (4.1/3) * 8 * 4.9886 / 3.5, 8))
            else:
                
                plt.figure(figsize=(8 * 4.9886 / 3.5, 8))
                
            # plt.semilogy()
            ax=sns.boxplot(x='Concentration Noise', 
                           y='JSD', 
                           hue='Group', 
                           data=all_dfs, 
                           fliersize=0, 
                           hue_order=order,
                           order=ordd, 
                           palette=palette
                           )
            
            handles, labels = ax.get_legend_handles_labels()

            sns.swarmplot(x='Concentration Noise', 
                          y='JSD', 
                          hue='Group', 
                          data=all_dfs, 
                          dodge=True,
                          color='black', 
                          s=2.5,
                          ax=ax,
                          hue_order=order,
                          order=ordd
                          )

            plt.ylabel('Median JSD')
            plt.xlabel(x_val_map[input_path])

            ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
            plt.semilogy()

            ax.set_yticks([0.01, .1, .25])
            ax.set_ylim(top=0.9)

            if hide_axes:
                plt.xlabel(None)
                plt.ylabel(None)
                plt.title(None)
                ax.set_xticks([])
                ax.set(yticks=[0.01, .1, .25],
                       yticklabels=[])
                plt.legend().remove()
                out_path=out_path[:-4] + '-no-axes'+out_path[-4:]

            plt.savefig(out_path, 
                        format='pdf',
                        bbox_inches='tight', 
                        dpi=900)
            

if __name__=='__main__':
    main()
    
    
    
            
