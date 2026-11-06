import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.compose import TransformedTargetRegressor

sns.set(rc={'axes.facecolor':'white', 
            'figure.facecolor':'white', 
            'axes.edgecolor':'black', 
            'grid.color': 'black'
            }, 
        
       font_scale=2)

sns.set_style('ticks')


def hivrc():
    hivrc_bcfs = pd.read_csv('../../../COMPLETED_PUBLISHED_PAPER_STUFFF/DEBIAS-M nat micro/DEBIAS-M-revisions/results/HIVRC/processed_hiv_bcfs.csv', 
                             index_col=0)
    ww=hivrc_bcfs.values.ravel()

    for hide_axes in [False, True]:

        out_path='../results/plots/HIVRC_BCF_histogram.pdf'

        plt.figure(figsize=(8,8))
        plt.semilogy()
        ax = sns.histplot( x=ww, 
                          bins=50, 
                          color='#'+'3'*6
                           )
        plt.yticks([10, 100, 1000])
        plt.minorticks_off()  



        if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.title(None)
            ax.set(yticks=[10, 100, 1000],
                   yticklabels=[], 
                   xticks=[ 0, 5, 10 ], 
                   xticklabels=[]
                   )
            plt.legend().remove()
            out_path=out_path[:-4] + '-no-axes'+out_path[-4:]



        plt.savefig(out_path, 
                    format='pdf',
                    bbox_inches='tight', 
                    dpi=900)

def ex_3c():        

    df = pd.read_csv('../results/Example_for_3c.csv', index_col=0)
    df_long = pd.melt(
        df, 
        id_vars=['conc'],  
        value_vars=['contam_level', 'biased_contam_level'],   
        var_name='Group',   
        value_name='contam'  
    )        

    for hide_axes in [False, True]:

        out_path='../results/plots/Example_for_3c.pdf'

        plt.figure(figsize=( ( 3.4776 / 4.9886 )   * 8 * 4.9886 / 3.5, \
                             8 * 4.3242 / 3.5 ))
        plt.semilogx()
        ax=sns.scatterplot(data=df_long,
                        x='conc', 
                        y='contam', 
                        hue='Group', 
                        palette={'contam_level':'#D988B9', 
                                 'biased_contam_level':'#F5077A'},
                        s=75
                       )
        ax.minorticks_off()

        plt.yticks([0, .25, .5])
        plt.xticks([10**7, 10**8, 10**9])


        if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.title(None)
            ax.set(yticks=[0, .25, .5],
                   yticklabels=[], 
                   xticks=[10**7, 10**8, 10**9], 
                   xticklabels=[]
                   )
            plt.legend().remove()
            out_path=out_path[:-4] + '-no-axes'+out_path[-4:]


        plt.savefig(out_path, 
                    format='pdf',
                    bbox_inches='tight', 
                    dpi=900)
    return(None)
    
        
 def karstens():

    m1 = pd.read_csv('../../Fig1/dilution_summary.csv', index_col=0)
    m2 = pd.read_csv('../../Fig1/msystems.00290-19-st001.csv', index_col=0)

    x = m2['Concentration of samples (cells/ml)']
    y = m1['PercentContaminants'].loc[x.index]
    inverse_sigmoid = lambda p: np.log(p / (1-p))
    import math
    def sigmoid(x):
        return 1 / (1 + np.exp(-x))
    lr=LinearRegression()
    lr.fit(np.log(x)[:, np.newaxis],
           inverse_sigmoid(y/100))

    for hide_axes in [False, True]:

        out_path = '../results/plots/background-contam-vs-conc.pdf'

        plt.figure(figsize=(6,6))
        plt.semilogx()
        ax = sns.scatterplot(x=x, 
                             y=y, 
                             s=250, 
                             color='black')

        ax = sns.lineplot(x=np.logspace(6, 11, 100),
                          y = 100*sigmoid( lr.predict(np.log(np.logspace(6, 11, 100))[:, np.newaxis])), 
                          color='red', 
                          linewidth=2.5
                            )

        plt.yticks([0, 40, 80])
        plt.xticks([ 10**6, 10**8, 10**10])


        if hide_axes:
            plt.xlabel(None)
            plt.ylabel(None)
            plt.title(None)
            ax.set(yticks=[0, 40, 80],
                   yticklabels=[], 
                   xticks=[ 10**6, 10**8, 10**10], 
                   xticklabels=[]
                   )
            plt.legend().remove()
            out_path=out_path[:-4] + '-no-axes'+out_path[-4:]



        plt.savefig(out_path, 
                    format='pdf',
                    bbox_inches='tight', 
                    dpi=900)

        
def ex_3b():
    true_data = np.array([[.38, .38, .24], 
                          [.3, .25, .33], 
                          [.75, .05, .2],
                          [0, .5, .5]
                          ]) 

    biases = np.array([.5, 3, 2])
    concentrations = 1 / ( 0.1 + true_data[:, -1] )

    biased_data = rescale( true_data * biases )
    n_linsp_pts=1000
    for hide_axes in [False, True]:
        out_path = '../results/plots/demo_example.pdf'
  

        plt.figure(figsize=( ( 3.4776 / 4.9886 )   * 8 * 4.9886 / 3.5, \
                             8 * 4.3242 / 3.5 ))
        ax = sns.scatterplot(x = concentrations,
                             y = true_data[:, -1], 
                             s = 500, 
                             color = '#D988B9', 
                             marker='D', 
                             edgecolor="none"
                             )

        ax = sns.lineplot( x = np.linspace(concentrations.min(), 
                                                          concentrations.max(), 
                                                          n_linsp_pts),
                           y = TransformedTargetRegressor(
                                            regressor=LinearRegression(), 
                                            func=np.log, 
                                            inverse_func=np.exp
                                           ).fit(concentrations[:, np.newaxis],
                                                 true_data[:, -1][:, np.newaxis],
                                                 )\
                                    .predict( np.linspace(concentrations.min(), 
                                                          concentrations.max(), 
                                                          n_linsp_pts)[:, np.newaxis]
                                            )[:, 0], 
                           ax = ax, 
                           color='black',
                           linewidth=2.5, 
                           linestyle='--'
                          )


        ax = sns.scatterplot(x = concentrations,
                             y = biased_data[:, -1], 
                             s = 500, 
                             color = '#D988B9', 
                             ax = ax, 
                             marker='x', 
                             linewidth=5, 
                             edgecolor="none"
                             )

        ax = sns.lineplot( x = np.linspace(concentrations.min(), 
                                            concentrations.max(), 
                                            n_linsp_pts),
                           y = TransformedTargetRegressor(
                                           regressor=LinearRegression(), 
                                            func=np.log, 
                                            inverse_func=np.exp
                                           ).fit(concentrations[:, np.newaxis],
                                                 biased_data[:, -1][:, np.newaxis],
                                                 )\
                                    .predict( np.linspace(concentrations.min(), 
                                                          concentrations.max(), 
                                                          n_linsp_pts)[:, np.newaxis]
                                            )[:, 0], 
                           ax = ax, 
                           color='black',
                           linewidth=2.5, 
                           linestyle='--'
                          )


        plt.yticks([0.2, .3, .4, .5])

        if hide_axes: 
            ax.set_xlabel(None)
            ax.set_ylabel(None)
            ax.set(yticks=ax.get_yticks(),
                   xticks=ax.get_xticks(),
                   yticklabels=[], 
                   xticklabels=[]
                   )
            plt.title(None)
            out_path=out_path[:-4] + '-no-axes'+  out_path[-4:]

        plt.savefig(out_path, 
                    format='pdf',
                    bbox_inches='tight', 
                    dpi=900
                    )


if __name__=='__main__':
    
    hivrc()
    ex_3c()
    ex_3b()
    karstens()












        
        