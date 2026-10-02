"""Plot the PK covariate checks and grouped exposure-response results."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.facecolor':'white','pdf.fonttype':42})
NAVY='#244665';BRONZE='#a65c27'


def save(fig,name):
    for suffix in ['png','pdf']:
        fig.savefig(ROOT/'figures'/f'{name}.{suffix}',dpi=200,bbox_inches='tight')
    plt.close(fig)


def pk_checks():
    data=pd.read_csv(ROOT/'results/pk_covariate_checks.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,6),sharey=True)
    y=np.arange(len(data))
    for ax,key,title in zip(axes,['AUC_tau','Cmax'],['AUC ratio','Peak concentration ratio']):
        ax.scatter(data[key+'_ratio_published'],y+.10,label='Published point',color=NAVY,s=42,marker='o')
        ax.scatter(data[key+'_ratio_model'],y-.10,label='Implementation',color=BRONZE,s=42,marker='x')
        ax.axvline(1,color='#aaa',linewidth=.8)
        ax.set_xlabel(title);ax.set_xlim(.45,1.5);ax.grid(axis='x',alpha=.18)
    axes[0].set_yticks(y,data.contrast);axes[0].invert_yaxis()
    axes[1].legend(loc='lower right',fontsize=9)
    fig.suptitle('Selected PK covariate point contrasts',fontsize=14)
    fig.text(.5,.012,'960 mg once daily, day 30. Point agreement does not establish clinical validity.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.045,1,.95]);save(fig,'pk_covariate_checks')


def grouped_sensitivity():
    data=pd.read_csv(ROOT/'data/fda_figure26_quartiles.csv')
    fit=pd.read_csv(ROOT/'results/grouped_coordinate_sensitivity.csv')
    p=data.responders/data.total
    low=beta.ppf(.025,data.responders,data.total-data.responders+1)
    high=beta.ppf(.975,data.responders+1,data.total-data.responders)
    fig,(left,right)=plt.subplots(1,2,figsize=(10,4.9),gridspec_kw={'width_ratios':[1,1.5]})
    x=np.arange(4)
    left.errorbar(x,100*p,yerr=[100*(p-low),100*(high-p)],fmt='o',color=NAVY,capsize=4)
    left.set_xticks(x,['Q1','Q2','Q3','Q4']);left.set_ylim(0,65)
    left.set_xlabel('Predicted AUC quartile');left.set_ylabel('Objective response (%)')
    left.set_title('Published grouped counts',fontsize=11)
    for i,row in enumerate(data.itertuples()):left.text(i,100*p.iloc[i]+13,f'{row.responders}/{row.total}',ha='center',fontsize=9)
    names=['Displayed labels','Geometric midpoints','Arithmetic midpoints']
    y=np.arange(3)
    right.errorbar(fit.OR_per_doubling,y,
        xerr=[fit.OR_per_doubling-fit.conditional_95CI_lower,fit.conditional_95CI_upper-fit.OR_per_doubling],
        fmt='o',color=BRONZE,capsize=4)
    right.set_yticks(y,names);right.invert_yaxis();right.set_ylim(2.5,-.5)
    right.axvline(1,color='#777',ls=':',lw=1);right.set_xlim(.30,1.07)
    right.set_xlabel('OR per doubling of assigned AUC coordinate')
    right.set_title('Grouped-model representation sensitivity',fontsize=11)
    fig.suptitle('Coordinate sensitivity in grouped exposure-response analysis',fontsize=13)
    fig.text(.5,.035,'Left: exact binomial 95% intervals. Right: conditional model 95% intervals.\nWithin-bin exposure uncertainty and confounding are not removed.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.105,1,.93],w_pad=2.5);save(fig,'grouped_coordinate_sensitivity')


if __name__=='__main__':
    (ROOT/'figures').mkdir(exist_ok=True)
    pk_checks();grouped_sensitivity()
    print('Generated two retained figures, each as PNG and PDF.')
