"""Deterministic publication figures from executed endpoint and prediction tables."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/pktk-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/pktk-fontconfig')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np
import pandas as pd

def dose_ticks(ax, doses):
    """Label measured doses compactly without overlapping log-minor labels."""
    ticks = np.unique(np.asarray(doses, dtype=float))
    if len(ticks) > 1:
        span = np.log(ticks[-1] / ticks[0])
        chosen = [ticks[0]]
        for value in ticks[1:-1]:
            if np.log(value/chosen[-1]) >= .22*span and np.log(ticks[-1]/value) >= .22*span:
                chosen.append(value)
        ticks = np.asarray(chosen + [ticks[-1]])
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.3g}'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis='x', labelsize=7)

COLORS={'auc':'#0072B2','cmax':'#D55E00','linear':'#0072B2','power':'#D55E00','pooled_power':'#009E73','shrinkage_power':'#CC79A7'}

def figures(root,table,episodes,predictions,metrics,nonlinear):
    folder=root/'results/figures';folder.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.hashsalt':'pktk-research-v1','figure.dpi':120,'savefig.dpi':300})
    def save(fig,name):
        fig.savefig(folder/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(folder/(name+'.svg'),bbox_inches='tight',metadata={'Date':None})
        fig.savefig(folder/(name+'.png'),bbox_inches='tight',metadata={'Software':'pktk'})
        plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(11,8))
    unique=table.drop_duplicates('series_id')
    for ax,col,title in zip(axs.ravel(),['source','species','study_id','compound_name'],['Validated series by source','Validated series by species','Endpoint series by study','Endpoint series by compound']):
        if col in ['source','species']:
            counts=unique.groupby(col).size();counts.plot.bar(ax=ax,color='#0072B2');ax.set_ylabel('Series')
        else:
            counts=table.drop_duplicates('endpoint_series_id').groupby([col,'endpoint']).size().unstack(fill_value=0)
            if col=='study_id':counts.index=[s.replace('doi:','') for s in counts.index]
            counts.plot.barh(stacked=True,ax=ax,color=[COLORS[c] for c in counts.columns]);ax.tick_params(axis='y',labelsize=6);ax.set_xlabel('Endpoint-specific series')
        ax.set_title(title);ax.set_ylabel('' if col in ['study_id','compound_name'] else 'Series')
    fig.suptitle('Validated cohort; correlated strata are not independent studies',fontsize=12);fig.tight_layout();save(fig,'01_dataset_composition')
    names=sorted(table.compound_name.unique(),key=str.lower)
    for normalized,number,label in [(False,'02','dose_exposure'),(True,'03','dose_normalized_exposure')]:
        for endpoint in ['auc','cmax']:
            data=table[table.endpoint==endpoint];names_ep=[n for n in names if n in set(data.compound_name)]
            ncols=4;nrows=int(np.ceil(len(names_ep)/ncols));fig,axes=plt.subplots(nrows,ncols,figsize=(12,2.7*nrows),squeeze=False)
            for ax,name in zip(axes.ravel(),names_ep):
                d=data[data.compound_name==name]
                for sid,g in d.groupby('endpoint_series_id'):
                    g=g.sort_values('dose');y=g.exposure/g.dose if normalized else g.exposure
                    legend_label=f"{g.species.iloc[0][:3]}/{g.sex.iloc[0][0].upper()}/{'IV' if g.route.iloc[0]=='intravenous' else 'PO'}"
                    if d.vehicle.nunique()>1:legend_label+=' / '+('Cremophor/water' if 'Cremophor' in str(g.vehicle.iloc[0]) else str(g.vehicle.iloc[0])[:20])
                    if d.regimen.nunique()>1:legend_label+=' / '+str(g.regimen.iloc[0])
                    ax.plot(g.dose,y,'o-',ms=3,lw=.9,alpha=.8,label=legend_label)
                ax.set_xscale('log');ax.set_yscale('log');ax.set_title(name.strip(),fontsize=8);ax.grid(alpha=.15)
                dose_ticks(ax,d.dose)
                ax.set_xlabel('Dose (mg/kg)');ax.set_ylabel(('Exposure / dose' if normalized else endpoint.upper())+' (log scale)')
                ax.legend(fontsize=5,frameon=False,loc='best')
            for ax in axes.ravel()[len(names_ep):]:ax.set_visible(False)
            units='AUC: µg·h/mL' if endpoint=='auc' else 'Cmax: µg/mL'
            if normalized:units += ' per (mg/kg)'
            fig.suptitle(f'{endpoint.upper()}: '+('dose-normalized ' if normalized else '')+f'source group endpoints; {units}\nLines connect separate comparable strata; differing AUC definitions/matrices are not pooled',fontsize=10)
            fig.tight_layout(rect=(0,0,1,.95));save(fig,f'{number}_{label}_{endpoint}')
    good=predictions[(predictions.status=='ok') & (predictions.policy=='leave_one_dose_out') & (predictions.n_support==2)]
    fig,axes=plt.subplots(2,2,figsize=(10,8))
    for i,endpoint in enumerate(['auc','cmax']):
        for j,model in enumerate(['linear','power']):
            ax=axes[i,j];g=good[(good.endpoint==endpoint)&(good.model==model)]
            for task,marker in [('low_to_high','^'),('high_to_low','v'),('interpolation','o')]:
                v=g[g.task==task];ax.scatter(v.observed,v.predicted,s=18,alpha=.55,marker=marker,label=task.replace('_',' '))
            low=min(g.observed.min(),g.predicted.min());high=max(g.observed.max(),g.predicted.max())
            ax.plot([low,high],[low,high],color='black',lw=.8);ax.set_xscale('log');ax.set_yscale('log')
            unit='µg·h/mL' if endpoint=='auc' else 'µg/mL'
            ax.set_xlabel('Observed '+unit);ax.set_ylabel('Predicted '+unit);ax.set_title(f'{endpoint.upper()} · {model} · 2 supports');ax.legend(fontsize=6)
    fig.tight_layout();save(fig,'04_predicted_observed')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,endpoint in zip(axes,['auc','cmax']):
        for model in ['linear','power']:
            g=good[(good.endpoint==endpoint)&(good.model==model)&(good.task!='interpolation')]
            ax.scatter(g.dose_ratio,g.absolute_log_error,s=18,alpha=.5,color=COLORS[model],label=model)
        ax.axhline(np.log(2),c='gray',ls='--',lw=.8,label='twofold error');ax.set_xscale('log');ax.set_xlabel('Target / nearest support dose ratio (≥1)');ax.set_ylabel('Absolute log error');ax.set_title(endpoint.upper());ax.legend()
    fig.tight_layout();save(fig,'05_error_dose_ratio')
    for dimension,num in [('compound_name','06'),('study_id','07')]:
        fig,axes=plt.subplots(1,2,figsize=(12,6))
        for ax,endpoint in zip(axes,['auc','cmax']):
            g=good[(good.endpoint==endpoint)&(good.task=='low_to_high')&good.model.isin(['linear','power'])]
            scores=g.groupby([dimension,'model']).absolute_log_error.mean().unstack();scores= scores.sort_values('linear')
            if dimension=='study_id':scores.index=[s.replace('doi:','') for s in scores.index]
            scores.plot.barh(ax=ax,color=[COLORS[c] for c in scores.columns]);ax.set_ylabel('');ax.set_xlabel('Mean absolute log error');ax.set_title(endpoint.upper()+' · low→high · 2 supports');ax.tick_params(axis='y',labelsize=7)
        fig.tight_layout();save(fig,num+'_performance_by_'+dimension)
    # Post-hoc illustrative extremes, explicitly descriptive, one series per panel.
    chosen=[]
    for endpoint in ['auc','cmax']:
        n=nonlinear[(nonlinear.endpoint==endpoint)&(nonlinear.n_doses>=3)].copy();n['deviation']=abs(np.log(n.normalized_high_low_ratio))
        chosen.extend([n.sort_values(['deviation','endpoint_series_id']).iloc[0],n.sort_values(['deviation','endpoint_series_id']).iloc[-1]])
    fig,axes=plt.subplots(2,2,figsize=(10,7))
    for ax,r in zip(axes.ravel(),chosen):
        g=table[table.endpoint_series_id==r.endpoint_series_id].sort_values('dose');norm=g.exposure/g.dose
        ax.plot(g.dose,norm,'o-',c=COLORS[r.endpoint]);ax.axhline(norm.iloc[0],c='gray',ls='--',label='Low-dose proportional reference')
        ax.set_xscale('log');ax.set_yscale('log');dose_ticks(ax,g.dose)
        unit='µg·h/mL' if r.endpoint=='auc' else 'µg/mL'
        ax.set_xlabel('Dose (mg/kg)');ax.set_ylabel(r.endpoint.upper()+' / dose\n'+unit+' per (mg/kg)')
        ax.set_title(f'{r.compound_name.strip()} · {r.species}/{r.sex}\n{r.route} · {r.study_id} · β={r.power_exponent:.2f}, high/low normalized={r.normalized_high_low_ratio:.2f}',fontsize=9);ax.legend(fontsize=6)
    fig.suptitle('Post-hoc examples: near-proportional versus strong observed deviation',fontsize=11);fig.tight_layout(rect=(0,0,1,.95));save(fig,'08_proportionality_examples')
