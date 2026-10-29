"""Render a deterministic two-page research brief from independently verified results.

No collection or model fitting. Source figures retain their original data licenses.
"""
import json
import os
from pathlib import Path
import textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/pktk-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader
from verify_final_results import verify, ROOT

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False})
BLUE='#215678'
GRAY='#444444'

def text(fig,x,y,value,size=10,bold=False,width=100,color=GRAY):
    lines='\n'.join(textwrap.fill(p,width=width) if p else '' for p in value.split('\n'))
    fig.text(x,y,lines,ha='left',va='top',fontsize=size,color=color,
             weight='bold' if bold else 'normal',linespacing=1.45)


def page(number,title,subtitle):
    fig=plt.figure(figsize=(8.27,11.69),facecolor='white')
    text(fig,.08,.955,'PRECLINICAL PK/TK  |  PUBLIC-DATA RESEARCH',9,True,color=BLUE)
    text(fig,.08,.913,title,20,True,width=60,color='black')
    text(fig,.08,.867,subtitle,10,width=95)
    fig.add_artist(plt.Line2D([.08,.92],[.835,.835],color=BLUE,lw=1))
    text(fig,.08,.025,'Independent academic exploration • Verified snapshot 2026-10-29',8)
    text(fig,.88,.025,str(number),9)
    return fig


def main():
    v=verify();a=v['dataset_a'];b=v['dataset_b']
    out=ROOT/'results/finalization';out.mkdir(parents=True,exist_ok=True)
    meta={'Title':'Preclinical dose–exposure prediction: research brief',
          'Author':'PK/TK public-data research project','CreationDate':None,'ModDate':None}
    with PdfPages(ROOT/'reports/pk_tk_research_brief.pdf',metadata=meta) as pdf:
        fig=page(1,'From public endpoints to comparable doses',
                 'When can initial PK measurements predict exposure at an unobserved dose?')
        text(fig,.08,.812,'Research motivation',12,True,color=BLUE)
        text(fig,.08,.784,'A useful cross-dose benchmark needs more than molecular coverage: it needs matched experimental conditions, explicit dose levels, and a genuinely held-out exposure. This study tests how far simple empirical scaling can go.',width=98)
        text(fig,.08,.703,'Dataset A • molecular PK coverage',12,True,color=BLUE)
        text(fig,.08,.674,f"{a['unique_structured_compounds']:,} structure-resolved compounds • {a['eligible_endpoint_records']:,} eligible endpoint records\n5,133 excluded records retained for audit; 355 compounds have both AUC and Cmax.",width=102)
        ax=fig.add_axes([.24,.445,.63,.17])
        labels=['AUC','Cmax','Clearance','Volume','Half-life']
        vals=[a[k] for k in ['compounds_with_auc','compounds_with_cmax','compounds_with_clearance','compounds_with_volume','compounds_with_half_life']]
        ax.barh(labels[::-1],vals[::-1],color=BLUE,height=.6)
        for i,n in enumerate(vals[::-1]):ax.text(n+20,i,f'{n:,}',va='center',fontsize=9)
        ax.set_xlim(0,1450);ax.set_xlabel('Unique structured compounds (overlapping endpoint coverage)',fontsize=9)
        ax.tick_params(labelsize=9);ax.spines[['top','right']].set_visible(False)
        text(fig,.08,.397,'Dataset B • strict cross-dose benchmark',12,True,color=BLUE)
        text(fig,.08,.368,'59 comparable series • 21 compounds • 19 independent source studies\n42 AUC-eligible series • 33 Cmax-eligible series • 16 eligible for both\n19 compounds have at least three comparable doses; 19 have resolved structures.',width=105)
        text(fig,.08,.282,'Provenance and comparability',12,True,color=BLUE)
        text(fig,.08,.253,'A: ChEMBL 37 rat ADME tranche, measured-animal PKSmart Table S2, and structured B endpoints. B: CEBS (47 series) and reviewed publications (12 series). Source IDs, hashes, original units and QC decisions are preserved. Route, matrix, AUC window and protocol differences remain separate.',width=101)
        text(fig,.08,.159,'The distinction matters',12,True,color=BLUE)
        text(fig,.08,.131,'Dataset A is not cross-dose validated. Its 627 publication/study proxies are not verified independent experiments. Dataset B supports 513 episodes, but only 19 linked compound–study evaluation units overall.',width=101)
        text(fig,.08,.071,'Sources: ChEMBL (CC BY-SA 3.0); PKSmart, doi:10.1186/s13321-025-01066-5 (CC BY 4.0).\nCEBS/PMC citations and article-specific terms: reports/data_inventory.md. No raw dataset is included.',7.5,width=120)
        pdf.savefig(fig);fig.savefig(out/'brief_page_1.png',dpi=160);plt.close(fig)
        fig=page(2,'Empirical extrapolation: what the data show',
                 'Two support doses → a higher held-out dose • AUC and Cmax evaluated separately')
        text(fig,.08,.812,'Verified baseline comparison',12,True,color=BLUE)
        ax=fig.add_axes([.08,.651,.84,.137]);ax.axis('off')
        rows=[]
        for r in v['baseline_results']:
            rows.append(['AUC' if r['endpoint']=='auc' else 'Cmax','Linear' if r['model']=='linear' else 'Power law',
                         f"{r['mae_ln']:.3f}",f"{r['rmse_ln']:.3f}",f"{100*r['mean_absolute_relative_error']:.1f}%",str(r['n_components'])])
        table=ax.table(cellText=rows,colLabels=['Endpoint','Model','Log MAE','Log RMSE','Mean ARE','Units'],loc='center',cellLoc='center',colWidths=[.15,.21,.16,.17,.18,.13])
        table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.5)
        for (i,j),cell in table.get_celld().items():
            cell.set_edgecolor('#dddddd');cell.set_linewidth(.5)
            if i==0:cell.set_facecolor('#eaf0f4');cell.set_text_props(weight='bold')
        text(fig,.08,.643,'Natural logs; ARE = absolute relative error. Equal series/component weighting.\nOnly 11 independent AUC units and 12 Cmax units; repeated episodes are not replicates.',8,width=116)
        ax=fig.add_axes([.04,.303,.92,.302]);ax.imshow(plt.imread(ROOT/'results/figures/06_performance_by_compound_name.png'));ax.axis('off')
        text(fig,.08,.301,'Compound-level errors vary substantially. Existing figure 06; selected, correlated strata.\nMethyl eugenol is favorable; bromodichloromethane AUC and HY-071085 Cmax are challenging.',8,width=117)
        text(fig,.08,.245,'Interpretation and limits',12,True,color=BLUE)
        text(fig,.08,.216,'Power-minus-linear log-MAE differences: AUC +0.058 (95% exploratory interval −0.176 to +0.308); Cmax −0.072 (−0.391 to +0.231). Both include zero. Power law does not consistently outperform linear scaling; nonlinear trends do not identify a mechanism.',9.5,width=103)
        text(fig,.08,.135,'Next steps',12,True,color=BLUE)
        text(fig,.08,.108,'Acquire more independent, comparable dose series and animal-level uncertainty before evaluating molecular or mechanistic hybrid models. Complex molecular ML was not evaluated; robust scaffold tests and calibrated prediction intervals remain unsupported.',9,width=111)
        text(fig,.08,.057,'Figure data: CEBS and source publications listed in data_inventory.md; includes lorlatinib\n(PMC12570999, CC BY-NC 4.0). Retain source-specific terms; figure is not covered by the code MIT license.',7,width=129)
        pdf.savefig(fig);fig.savefig(out/'brief_page_2.png',dpi=160);plt.close(fig)
    pdf=PdfReader(ROOT/'reports/pk_tk_research_brief.pdf')
    assert len(pdf.pages)==2
    assert '2,713' in pdf.pages[0].extract_text() and '0.525' in pdf.pages[1].extract_text()
    print('Rendered and checked two-page PDF; previews in results/finalization/')

if __name__=='__main__':main()
