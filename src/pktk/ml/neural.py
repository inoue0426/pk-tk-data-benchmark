"""Small CPU graph networks with explicit validation-only early stopping."""
import copy,random,time
import numpy as np
import torch
from torch import nn
from torch_geometric.nn import GCNConv,GINConv,GINEConv,global_mean_pool,global_add_pool
from torch_geometric.loader import DataLoader


def seed_all(seed):
    random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)


class GraphModel(nn.Module):
    def __init__(self,kind,nodes,edges,hidden=32,ablate=False):
        super().__init__();self.kind=kind;self.ablate=ablate
        self.input=nn.Linear(nodes,hidden);self.layers=nn.ModuleList()
        for _ in range(2):
            mlp=nn.Sequential(nn.Linear(hidden,hidden),nn.ReLU(),nn.Linear(hidden,hidden))
            self.layers.append(GCNConv(hidden,hidden) if kind=='gcn' else GINConv(mlp,train_eps=True) if kind=='gin' else GINEConv(mlp,edge_dim=edges,train_eps=True))
        self.head=nn.Sequential(nn.Linear(hidden*2+1,hidden),nn.ReLU(),nn.Dropout(.1),nn.Linear(hidden,1))
    def forward(self,b):
        x=b.x
        edge=b.edge_attr
        if self.ablate:
            x=x.clone();x[:,13:]=0;edge=torch.ones_like(edge)
        x=torch.relu(self.input(x))
        for layer in self.layers:
            x=torch.relu(layer(x,b.edge_index,edge) if self.kind=='gine' else layer(x,b.edge_index))
        h=torch.cat([global_mean_pool(x,b.batch),global_add_pool(x,b.batch)/20,b.dose.reshape(-1,1)],dim=1)
        return self.head(h).flatten()


class MLP(nn.Module):
    def __init__(self,features,hidden):
        super().__init__();self.net=nn.Sequential(nn.Linear(features,hidden),nn.ReLU(),nn.Dropout(.1),nn.Linear(hidden,hidden//2),nn.ReLU(),nn.Linear(hidden//2,1))
    def forward(self,x):return self.net(x).flatten()


def fit_neural(kind,graphs,x,table,indices,params,seed,cfg,checkpoint,ablate=False):
    seed_all(seed);torch.set_num_threads(cfg['threads'])
    tr,va=indices['train'],indices['validation'];target=table.y.to_numpy(float)
    counts=table.iloc[tr].family.value_counts();w=np.array([1/counts[f] for f in table.iloc[tr].family]);w/=w.mean()
    ymean=float(np.average(target[tr],weights=w));yscale=max(float(np.sqrt(np.average((target[tr]-ymean)**2,weights=w))),1e-6)
    dose=table.log_dose.to_numpy(float);dm=float(dose[tr].mean());ds=max(float(dose[tr].std()),1e-6)
    items=[]
    for i,g in enumerate(graphs):
        g=g.clone();g.dose=torch.tensor([(dose[i]-dm)/ds],dtype=torch.float32);g.y=torch.tensor([(target[i]-ymean)/yscale],dtype=torch.float32);g.w=torch.ones(1);items.append(g)
    for i,weight in zip(tr,w):items[i].w=torch.tensor([weight],dtype=torch.float32)
    if kind=='mlp':
        model=MLP(x.shape[1],params['hidden']);features=torch.tensor(x,dtype=torch.float32)
        trainloader=None
    else:
        model=GraphModel(kind,graphs[0].num_node_features,graphs[0].num_edge_features,params['hidden'],ablate)
        trainloader=DataLoader([items[i] for i in tr],batch_size=cfg['batch_size'],shuffle=True,num_workers=0)
    optim=torch.optim.AdamW(model.parameters(),lr=params['lr'],weight_decay=.001)
    def infer(ids):
        model.eval()
        with torch.no_grad():
            if kind=='mlp':out=model(features[ids]).numpy()
            else:out=np.concatenate([model(b).numpy() for b in DataLoader([items[i] for i in ids],batch_size=128,shuffle=False,num_workers=0)])
        return out*yscale+ymean
    def val_error(pred):
        import pandas as pd
        return float(pd.DataFrame({'family':table.iloc[va].family.to_numpy(),'error':abs(target[va]-pred)}).groupby('family').error.mean().mean())
    best=float('inf');best_state=None;wait=0;history=[];best_epoch=0
    for epoch in range(cfg['max_epochs']):
        model.train()
        if kind=='mlp':
            ids=np.random.permutation(len(tr))
            for start in range(0,len(tr),cfg['batch_size']):
                loc=ids[start:start+cfg['batch_size']];idx=tr[loc]
                pred=model(features[idx]);y=torch.tensor((target[idx]-ymean)/yscale,dtype=torch.float32)
                loss=((pred-y)**2*torch.tensor(w[loc],dtype=torch.float32)).mean()
                optim.zero_grad();loss.backward();nn.utils.clip_grad_norm_(model.parameters(),5);optim.step()
        else:
            for b in trainloader:
                pred=model(b);loss=((pred-b.y)**2*b.w).mean()
                optim.zero_grad();loss.backward();nn.utils.clip_grad_norm_(model.parameters(),5);optim.step()
        score=val_error(infer(va));history.append(dict(epoch=epoch+1,validation_mae_log10=score))
        if score<best-1e-6:
            best=score;best_state=copy.deepcopy(model.state_dict());best_epoch=epoch+1;wait=0
        else:wait+=1
        if wait>=cfg['patience']:break
    if best_state is None:raise ValueError('Neural training produced no finite validation score')
    model.load_state_dict(best_state)
    torch.save(dict(state_dict=best_state,kind=kind,params=params,ymean=ymean,yscale=yscale,dose_mean=dm,dose_scale=ds,ablate=ablate),checkpoint)
    predictions={role:infer(ids) for role,ids in indices.items() if role!='train'}
    return predictions,dict(validation_mae_log10=best,epochs=len(history),best_epoch=best_epoch,parameters=sum(p.numel() for p in model.parameters()),history=history)
