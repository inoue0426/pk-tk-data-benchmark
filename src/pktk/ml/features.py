"""Cached fixed molecular features; no learned preprocessing or labels in cache."""
import json
import numpy as np
import torch
from rdkit import Chem,DataStructs
from rdkit.Chem import Descriptors,Crippen,Lipinski,rdMolDescriptors,rdFingerprintGenerator
from torch_geometric.data import Data
from .data import ROOT,OUT,sha,dump

DESCRIPTORS=[('MW',Descriptors.MolWt),('logP',Crippen.MolLogP),('TPSA',rdMolDescriptors.CalcTPSA),
 ('HBD',Lipinski.NumHDonors),('HBA',Lipinski.NumHAcceptors),('rotatable',Lipinski.NumRotatableBonds),
 ('rings',Lipinski.RingCount),('aromatic_rings',Lipinski.NumAromaticRings),('fraction_sp3',rdMolDescriptors.CalcFractionCSP3),
 ('heavy_atoms',Lipinski.HeavyAtomCount),('NHOH',Lipinski.NHOHCount),('NO',Lipinski.NOCount),
 ('MR',Crippen.MolMR),('heteroatoms',Lipinski.NumHeteroatoms),('valence_electrons',Descriptors.NumValenceElectrons),
 ('formal_charge',Chem.GetFormalCharge),('aliphatic_rings',Lipinski.NumAliphaticRings)]

def onehot(x,values):return [float(x==v) for v in values]+[float(x not in values)]

def graph(smiles):
    m=Chem.MolFromSmiles(smiles);Chem.AssignStereochemistry(m,cleanIt=True,force=True)
    atoms=[]
    for a in m.GetAtoms():
        atoms.append(onehot(a.GetAtomicNum(),[1,5,6,7,8,9,14,15,16,17,35,53])+onehot(a.GetDegree(),list(range(7)))+
          onehot(a.GetFormalCharge(),[-3,-2,-1,0,1,2,3])+onehot(str(a.GetHybridization()),['SP','SP2','SP3','SP3D','SP3D2'])+
          onehot(a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else 'none',['R','S','none'])+
          [float(a.GetIsAromatic()),float(a.IsInRing()),a.GetTotalNumHs()/4,a.GetMass()/100])
    edges=[];bonds=[]
    for b in m.GetBonds():
        v=onehot(str(b.GetBondType()),['SINGLE','DOUBLE','TRIPLE','AROMATIC'])+onehot(str(b.GetStereo()),['STEREONONE','STEREOZ','STEREOE','STEREOCIS','STEREOTRANS'])+[float(b.GetIsConjugated()),float(b.IsInRing())]
        edges.extend([[b.GetBeginAtomIdx(),b.GetEndAtomIdx()],[b.GetEndAtomIdx(),b.GetBeginAtomIdx()]]);bonds.extend([v,v])
    return Data(x=torch.tensor(atoms,dtype=torch.float32),edge_index=torch.tensor(edges,dtype=torch.long).T.contiguous().reshape(2,-1),edge_attr=torch.tensor(bonds,dtype=torch.float32).reshape(-1,13))


def load_features(table):
    cache=ROOT/'results/ml/cache';cache.mkdir(parents=True,exist_ok=True)
    key=sha(OUT/'tasks.csv')+'-'+sha(__file__)
    if (cache/'manifest.json').exists() and json.loads((cache/'manifest.json').read_text())['key']==key:
        a=np.load(cache/'features.npz');g=torch.load(cache/'graphs.pt',weights_only=False)
        return a['descriptors'],a['fingerprints'],g
    generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
    descriptors=[];fps=[];graphs=[]
    for r in table.itertuples():
        m=Chem.MolFromSmiles(r.canonical_smiles)
        descriptors.append([fn(m) for _,fn in DESCRIPTORS])
        fps.append(generator.GetFingerprintAsNumPy(m));graphs.append(graph(r.canonical_smiles))
    desc=np.asarray(descriptors,dtype=np.float64);desc[~np.isfinite(desc)]=np.nan
    fps=np.asarray(fps,dtype=np.float32)
    np.savez_compressed(cache/'features.npz',descriptors=desc,fingerprints=fps)
    torch.save(graphs,cache/'graphs.pt');dump(cache/'manifest.json',dict(key=key,descriptor_names=[n for n,_ in DESCRIPTORS],fp='Morgan radius2 2048bits with chirality',atom_features=graphs[0].num_node_features,bond_features=graphs[0].num_edge_features))
    return desc,fps,graphs


def nearest_similarity(fps,train,query):
    a=fps[query];b=fps[train];intersection=a@b.T
    union=a.sum(1)[:,None]+b.sum(1)[None,:]-intersection
    return (intersection/np.maximum(union,1)).max(axis=1)
