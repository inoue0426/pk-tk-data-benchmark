"""Positive empirical exposure models; all inputs are permitted support observations."""
import numpy as np
from scipy.optimize import least_squares

class NotIdentifiable(ValueError):
    pass

def validate(doses,exposures,target):
    d=np.asarray(doses,float);y=np.asarray(exposures,float)
    if d.ndim!=1 or d.shape!=y.shape or not len(d) or len(set(d))!=len(d):
        raise ValueError('Distinct paired support doses required')
    if not np.isfinite(np.r_[d,y,target]).all() or np.any(d<=0) or np.any(y<=0) or target<=0:
        raise ValueError('Positive finite dose/exposure required')
    if target in d:raise ValueError('Target dose must be held out')
    return d,y

def exponent(doses,exposures):
    if len(set(doses))<2:raise NotIdentifiable('At least two support doses required')
    return float(np.polyfit(np.log(doses),np.log(exposures),1)[0])

def predict(model,doses,exposures,target,pooled=None,weight=.5,max_condition=1000,minimum_support=3):
    d,y=validate(doses,exposures,target);x=np.log(d);z=np.log(y)
    params={}
    if model=='linear':b=1.
    elif model=='power':b=exponent(d,y)
    elif model=='pooled_power':
        if pooled is None:raise NotIdentifiable('Insufficient disjoint reference components')
        b=float(pooled)
    elif model=='shrinkage_power':
        if pooled is None:raise NotIdentifiable('Insufficient disjoint reference components')
        b=weight*exponent(d,y)+(1-weight)*pooled
    elif model in ['saturating','superlinear']:
        if minimum_support<3:raise ValueError('Curvature support threshold must be at least three')
        if len(d)<minimum_support:raise NotIdentifiable('Curvature requires at least three support doses')
        sign=-1 if model=='saturating' else 1
        # Positive two-parameter rational/superlinear curves, fitted on log exposure.
        def residual(p):return p[0]+x+sign*np.logaddexp(0,x-p[1])-z
        lo=np.log(d.min()/100);hi=np.log(d.max()*100)
        starts=[np.log(d.min()),np.mean(x),np.log(d.max())]
        fits=[least_squares(residual,[np.mean(z-x),k],bounds=([-np.inf,lo],[np.inf,hi]),max_nfev=2000) for k in starts]
        fit=min(fits,key=lambda f:float(np.sum(f.fun**2)))
        k=float(np.exp(fit.x[1]));condition=float(np.linalg.cond(fit.jac))
        if not fit.success or not d.min()<k<d.max() or condition>max_condition:
            raise NotIdentifiable('Curvature scale not supported within observed support range')
        value=float(np.exp(fit.x[0]+np.log(target)+sign*np.logaddexp(0,np.log(target)-fit.x[1])))
        params=dict(a=float(np.exp(fit.x[0])),k=k,jacobian_condition=condition,support_log_sse=float(np.sum(fit.fun**2)))
        if not np.isfinite(value) or value<=0:raise NotIdentifiable('Invalid numerical prediction')
        return value,params
    else:raise ValueError('Unknown model')
    a=float(np.mean(z-b*x));value=float(np.exp(a+b*np.log(target)))
    if not np.isfinite(value) or value<=0:raise NotIdentifiable('Invalid numerical prediction')
    return value,dict(log_a=a,b=float(b))
