import numpy as np

def linear_scaling(doses, exposures, target_dose):
    doses=np.asarray(doses,float); exposures=np.asarray(exposures,float)
    if len(doses)==0 or np.any(doses<=0) or np.any(exposures<=0) or target_dose<=doses.max():
        raise ValueError('Positive support and strictly extrapolative target required')
    # Log-space intercept constrained to slope one; every support dose contributes.
    return float(target_dose*np.exp(np.mean(np.log(exposures/doses))))

def power_law(doses,exposures,target_dose):
    if len(set(doses))<2: raise ValueError('Power law requires two distinct support doses')
    linear_scaling(doses,exposures,target_dose)  # domain validation
    slope,intercept=np.polyfit(np.log(doses),np.log(exposures),1)
    return float(np.exp(intercept+slope*np.log(target_dose)))
