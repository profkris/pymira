import pickle
import os 
import numpy as np
from matplotlib import pyplot as plt

dirs = '/mnt/data2/Sahai/export/nifti'
ko_prefs = ['KO9','KO7']
ctrl_prefs = ['CTRL4','CTL10']
het_prefs = ['HET7']

def load_data(prefs,pname='radii'):
    param = []
    for pref in prefs:
        path = os.path.join(dirs,pref)
        with open(os.path.join(path,'{}.p'.format(pname)),'rb') as fo:
            param.extend(pickle.load(fo))
    return np.asarray(param)

nbins = 10
pname, param_label,range = 'radii','Vessel radius (um)',[0.,80.]
#pname, param_label,range = 'vessel_length','Vessel length (um)',[0.,200.]
#pname, param_label,range,nbins = 'branching_angle','Branching angle',[0.,180.],20
#pname, param_label,range = 'vessel_volume','Vessel volume (um3)',[0.,5000.]
#pname, param_label,range = 'nconn','Branch node connections',[2.,5.]
ko_param = load_data(ko_prefs,pname=pname)
ctrl_param = load_data(ctrl_prefs,pname=pname)
het_param = load_data(het_prefs,pname=pname)

def histogram(v, nbins=30,range=None,xlabel=None,color='green',clear=True,labels=None,density=True):
    if clear:
        plt.clf()
    #range = [v.min(),v.max()]
        
    # the histogram of the data
    n, bins, patches = plt.hist(v, nbins, density=density,label=labels,range=range)#, facecolor=color, alpha=0.5)
    plt.legend(prop={'size': 10})
    if xlabel is not None:
        plt.xlabel(xlabel)
        
    #plt.legend(handles, labels)

#histogram(ko_param, nbins=30,range=None,color='red',clear=True)
#histogram(ctrl_param, nbins=30,range=None,color='red',clear=False)        
#histogram(het_param, nbins=30,range=None,color='blue',clear=False)  
labels = ['KO','CTRL','HET']              
histogram([ko_param,ctrl_param,het_param], nbins=nbins,range=range,color='red',clear=True,labels=labels,xlabel=param_label,density=False)
plt.show()
    
tots = np.asarray([np.sum(x) for x in [ko_param,ctrl_param,het_param]])
nseg = np.asarray([x.shape[0] for x in [ko_param,ctrl_param,het_param]])
    
import pdb
pdb.set_trace()
