import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

##This is taken verbatim from the chong data, going to check today to see if there's any abnormalities due to perhaps
##the response code being different. The relevant functions here are err_plots and its depenencies.

def quantile_list(data,nquantiles):
  """Function that finds the value of n quantiles in an array. \n
  data: 1-D array \n
  nquantiles: number of quantiles to return"""
  quantiles_list = []
  for i in range(nquantiles):
    q = i/nquantiles
    if len(data) == 0:
      data = [0,0]
    quantiles_list.append(np.quantile(data,q))
  return np.array(quantiles_list)

def quantile_data(data,nquantiles):
    """Function that returns a list of the data split into quantiles. \n
  data: 1-D array \n
  nquantiles: number of quantiles to split the data into"""
    data_list = []
    l = np.linspace(0,1,nquantiles+1)
    l = np.quantile(data.rxtime,l)
    for i in range(nquantiles):
        #print((data.rt>=l[i])*(data.rt<l[i+1]))
        data_list.append(data[(data.rxtime>=l[i])*(data.rxtime<l[i+1])])
    return data_list

def rate_by_quantile(data,response,nquants):
  """Returns the rate of occurence for a certain response within each n quantile(s) \n
  data: dataframe of hddm-ready data \n
  response: The specifc response to provide the rate of occurence for \n
  nquants: number of quantiles to return rate for"""
  subs = np.unique(data.subj)
  out = np.zeros((len(subs),nquants))
  for i in range(len(subs)):
      d = data[data.subj==subs[i]]
      qd = quantile_data(d,nquants)
      for j in range(nquants):
          out[i,j] = np.mean(qd[j].response==response)
  print(np.mean(out,axis = 0))
  return np.mean(out,axis=0)

def load_posterior_predictive(model, task = 0, coherence = 0, group = 0):
  """Loads the posterior predictive of a model with the specified paramters"""
  files = os.listdir("data/posterior_predictives/{}".format(model))
  isempty = True
  for file in files:
    if (model in file) and ("task_{}".format(task) in file) and ("coh_{}".format(coherence) in file) and ("group_{}".format(group) in file):
      data = pd.read_pickle("data/posterior_predictives/{}/{}".format(model,file))
      isempty = False
  if isempty:
    print("model with these arguments does not exist")
  else:
    return data


### Plotting functions ###

def err_plot(data,nquants):
  """Plots the rate of high and low dimension errors in each quantile"""
  plt.plot(rate_by_quantile(data,1,nquants) + rate_by_quantile(data,0,nquants),'-o')
  plt.plot(rate_by_quantile(data,2,nquants) + rate_by_quantile(data,0,nquants),'-o')
  plt.ylim((0,0.3))

def single_qpp(data_table, nquants, coh_level='highDim'):
  """Plots a single quantile probability plot. If coherence level is 'highDim', 
  then errors in the high dimension are used for the error rates. Otherwise, low dimension 
  errors are used."""
  if coh_level == 'highDim':
    for i in [-1,1]:
        data = data_table[data_table.highDimCoh == i]
        rate = np.mean(data.response==1)
        print(rate)
        plt.plot(np.repeat(rate,nquants),quantile_list(data.rxtime[data.response==1],nquants),'-o')

        rate = sum(data.response==3)/len(data.rxtime)
        plt.plot(np.repeat(rate,nquants),quantile_list(data.rxtime[data.response==3],nquants),'-o')

  else:
    for i in [-1,1]:
        data = data_table[data_table.lowDimCoh == i]
        rate = sum(data.response==2)
        plt.plot(np.repeat(rate,nquants),quantile_list(data.rxtime[data.response==2],nquants),'-o')

        rate = sum(data.response==3)/len(data.rxtime)
        plt.plot(np.repeat(rate,nquants),quantile_list(data.rxtime[data.response==3],nquants),'-o')


def single_coh_hist(data,nquants):
  """Plots a histogram of the data and idenifies n quantiles in the plot."""
  plt.hist(data.rxtime,bins=50,density=True,range=(0,8))

  q_list = quantile_list(data.rxtime,nquants)
  plt.vlines(q_list,0,1,linestyles='--',colors='r')

  lines = q_list
  lines = np.append(lines,5)

  #locs = []
  for i in range(nquants):
    #locs.append()
    plt.text(np.mean([lines[i],lines[i+1]])-0.1,0.95,'{}'.format(i+1))

    

def descriptive_plts(chong_data):
  """Function plotting plots that only require the data!"""
    
def all_plots(chong_data):
  """Function that plots all plots relevant to chong data study"""
  plt.subplot(2,5,1)
  single_coh_hist(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==-1)], 5)
  plt.title('Reaction Time Quantiles \n by Coherence: LL',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  plt.subplot(2,5,2)
  plt.title('Reaction Time Quantiles \n by Coherence: LH',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  single_coh_hist(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==1)],5)
  plt.subplot(2,5,3)
  plt.title('Reaction Time Quantiles \n by Coherence: HL',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  single_coh_hist(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==-1)],5)
  plt.subplot(2,5,4)
  plt.title('Reaction Time Quantiles \n by Coherence: HH',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  single_coh_hist(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==1)],5)

  plt.subplot(2,5,6)
  err_plot(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==-1)],5)
  plt.title('Error Rate in RT Quantiles \n by Coherence: LL',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['High Dimension Error','Low Dimension Error'])

  plt.subplot(2,5,7)
  err_plot(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==1)],5)
  plt.title('Error Rate in RT Quantiles \n by Coherence: LH',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['High Dimension Error','Low Dimension Error'])


  plt.subplot(2,5,8)
  err_plot(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==-1)],5)
  plt.title('Error Rate in RT Quantiles \n by Coherence: HL',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['High Dimension Error','Low Dimension Error'])


  plt.subplot(2,5,9)
  err_plot(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==1)],5)
  plt.title('Error Rate in RT Quantiles \n by Coherence: HH',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['High Dimension Error','Low Dimension Error'])


  plt.subplot(2,5,5)
  plt.title('High Dimension Quantile\nProbability Plot', fontsize = 18)
  single_qpp(chong_data,5)
  plt.legend(['Low coherence,\n incorrect','Low coherence,\n correct','High coherence,\n incorrect','High coherence,\n correct'])
  plt.xlabel('Response Probability',fontsize=14)
  plt.ylabel('Reaction Time (seconds)',fontsize=14)

  plt.subplot(2,5,10)
  plt.title('Low Dimension Quantile\nProbability Plot', fontsize = 18)
  plt.ylabel('Reaction Time (seconds)')
  single_qpp(chong_data,5,'lowDim')
  plt.legend(['Low coherence,\n incorrect','Low coherence,\n correct','High coherence,\n incorrect','High coherence,\n correct'])
  plt.xlabel('Response Probability',fontsize=14)
  plt.ylabel('Reaction Time (seconds)',fontsize=14)
    
    

