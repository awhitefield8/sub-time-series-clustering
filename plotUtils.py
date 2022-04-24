import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from random import sample
from sklearn.cluster import KMeans, AgglomerativeClustering
import numpy as np

from clusteringAlgos import *
from models import *
import pandas as pd


### plotting for stc objects

def pth_to_df(pth_list):
    """ turn trajectories into plot"""
    df_output = pd.DataFrame(columns = ['new_id', 'date_id', 'value','id','cluster'])
    
    counter = 1
    pth_id = 1
    for pth in pth_list:
        for time_id in range(len(pth.pt_array)):
            c = 1 #add later
            startdate = pth.bounds[0]
            new_row = {'new_id':counter, 'date_id':time_id + startdate , 'value': pth.pt_array[time_id],'id': pth_id,'cluster':c}
            df_output = df_output.append(new_row, ignore_index = True)
            counter = counter + 1  
        pth_id = pth_id + 1
    return(df_output)
            


def trajs_to_df(traj_obj_list):
    """ turn trajectories into plot"""
    df_output = pd.DataFrame(columns = ['new_id', 'date_id', 'value','id','cluster'])
    
    counter = 1
    for traj in traj_obj_list:
        for time_id in range(len(traj.pt_array)):
            c = 1 #add later
            new_row = {'new_id':counter, 'date_id':time_id , 'value': traj.pt_array[time_id],'id': traj.trajID,'cluster':c}
            df_output = df_output.append(new_row, ignore_index = True)
            counter = counter + 1  
    return(df_output)
            
        

def plot_stc(stc_obj):
    """
    plots an stc object
    """
    
    trajs = stc_obj.panel.trajs #list of trajectory objects
    pathlets = stc_obj.superPths #list of superPths
    
    # add dfs
    traj_df = trajs_to_df(trajs)
    pth_df = pth_to_df(pathlets)
    
    # identify breaks
    breaks1 = list(set( [pth.bounds[1] for pth in stc_obj.superPths]))
    breaks2 = list(set( [pth.bounds[0] for pth in stc_obj.superPths]))
    breaks = breaks1 + breaks2
    return({"traj_df":traj_df,
            "pth_df":pth_df,
            "breaks": breaks})
    
    

### for dfs


def df_to_plotdf(df,colnames):
    """
    takes in a df and converts to one that is easy to plot
    
    """

    counter = 1
    df_plot = pd.melt(df,
                  id_vars=["date_id"],
                  value_vars=colnames,
                  var_name = "id")
    df_plot["group"] = 1
        
    return(df_plot)
        





### testing
if __name__ == "__main__":
    df_abnormal_returns = pd.read_csv("data/abnormal_returns_data.csv")
    start_date = 3050
    df_preproc = df_abnormal_returns
    df = df_preproc[df_preproc["date_id"] >= start_date]
    company_names = list(df.columns)[2:]
    trajectories = []
    for i in company_names:
        trajectories.append(list(df[i]))
    traj_panel = Panel(trajObj_list = trajectories,start_time=1)

    C1,C2,C3,C4 = 10,1,1,1
    K_LIST = [5]
    SEG_LIST = [4]

    res_stc_g1 = greedy1(PANEL = traj_panel,
                        k_list= K_LIST,
                        seg_list = SEG_LIST, 
                        c1=C1,
                        c2=C2,
                        c3=C3,
                        c4=C4)
    
    plot_dfs = plot_stc(res_stc_g1)
    print(plot_dfs["pth_df"])




### plotting for list of list trajectories


def traj_list_to_df(trajs,start_time=0):
    """ list of list trajectories to df """

    df_output = pd.DataFrame(columns = ['new_id', 'date_id', 'value','id','cluster'])
    
    counter = 1
    pth_id = 1
    for traj in trajs:
        for time_id in range(len(traj)):
            new_row = {'new_id':counter, 'date_id':time_id  + start_time , 'value': traj[time_id],'id': pth_id,'cluster':1}
            df_output = df_output.append(new_row, ignore_index = True)
            counter = counter + 1  
        pth_id = pth_id + 1
    return(df_output)