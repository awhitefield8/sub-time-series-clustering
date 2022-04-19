from sklearn.cluster import KMeans, AgglomerativeClustering
import numpy as np
import math
from queue import PriorityQueue


from models import *
import copy




def naive_0(PANEL,
            k,
            c1=1,
            c2=1,
            c3=1,
            c4=0):
    """ perform k-means++ with an L1 distance
    Args:
        trajs: trajectories list
    Returns:
        sub-time series clustering object
    """
    X = np.array(PANEL.raw_trajs)

    kmeans_res = KMeans(n_clusters=k,
                init='k-means++',
                n_init=10,
                max_iter=300).fit(X)
    centers = kmeans_res.cluster_centers_
    labels = kmeans_res.labels_

    # update centers to form: [pth1,pth2,...]
    centers_proc = []
    counter = 0
    for c in centers: #generate pathlet list
        trajIDs_iter = [i for i in range(len(labels)) if labels[i] == counter]
        pth_iter = pathlet2(bounds = (0,len(c)) ,pt_array = c, id=counter,trajIDs=trajIDs_iter)
        centers_proc.append(pth_iter)
        counter = counter + 1

    # update lables dict of form [trajID: [pth1,pth2,...]]
    labels_proc = {}
    for i in range(len(labels)):
        pth = centers_proc[labels[i]]
        labels_proc[i] = [pth]

    stc_iter = Stc(c1=c1,c2=c2,c3=c3,c4=c4,panel=PANEL,pathlets=centers_proc,assignments=labels_proc) 

    stc_iter.genSuperPathlets() #generate superpathlets

    return(stc_iter)

    


def naive_1(PANEL,
            k_list,
            c1=1,
            c2=1,
            c3=1,
            c4=0):
    """ perform k-means++ with an L1 distance on various values of k, and pick the value that minimises the objective function
    Args:
        trajs: trajectories list
        k_list: list of ks to try
    Returns:
        sub-time series clustering object
    """

    current_loss = math.inf
    best_stc = []

    for k in k_list:
        
        stc_iter = naive_0(copy.deepcopy(PANEL),
                           k=k,
                           c1=c1,
                           c2=c2,
                           c3=c3,
                           c4=c4) 
        
        if stc_iter.loss() < current_loss:
            current_loss = stc_iter.loss() 
            best_stc = stc_iter
    
    return(best_stc)




def naive_2(PANEL,
            k_list,
            c1=1,
            c2=1,
            c3=1,
            c4=0):
    """ Naive 1, but then greedily drop trajectories from cluster assignments.
    Args:
        trajs: trajectories list
        k_list: list of ks to try
    Returns:
        sub-time series clustering object
    """

    n1_res = naive_1(PANEL=copy.deepcopy(PANEL),
                     k_list=k_list,
                     c1=c1,
                     c2=c2,
                     c3=c3,
                     c4=c4)
    
    current_res = n1_res
    current_loss = n1_res.loss()

    
    # loop through each trajectory, and check marginal benefit of dropping (if benefit, then drop)

    for traj_id in range(len(n1_res.panel.trajs)): 
        new_res = copy.deepcopy(current_res)
        new_res.updateSingleAssignment(updated_assignment=[],traj_id=traj_id)
        new_loss = new_res.loss()
        #print("new loss: ", str(new_loss))
        #print("current loss: ", str(current_loss))
        if new_loss < current_loss:
            #print("successful drop!")
            current_res = new_res
            current_loss = new_loss

    return(current_res)
        #turn into stc_obj


def naive_3(PANEL,
            k_list,
            seg_list,
            c1=1,
            c2=1,
            c3=1,
            c4=0,
            greedy_in_period = False): #if true, only apply c1 within periods, not that the end
    """
    Naive 3 - clusters different segments. Does not mwerge across segments

    Args:
        trajs: trajectories list
        k_list: list of ks to try
        seg_list: list of segments (e.g. [1,2,4], where 4 would divide into 4 regions)
    Returns:
        sub-time series clustering object
    """
    init_res = naive_1(PANEL=copy.deepcopy(PANEL),
                     k_list=k_list,
                     c1=c1,
                     c2=c2,
                     c3=c3)
    best_stc_ = init_res
    if greedy_in_period == True:
        best_loss_ = math.inf
    else:
        best_loss_ = init_res.loss()

    trajs_raw = PANEL.raw_trajs

    traj_ids = [ traj.trajID for traj in PANEL.trajs]

    for seg in seg_list: #for each segment number option
        #print("starting segments count: " + str(seg))
        seg_legth_iter = math.floor(PANEL.T/seg)
        #initialse pathlet assignments 
        pathlet_seg_iter = []
        assign_seg_iter = dict(zip(traj_ids,[[]]*len(traj_ids)))

        for i in range(seg):
            if i==(seg-1):
                bound_iter = (i*seg_legth_iter,PANEL.T  )
            else:
                bound_iter = (i*seg_legth_iter,(i*seg_legth_iter)+seg_legth_iter )
            #print("seg: " + str(i))
            #print("bounds: %d %d" % (bound_iter[0], bound_iter[1] ))

            trajs_raw_iter = [ traj[ bound_iter[0]:bound_iter[1]  ] for traj in trajs_raw ]

            panel_iter = Panel(trajObj_list = trajs_raw_iter,
                               start_time=bound_iter[0])

            res_iter = naive_1(PANEL=panel_iter,
                                k_list=k_list,
                                c1=c1,
                                c2=c2,
                                c3=c3,
                                c4=c4)
            
            for pth in res_iter.pathlets:
                pth.shiftBounds(bound_iter[0])

            pathlet_seg_iter.extend(res_iter.pathlets)
            for k ,v in res_iter.assignments.items():
                assign_seg_iter[k] = assign_seg_iter[k] + res_iter.assignments[k]
        

        #calculate assignment for iteration

        if greedy_in_period == True: #if greedy in period is true, we don't want to double penalise lots of pathlets across rounds...
        # ... as it is already taken into account by naive 1 in period. It will be set to false if naive 3 is called alone 
            c4_aug = 0
        else:
            c4_aug=c4

        stc_iter = Stc(c1=c1,c2=c2,c3=c3,c4=c4_aug,panel=PANEL,pathlets=pathlet_seg_iter,assignments=assign_seg_iter) 
        loss_iter_ = stc_iter.loss()
        #print(loss_iter_)
        if loss_iter_ < best_loss_:
            #print("improvement to %d" %  (loss_iter_))
            best_stc_ = stc_iter
            best_loss_ = loss_iter_

    best_stc_.genSuperPathlets()

    return(best_stc_)
        #turn into stc_obj
        





def greedy1(PANEL,
            k_list,
            seg_list,
            c1=1,
            c2=1,
            c3=1,
            c4=0):
    """
    Naive 3, with greedy merges of pathlets

    Only works for single input of K

    Args:
        trajs: trajectories list
        k_list: list of ks to try
        seg_list: list of segments (e.g. [1,2,4], where 4 would divide into 4 regions)
    Returns:
        sub-time series clustering object
    """
    if len(k_list) > 1:
        return("please specify a single number of cluster given a period in time")

    res = naive_3(copy.deepcopy(PANEL),
            k_list, #note that this will only output a single k per round
            seg_list,
            c1=0,
            c2=c2,
            c3=c3,
            c4=0) #run naive 3, but set c1 (number of paths) and c4 (switches) to zero

    res.updateParams(c1_=c1,c2_=c2,c3_=c3,c4_=c4) # update parameters 

    pathlet_bounds = copy.deepcopy(list(set([ i.bounds for i in res.pathlets]))) #list of bounds for pathlets to consider, deep copy and bounds may change
    pathlet_bounds.sort() # sort list
    n_segs = len(pathlet_bounds)

    for seg_id in range(n_segs):
        if seg_id==n_segs-1: #nothing to do at the final segment
            pass
        else: 
            ### for each break - 2 steps
            # 1) pick which clusters we merge
            # 2) allow trajs to switch clusters
            # 3) update res object
            
            bounds_start = pathlet_bounds[seg_id]
            bounds_end = pathlet_bounds[seg_id+1]
            meet_point = bounds_start[1]
            if meet_point != bounds_end[0]:
                print("error, not valid breakpoint")
            print("bounds: (%d,%d)" % (bounds_start[0],bounds_start[1]))
            pathlets_start = [ i for i in res.pathlets if i.bounds[1] == meet_point ]
            pathlets_end = [ i for i in res.pathlets if i.bounds[0] == meet_point ]

            #for each pair of clusters between segments, find the similarity
            sim_matrix = np.array( [[1.0] * len(pathlets_end) for i in range(len(pathlets_start))] )
            for s in range(len(pathlets_start)):
                for e in range(len(pathlets_end)):
                    sim_matrix[s,e] = round(comp_diff(pathlets_start[s].trajIDs ,pathlets_end[e].trajIDs),2)

            #add clusters to priority queue in order how good of the match they have
            matches = list()
            for i in range(len(pathlets_start)):
                closest_match_score = copy.deepcopy(min(sim_matrix[i]))
                closest_match_id = sim_matrix[i].argmin()
                if closest_match_score<0.90: #only merge if close enough match <<<<<< this is an important parameter
                    #print("sucess: match score: " + str(closest_match_score))
                    matches.append( (i,closest_match_id) )
                    sim_matrix[:,closest_match_id] = 1   #add weights to column just used
                    #print(sim_matrix)
            
            #merge clusters if the match is good enough
            #print(matches)
            if len(matches) == 0 : 
                pass #if no matches, leave the pathlets as is
            else:
                PTH_TUPLE_LIST = list()
                for m in matches:
                    PTH_TUPLE_LIST.append( (pathlets_start[m[0]],pathlets_end[m[1]]) )
                #print(PTH_TUPLE_LIST)
                res.mergePathlets(pth_tuple_list=PTH_TUPLE_LIST,bounds=bounds_end)
            
        #finally, generate superpaths
        res.genSuperPathlets()
                
    return(res)


    






def greedy2(PANEL,
            k_list,
            seg_list,
            c1=1,
            c2=1,
            c3=1,
            c4=0):
    """ Greedy 1, but with multiple splits
    Works for multiple inputs of K

    Args:
        trajs: trajectories list
        k_list: list of ks to try
        seg_list: list of segments (e.g. [1,2,4], where 4 would divide into 4 regions)
    Returns:
        sub-time series clustering object
    """

    # we want to greedily split the space in each segment - call naive 3 with c4=0. c1 is ok to keep, as we have greedy_in_period set to true
    res = naive_3(copy.deepcopy(PANEL),
            k_list, #allows multiple k
            seg_list,
            c1=c1,
            c2=c2,
            c3=c3,
            c4=0,
            greedy_in_period = True) #run naive 3

    res.updateParams(c1_=c1,c2_=c2,c3_=c3,c4_=c4) # update parameters 

    pathlet_bounds = copy.deepcopy(list(set([ i.bounds for i in res.pathlets]))) #list of bounds for pathlets to consider, deep copy and bounds may change
    pathlet_bounds.sort() # sort list
    n_segs = len(pathlet_bounds)

    for seg_id in range(n_segs):
        if seg_id==n_segs-1: #nothing to do at the final segment
            pass
        else: 
            ### for each break - 2 steps
            # 1) pick which clusters we merge
            # 2) allow trajs to switch clusters
            # 3) update res object
            
            bounds_start = pathlet_bounds[seg_id]
            bounds_end = pathlet_bounds[seg_id+1]
            meet_point = bounds_start[1]
            if meet_point != bounds_end[0]:
                print("error, not valid breakpoint")
            print("bounds: (%d,%d)" % (bounds_start[0],bounds_start[1]))
            pathlets_start = [ i for i in res.pathlets if i.bounds[1] == meet_point ]
            pathlets_end = [ i for i in res.pathlets if i.bounds[0] == meet_point ]

            #for each pair of clusters between segments, find the similarity

            #a) calculate goodness of each new path
            #b) if goodness is good enough - add it 
            #c) next recalc goodness other new paths - if good enough, add. Continue
            #d) add new paths to res. may want to kill existing paths

            #e) new error term: total_path costs
            
            
            
            
        #finally, generate superpaths
        res.genSuperPathlets()
                
    return(res)






def greedy3(PANEL,
            k_list,
            seg_list,
            c1=1,
            c2=1,
            c3=1,
            c4=0):
    """ Greedy 1, but with multiple splits
    Works for multiple inputs of K
    And iterates for different seg splits

    Args:
        trajs: trajectories list
        k_list: list of ks to try
        seg_list: list of segments (e.g. [1,2,4], where 4 would divide into 4 regions)
    Returns:
        sub-time series clustering object
    """


    pass


### things to add
# option to recalibrate greedy centers (this won't actually matter for grouping, but will impact objective)
