from genUtils import *
import copy
import math
import itertools

import numpy as np




class pt(object):
    """ A class representing a point in a trajectory.
    
    Attributes:
        lat (float): latitude.
        lon (float): longitude.
        trajID (int): ID of the point's trajectory.
        t (float): timestamp associated with the point.
    """
    
    def __init__(self):
        """ Default constructor with dummy initialization values. """
        self.lat = 0.0
        self.lon = 0.0
        self.trajID = -1
        self.t = -1.0
        
    def __hash__(self):
        """ Computes a hash so that pt objects can be used as dict keys. """
        return hash((self.lat, self.lon, self.trajID, self.t))
        
    def __eq__(self, other):
        """ Define == operator for pt objects. """
        return (self.lat==other.lat and self.lon==other.lon and self.trajID==other.trajID and self.t==other.t)
        
    def __str__(self):
        """ Return string to be output while printing a pt object. """
        return "Point TrajID %d ; lat-long (%f,%f); time %f" % (self.trajID, self.lat, self.lon, self.t)
    


class traj(object):
    """ A class representing a trajectory. 
    
    Attributes:
        pts: list of points in the trajectory.
    """
    
    def __init__(self,trajID,pt_array=None):
        """ Initialize trajectory with empty list of points. """
        self.pts = []
        self.assignments = [] #list of pathlets that cover it
        self.pt_array = pt_array
        self.trajID = trajID
        self.assingmentBounds = [] #a list of bounds corresponding to the assingment  ### { pth1_hash: [(),(),()] }

    def addPt(self, lat, lon, trajID, t):
        """ Add a pt to the trajectory.
        
        Args:
            lat (float): latitude of point.
            lon (float): longitude of point.
            trajID (int): trajID of the point (all points of a traj will have same ID).
        """
        p = pt()
        p.lat = lat
        p.lon = lon
        p.trajID = trajID
        p.t = t
        self.pts.append(p)

    def sortPts(self):
        """ Sort points of trajectory in ascending order of timestamp. """
        self.pts = sorted(self.pts, key = lambda x: x.t)

    def replaceAssignment(self,a):
        """
        replace whole assignment with a new assignment list a in form: [pth1,pth2,...]
        """
        self.assignments = a 

    def removePathlet(self,pth):
        self.assignments.remove(pth)

    def addPathlet(self,pth):
        self.assignments.append(pth)

    def switches(self):
        """ count switches """
        s = 0
        prev_pth = None
        for pth in self.assignments:
            if prev_pth is None:
                pass #no switches if always in a cluster
            else:
                if pth.parent != prev_pth:
                    s = s + 1
            prev_pth = pth
        return(s)

    
    def coverage(self):
        """ counts the proportion of the trajectory that is covered """
        cover = 0
        if len(self.assignments) == 0:
            return(cover)
        else:
            for pth in self.assignments:
                cover = cover + pth.bounds[1] - pth.bounds[0] 

            if cover > len(self.pts):
                print("printing bounds")
                print("len: " + str(len(self.pts)))
                for pth in self.assignments:
                    print(pth)
                    

            return(cover/len(self.pts))



class pathlet2(object):
    """ A class representing a pathlet
        
        Attributes:
            bounds (int,int): start and end indices of points in the timestamp-sorted list of
                              trajectory points.
            pt_array
            id: pathlet id
            trajIDs: list of traj ids in form: [trajID1,trajID2,...]
            parent=None
    """
    
    def __init__(self, bounds, pt_array,id,trajIDs=[],parent=None,child=None):
        """ Initialize pathlet with points from a trajectory.
        
            Args:
                bounds (int,int): start and end indices.
                id: pathlet id
                pt_array: np.array
                trajIDs (int): ID of trajectories that map to it
                parent: parent pathlet
        """
        self.bounds = bounds
        self.pt_array = pt_array
        self.id=id
        self.trajIDs=trajIDs
        self.parent = parent
        self.child = child #currently only single children

        
    def __str__(self):
        """ Return string to be output while printing a pathlet object. """
        return "Pathlet id %d ; bounds (%d, %d); Trajs %d: " % (self.id, self.bounds[0], self.bounds[1],len(self.trajIDs))
        
    def __hash__(self):
        """ Define a hash function so that pathlets can be used as keys in a dict. """
        return hash((self.id, self.bounds[0], self.bounds[1]))

    def updateBounds(self,new_bounds):
        self.bounds = new_bounds

    def shiftBounds(self,shift):
        self.updateBounds( (self.bounds[0] + shift, self.bounds[1] + shift) )
        #self.bounds = (self.bounds[0] + shift, self.bounds[1] + shift)

    def addTrajIds(self,trajIDs):
        self.trajIDs.extend(trajIDs)

    def updatePointArray(self,pt_array):
        """ updates point array, and in turn bounds
        always adds points to the right
        """
        self.pt_array = pt_array
        self.updateBounds(  (self.bounds[0], self.bounds[0] +  len(pt_array)  ) )

    def updateParent(self,new_parent):
        """ add a parent to the pathlet """
        self.parent = self.new_parent

    def addChild(self,new_child):
        """ a child to pathlet """
        self.child = new_child #currently only single children



class superPath(object):
    """ pathlet combining a list of paths that are connected
    """
    def __init__(self, id,pt_array,bounds):
        self.id = id
        self.pt_array = pt_array
        self.bounds = bounds


    

    def extendPath(self,new_pt_array):
        """ extend path and bounds accordingly"""
        self.pt_array = np.concatenate((self.pt_array, new_pt_array), axis=0)
        self.bounds = (self.bounds[0],self.bounds[0] + len(self.pt_array))

    def __str__(self):
        """ Return string to be output while printing a pathlet object. """
        return "Pathlet id %d ; bounds (%d, %d)" % (self.id, self.bounds[0], self.bounds[1])
        


class Panel(object):
    """ A list of time series lists

    Attributes:

    """

    def __init__(self,trajObj_list,start_time=1):
        self.T = len(trajObj_list[0]) #length of first
        self.N = len(trajObj_list) #number of time series objects
        self.start_time = start_time
        
        #add structure to trajectories and points
        trajs = []
        counter_traj_id = 0
        for raw_traj in trajObj_list: #for each time series list
            trajs_it = traj(trajID=counter_traj_id,pt_array = raw_traj)
            counter_time = start_time
            for p in raw_traj:
                trajs_it.addPt(lat = start_time, lon = p, trajID = counter_traj_id, t = counter_time)
                counter_time = counter_time + 1
            
            counter_traj_id = counter_traj_id + 1
            trajs.append(trajs_it)
        
        self.raw_trajs = trajObj_list
        self.trajs = trajs

    def replaceAssignments(self,assignments):
        """updates its trajectories with assignments
        Args:
            assignments: dictionary in form: [trajID: [pth1,pth2,...,pthn], ...]
            """
            
        for k ,v in assignments.items():
            self.trajs[k].replaceAssignment(v)

class Stc(object):
    """ A class representing a sub-time series object
    
    Attributes:
        c1,c2,c3: cost of pathlets, cost of non coverage, cost of distances
        panel: panel object
        pathlets: pathlets list in form [pth1,pth2,...]
        assignments: allocation of trajectories to pathlets in form: [trajID: [pth1,pth2,...,pthn], ...]
        
    """
    
    def __init__(self,c1,c2,c3,c4,panel,pathlets,assignments,metric="l1",superPths=[],c5=0):
        """ Default constructor with dummy initialization values. """
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4
        self.c5 = c5
        
        self.pathlets = pathlets #pathlet list. A list of pathlets [pth1,pth2,...]. Each pathlet contains list of trajs attached to it
        self.assignments = assignments #for each traj - assign a list indicating subtrajectories assigned to each cluster: {trajID: [pth1, pth2, ...]}
        self.metric = metric
        self.superPths = superPths

        panel_new = copy.deepcopy(panel)
        panel_new.replaceAssignments(assignments) #initialise by adding assingments (helpful after kmean operation)
        self.panel = panel_new

        
    def __str__(self):
        """ Return string to be output while printing a pt object. """
        return "SuperPathlets: %a, Total cost: %a" % (self.totalSuperpaths(),self.loss())

    def updateParams(self,c1_,c2_,c3_,c4_,c5_=0):
        """ update model parameters """
        self.c1 = c1_
        self.c2 = c2_
        self.c3 = c3_
        self.c4 = c4_
        self.c5 = c5_

    
    def totalSuperpaths(self):
        """ return total number of superpaths"""
        total_superpaths = [pth for pth in self.pathlets if pth.parent is None] # remove those who have parents
        return(len(total_superpaths))

    def totalPaths(self):
        """ return total number of paths """
        return(len(self.pathlets))

    def switches(self):
        switches = 0
        for traj in self.panel.trajs:
            ### calculate number of pthlets each traj is in (want to punish lots of switches)
            switches = switches + traj.switches()
        return(switches)


    def apply_metric(self,a,b):
        if self.metric == "l1":
            return(l1_dist(a,b))
        else:
            print("unusable metric")

    def loss(self):
        
        ### total unique paths
        total_superpaths = self.totalSuperpaths()

        #frac uncovered points
        frac_uncovered_points = sum([  (1- traj.coverage()) for traj in self.panel.trajs])

        ### quality
        quality = 0
        for traj in self.panel.trajs:
            ### calculate loss for each traj
            for pth in traj.assignments: # iterates through list of form: [pth1,pth2,...]
                bounds = pth.bounds
                quality = quality + self.apply_metric(traj.pt_array[bounds[0]:bounds[1]],pth.pt_array) #add 1 to upper bound

        #switches
        switches = self.switches()

        # total pathlets
        total_paths = self.totalPaths()

        mu = (self.c1*(total_superpaths)) + (self.c2*(frac_uncovered_points)) + (self.c3*(quality)) + self.c4*switches + self.c5*total_paths
        #print("paths: %d"  %(total_paths))
        #print("frac uncovered: %d" % (frac_uncovered_points))
        #print("quality: %d"  %(quality)) 
        #print("switches: %d"  %(switches)) 
        return(mu)

    def updateSingleAssignment(self,updated_assignment,traj_id):
        """updates assignment for a single trajectorty """
        self.panel.trajs[traj_id].replaceAssignment(updated_assignment)

    def addPathlet(self,new_pathlet):
        """ add new pathlet object """
        self.pathlets.append(new_pathlet)

    def mergePathlets(self,pth_tuple_list,bounds):
        """ merges extend pths into new period - used in greedy 1

        Args:
            pth_tuple: list of pairs of paths to merge [(pth1.1,pth2.1),(pth1.2,pth2.2),...] 
            bounds: bounds of the update 

        for each pair
            0) extend points of pth1.1 to include all of pth1.1 and pth2.1. we do this by making new, but linked pathlets
            1) all traj in pth2.1, remove assignment
        
        we have extended pth1.1,pth1.2, ect... We have killed pth2.1,pth2.2 ...
        2)
        next we need to assign all unassigned trajectories in this segment
        for each traj - assign to nearst center remaining
        either assign to old cluster, or to new. assign greedily
        then update all pathlets to be the k-median of thier cluster

         """

        unassigned_trajs= []

        new_pth_id_counter = (self.totalPaths()*bounds[0]) #ensure no counting overlap
        for pth_tuple in pth_tuple_list:
            left_path = pth_tuple[0]
            right_path = pth_tuple[1]
            # 0) extend pth1 by creating linked path
            new_path = pathlet2(bounds = bounds ,pt_array = right_path.pt_array, id=new_pth_id_counter,trajIDs=[],parent=left_path)
            left_path.addChild(new_path) # <<< newline
            self.addPathlet(new_path)
            new_pth_id_counter = new_pth_id_counter + 1 #update counter so no dupe ids

            # 1) remove assignment from pth2
            for traj_id in right_path.trajIDs:
                traj = self.panel.trajs[traj_id] #these match (have confirmed)
                traj.removePathlet(right_path)
                unassigned_trajs.append(traj) #add to unassigned list
 
            # 2) remove from panel objects
            self.pathlets.remove(right_path)
            
        
        for traj in unassigned_trajs:
            #assign to closest trajectory (trading off proximity and previous assignment)
            traj_path = traj.pt_array[bounds[0]:bounds[1]]
            candidate_paths = [ i for i in self.pathlets if i.bounds[1] == bounds[1] ] #all candiate path end at end of the right period

            #print("candidate paths")
            #for i in candidate_paths:
            #    print(i)

            # next find if previous path exists in this period
            existing_assignments = traj.assignments
            prev_assignment = [pth for pth in existing_assignments if pth.bounds[1] == bounds[0] ] #prev assignment
            parent_obj = prev_assignment[0]
            prev_path = [pth for pth in self.pathlets if pth.parent == parent_obj]

            #print("prev paths")
            #for i in prev_path:
            #    print(i)

            if len(prev_path) > 0:
                current_path = prev_path[0]
                pth_path = prev_path[0].pt_array[(bounds[0] - current_path.bounds[0]): (bounds[0] - current_path.bounds[0]) + bounds[1]]
                current_cost = self.apply_metric(traj_path,pth_path) - self.c4 #benefit of no switching cost
            else:
                current_cost = math.inf
                current_path = 1
            for pth in candidate_paths:
                pth_path = pth.pt_array[(bounds[0] - pth.bounds[0]): (bounds[0] - pth.bounds[0]) + bounds[1]]
                cost = self.apply_metric(traj_path,pth_path)
                if cost < current_cost:
                    current_cost = cost
                    current_path = pth
            traj.addPathlet(current_path) #add pathlet to traj
            current_path.addTrajIds([traj.trajID]) #add traj to pathlet


    def mergePathlets2(self,left_paths,right_paths,bounds):
        """ - used in greedy 2
        1) extrapolate olds paths
        2) out of all paths, (old extended ones, and new ones, see which is best to add)
        3) update assingments
        
         """
        #1 and 2) extrapolate old paths, calculate marginal cost of new and old paths
        
        old_extrapolated_path_candidates = []
        marginal_benefits = []

        counter_id = (self.totalPaths()*bounds[0])
        for old_pth in left_paths:
            ### setup new path
            oe_pt_array = 1
            trajIDs = old_pth.trajIDs
            raw_trajs = [ self.panel.raw_trajs[i] for i in trajIDs] #extract raw trajectories
            oe_pt_array = np.array([ np.median(   [rt[i] for rt in raw_trajs]    ) for i in range(bounds[0],bounds[1])])
            oe_path = pathlet2(bounds = bounds ,pt_array = oe_pt_array,id=counter_id,trajIDs=trajIDs,parent=old_pth)
            old_extrapolated_path_candidates.append(oe_path)
            counter_id = counter_id + 1

            ### store cost
            raw_trajs_inbounds = [ self.panel.raw_trajs[i][bounds[0]:bounds[1]] for i in trajIDs] #extract raw trajectories
            total_dist = sum( [self.apply_metric(np.array(oe_path.pt_array),np.array(a)) for a in raw_trajs_inbounds  ] )
            cost = self.c3*total_dist +  self.c5 
            marginal_benefits.append([cost/len(trajIDs),oe_path])

        new_path_candidates = right_paths

        for new_pth in new_path_candidates:
            trajIDs = new_pth.trajIDs
            raw_trajs_inbounds = [ self.panel.raw_trajs[i][bounds[0]:bounds[1]] for i in trajIDs]
            total_dist = sum( [self.apply_metric(np.array(new_pth.pt_array),np.array(a)) for a in raw_trajs_inbounds  ] )
            cost = self.c1 + self.c3*total_dist + (self.c4)*(len(trajIDs)) + self.c5 #cost is cost of new superpath + points that move to it + new path
            marginal_benefits.append([cost/len(trajIDs),new_pth])

        marginal_benefits = sorted(marginal_benefits, key=lambda x: x[0])

        #3) add to pths until 
        # we don't want to add more clusters than the max of k and oe. 
        # reduce cluster size if we add all of k first

        #initial parameters
        chosen_pths = []
        left_lim = len(right_paths)
        left_count = 0
        right_lim = len(left_paths)
        right_count = 0
        total_lim = max(left_lim,right_lim)
        total_count = 0
        trajs_added = []
        overlap_threshold = (1/(total_lim)) # the more sets, the lower (and stricter) the threshold (may want to adapt this)

        for pick in marginal_benefits:
            pth = pick[1] #pick second object (the fisrt is score)

            iter_overlap = overlap(pth.trajIDs,trajs_added)    #comp similarity trajs added so far : prop of trajs already covered 
            if iter_overlap < overlap_threshold: #check if overlap is below the threshold
                chosen_pths.append(pth)
                trajs_added.extend(pth.trajIDs)                
                #add counts
                total_count = total_count + 1
                if pth in right_paths:
                    #print("in right")
                    right_count = right_count + 1
                else:
                    #print("in left")
                    left_count = left_count + 1
                    #if an old path - add child
                    parent_pth = pth.parent
                    if parent_pth.child is not None:
                        print("error - trying to assing child to parent with child")
                    parent_pth.addChild(pth)
                    #print(parent_pth)
                    #print(parent_pth.child)
            else:
                #print("too much overlap, skipping")
                pass

            #break if we reach any of the following limits
            if left_count >= left_lim:
                break
            elif right_count >= right_lim:
                break
            elif total_count >= total_lim:
                break
        
        #4) add new paths to res
        #a) remove existing pths from panel and trajs
        for right_path in right_paths:
            self.pathlets.remove(right_path)
            for traj_id in right_path.trajIDs:
                traj = self.panel.trajs[traj_id] #these match (have confirmed)
                traj.removePathlet(right_path)
        
        #add new paths
        for chosen_pth in chosen_pths:
            #print(chosen_pth)
            self.addPathlet(chosen_pth)

        #b) then for each traj, add to their nearest traj

        for traj in self.panel.trajs:
            current_cost = math.inf
            current_path = 1
            traj_path = traj_path = traj.pt_array[bounds[0]:bounds[1]]
            for pth in chosen_pths:
                pth_path = pth.pt_array
                cost = self.apply_metric(np.array(traj_path),np.array(pth_path))
                if pth.parent in traj.assignments: #cost fall if same as parent (to be consistent with greedy 1)
                    cost = cost - self.c4
                if cost < current_cost:
                    current_cost = cost
                    current_path = pth

            traj.addPathlet(current_path)
            current_path.addTrajIds([traj.trajID])
        
        #3) add each pathlet to assignments
            



    def genSuperPathlets(self):
        ### merge pathlets that are linked to each other (i.e. have the same parent)

        base_pths = [pth for pth in self.pathlets if pth.parent is None]
        
        super_pths = list()

        counter_id = 1
        for base_pth in base_pths:
            #initialise superpath
            super_pth = superPath(id=counter_id,pt_array=base_pth.pt_array,bounds=base_pth.bounds)
            #add on children
            child = base_pth.child
            while child is not None:
                super_pth.extendPath(child.pt_array)
                child = child.child #add new child
            counter_id = counter_id + 1
            super_pths.append(super_pth)
        
        self.superPths = super_pths
        

