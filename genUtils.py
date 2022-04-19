import numpy as np

def l1_dist(x,y):
    """
    Computes 1 distance between vectors x and y
    """
    return(np.linalg.norm(x - y, ord=1))



    
def comp_diff(s1,s2):
    """ compare the similarity between s1 and s2 (can approx with hashing if nessesary) 
    #0.5 if exactly the same, 1 if completely different

    Args:
        s1,s2: two sets
    """
    total_element = len(list(set(s1 + s2)))
    total_size = len(s1 + s2)
    return(total_element/total_size)
