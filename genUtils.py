import numpy as np

def l1_dist(x,y):
    """
    Computes 1 distance between vectors x and y
    """
    return(np.linalg.norm(x - y, ord=1))



    
def comp_diff(s1,s2,type="Jaccard"):
    """ compare the similarity between s1 and s2 (can approx with hashing if nessesary) 
    # 1 if exactly the same, 0 if completely different

    Args:
        s1,s2: two list objects
    """
    A = set(s1)
    B = set(s2)

    if type == "Jaccard":
        num = len(A.intersection(B))# a intersect B
        denom  = len(set(s1 + s2)) # a union b
        return(num/denom)
    else:
        print("not using jaccard")
        #this is a weird metric, lets not use it!
        total_element = len(list(set(s1 + s2)))
        total_size = len(s1 + s2)
        return( 2*(1 - total_element/total_size))


def overlap(a,b):
    """ exports what proportion of a is in b
    Args:
        a,b: two lists

    Return:
        proportion of a that is in b
    """
    a_clean = list(set(a))
    b_clean = list(set(b))

    a_in_b = [i for i in a_clean if i in b_clean]

    return(len(a_in_b)/len(a))


