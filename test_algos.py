from clusteringAlgos import *
from models import *

import pandas as pd
import random


df_abnormal_returns = pd.read_csv("data/normalised_price_data.csv")
start_date = 2800
df_preproc = df_abnormal_returns
df = df_preproc[df_preproc["date_id"] >= start_date]
company_names = list(df.columns)[2:]
trajectories = []
for i in company_names:
    trajectories.append(list(df[i]))
traj_panel = Panel(trajObj_list = trajectories,start_time=1)

C1,C2,C3,C4,C5 = 1000,1,1,1,1
K_LIST = [3]
SEG_LIST = [2]

random.seed(1)


res_stc_g1 = greedy1(PANEL = traj_panel,
                     k_list= K_LIST,
                     seg_list = SEG_LIST, 
                     c1=C1,
                     c2=C2,
                     c3=C3,
                     c4=C4,
                     c5=C5)

print("g1")
print(str(res_stc_g1))
#for i in res_stc_g1.pathlets:
#    print(i)


"""

res_stc_n0 = naive_0(PANEL = traj_panel,
              k= 2,
              c1=C1,
              c2=C2,
              c3=C3,
              c4=C4)

print("loss N1",str(res_stc_n0))




res_stc_n1 = naive_1(PANEL = traj_panel,
              k_list= K_LIST,
              c1=C1,
              c2=C2,
              c3=C3,
              c4=C4,
              c5=C5)

print("n1")
print(str(res_stc_n1))



res_stc_n2 = naive_2(PANEL = traj_panel,
              k_list= K_LIST,
              c1=C1,
              c2=C2,
              c3=C3)

print("loss N2" ,str(res_stc_n2.loss()))


res_stc_n3 = naive_3(PANEL = traj_panel,
                     k_list= K_LIST,
                     seg_list = SEG_LIST, 
                     c1=C1,
                     c2=C2,
                     c3=C3,
                     c4=C4,
                     greedy_in_period=False
                     )

print("n3")
print(str(res_stc_n3))

#for i in res_stc_n3.superPths:
#    print(i)

"""
res_stc_g1 = greedy1_1(PANEL = traj_panel,
                     k_list= K_LIST,
                     seg_list = SEG_LIST, 
                     c1=C1,
                     c2=C2,
                     c3=C3,
                     c4=C4,
                     c5=C5)

print("g1")
print(str(res_stc_g1))
#for i in res_stc_g1.pathlets:
#    print(i)



"""
res_stc_g2 = greedy2(PANEL = traj_panel,
                     k_list= K_LIST,
                     seg_list = SEG_LIST, 
                     c1=C1,
                     c2=C2,
                     c3=C3,
                     c4=C4,
                     c5=C5)

print("g2")
print(str(res_stc_g2))
"""

#for i in res_stc_g2.pathlets:
#    print(i)


random.seed(3)

res_stc_g3 = greedy3(PANEL = traj_panel,
                     k_list= K_LIST,
                     seg_list = SEG_LIST, 
                     c1=C1,
                     c2=C2,
                     c3=C3,
                     c4=C4,
                     c5=C5)

print("g3")
print(str(res_stc_g3))



"""

for i in res_stc_g3.superPths:
    print(i)


for i in res_stc_g3.pathlets:
    print(i)
    print(i.parent)
    print(i.child)

"""