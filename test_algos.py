from clusteringAlgos import *
from models import *

import pandas as pd


df_abnormal_returns = pd.read_csv("data/abnormal_returns_data.csv")
start_date = 3000
df_preproc = df_abnormal_returns
df = df_preproc[df_preproc["date_id"] >= start_date]
company_names = list(df.columns)[2:]
trajectories = []
for i in company_names:
    trajectories.append(list(df[i]))
traj_panel = Panel(trajObj_list = trajectories,start_time=1)

C1,C2,C3,C4 = 500,1,1,1
K_LIST = [5]
SEG_LIST = [4]

"""
res_stc_n0 = naive_0(PANEL = traj_panel,
              k= 2,
              c1=C1,
              c2=C2,
              c3=C3,
              c4=C4)

print("loss N1",str(res_stc_n0))
"""


"""
res_stc_n1 = naive_1(PANEL = traj_panel,
              k_list= K_LIST,
              c1=C1,
              c2=C2,
              c3=C3,
              c4=C4)


print("loss N1",str(res_stc_n1))
"""

"""
res_stc_n2 = naive_2(PANEL = traj_panel,
              k_list= K_LIST,
              c1=C1,
              c2=C2,
              c3=C3)

print("loss N2" ,str(res_stc_n2.loss()))
"""

"""
res_stc_n3 = naive_3(PANEL = traj_panel,
                     k_list= K_LIST,
                     seg_list = SEG_LIST, 
                     c1=C1,
                     c2=C2,
                     c3=C3,
                     c4=C4,
                     greedy_in_period=False
                     )

print(str(res_stc_n3))

for i in res_stc_n3.superPths:
    print(i)
"""

res_stc_g1 = greedy1(PANEL = traj_panel,
                     k_list= K_LIST,
                     seg_list = SEG_LIST, 
                     c1=C1,
                     c2=C2,
                     c3=C3,
                     c4=C4)


print(str(res_stc_g1))


