import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import kendalltau, spearmanr
import scipy.stats
from sklearn import preprocessing
import matplotlib.pyplot as plt
from pandas import DataFrame
from itertools import chain
import cv2
from PIL import Image
import os
from sklearn.model_selection import train_test_split
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import random
import h5py 
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
from sklearn.metrics import adjusted_rand_score
from tqdm import tqdm

test_mode=0
dataset_name='Human_DLPFC_151673'
rerandom=0

def cal_dist(loc1,loc2,dis_threshold):
    dis_mtx=[]
    neighbours=[]
    for i in range(len(loc1)):
        dis_mtx.append([])
        neighbours.append([])
        for j in range(len(loc1)):
            dis_mtx[i].append(np.sqrt((loc1[i]-loc1[j])**2+(loc2[i]-loc2[j])**2))
            if dis_mtx[i][j] < dis_threshold:
                neighbours[i].append(j)
    return dis_mtx, neighbours

def cal_value_label(cluster_num,cluster_label,expr,neighbours):
    value_label=[0]
    for i in range(1,cluster_num):
        cluster_now=[ind for ind, x in enumerate(cluster_label) if x == i]
        caculated=[0 for i in range(len(cluster_label))]
        high=0
        low=0
        for in_cluster_spot in cluster_now:
            for spot in neighbours[in_cluster_spot]:
                if cluster_label[spot] != i and caculated[spot] == 0:
                    if expr[spot] >= expr[in_cluster_spot]:
                        low+=1
                    else:
                        high+=1
                    caculated[spot]=1
        if high > low:
            value_label.append(1)
        elif high < low:
            value_label.append(-1)
        else:
            value_label.append(0)
        '''
        print(str(i)+':')
        print(low)
        print(high)
        '''
    return value_label

def my_DBScan(loc1,loc2,expr,value_threshold,dis_threshold=1.5,min_spot=4):
    global dis_mtx
    global neighbours
    
    n_spot=len(loc1)
    cluster_label=[0 for i in range(n_spot)]
    candidate=[i for i in range(n_spot)]
    cluster_id=1
    while(len(candidate)):
        now=random.sample(candidate,1)[0]
        cluster_now=[now]
        expand=[now]
        while(len(expand)):
            now=expand[0]
            n_r=0
            candi_join=[]
            candi_expand=[]
            for spot in neighbours[now]:
                if abs(expr[spot]-expr[now]) < value_threshold:
                    if spot not in cluster_now:
                        candi_join.append(spot)
                    if spot in candidate and spot not in expand:
                        candi_expand.append(spot)
                    n_r+=1
            #if(n_r+1 >= min_spot):
            for i in candi_join:
                cluster_now.append(i)
            for i in candi_expand:
                expand.append(i)
            candidate.remove(now)
            expand.pop(0)
        if len(cluster_now) >= min_spot:
            for i in cluster_now:
                cluster_label[i]=cluster_id
            cluster_id+=1
    
    value_label=cal_value_label(cluster_id,cluster_label,expr,neighbours)
    return cluster_label, value_label, cluster_id

def random_select(candidate, num, random_times):
    if num == candidate:
        return -1
    
    global rerandom
    theta = candidate / 10
    res=[]
    for i in range(random_times):
        prob=num/candidate
        l=np.random.binomial(1,prob,candidate)
        while abs(sum(l) - num) > theta:
    	    l=np.random.binomial(1,prob,candidate)
    	    rerandom+=1
        res.append([ind for ind, x in enumerate(l) if x == 1])
    return res
    
def random_cal(use_act,act_name,high_domain_list,low_domain_list,high_random_list,low_random_list,domain_labels,domain_label_set,loc1,loc2,random_times=1000,test_mode=0):
    sig_act=[]
    jaccard_list=[]
    expr_level_list=[]
    explained_domain_list=[]
    p_v_list=[]
    for ge in range(len(act_name)):
        high_domain=high_domain_list[ge]
        low_domain=low_domain_list[ge]
        set_high=set(high_domain)
        set_low=set(low_domain)
        #high_loc=[1 if x in set_high else 0 for x in range(len(loc1))]
        #low_loc=[1 if x in set_low else 0 for x in range(len(loc1))]
        if test_mode and ge==1:
            print(len(high_domain))
            print(len(low_domain))
            print(len(loc1))

        high_random=high_random_list[ge]
        low_random=low_random_list[ge]
            
        #high_loc_random=[]
        #low_loc_random=[]
        #for i in range(random_times):
        #    high_loc_random.append([1 if x in high_random else 0 for x in range(len(loc1))])
        #    low_loc_random.append([1 if x in low_random else 0 for x in range(len(loc1))])

        jaccard_max=-1
        expr_level=0
        p_v_now=1
        for i in range(len(domain_label_set)):
            domain_inds=[ind for ind, x in enumerate(domain_labels) if x == domain_label_set[i]]
            set_domain=set(domain_inds)
            #domain_loc=[1 if x in set_domain else 0 for x in range(len(loc1))]
            jaccard_high=len(set_domain.intersection(set_high))/len(set_domain.union(set_high))
            jaccard_low=len(set_domain.intersection(set_low))/len(set_domain.union(set_low))
            #jaccard_high=adjusted_rand_score(high_loc,domain_loc)
            #jaccard_low=adjusted_rand_score(low_loc,domain_loc)
            h_p=random_times
            l_p=random_times
            for j in range(random_times):
                set_r_h=set(high_random[j])
                set_r_l=set(low_random[j])
                jaccard_r_h=len(set_domain.intersection(set_r_h))/len(set_domain.union(set_r_h))
                jaccard_r_l=len(set_domain.intersection(set_r_l))/len(set_domain.union(set_r_l))
                #jaccard_r_h=adjusted_rand_score(high_loc_random[j],domain_loc)
                #jaccard_r_l=adjusted_rand_score(low_loc_random[j],domain_loc)
                if jaccard_high > jaccard_r_h:
                    h_p-=1
                if jaccard_low > jaccard_r_l:
                    l_p-=1
            h_p=(h_p+1)/(random_times+1)
            l_p=(l_p+1)/(random_times+1)
            if test_mode and ge==1:
                print(domain_label_set[i]+':')
                print('h_p:'+str(h_p))
                print('l_p:'+str(l_p))
                print('jaccard_high:'+str(jaccard_high))
                print('jaccard_low:'+str(jaccard_low)+'\n')
        
            jaccard=-1
            if h_p < 0.01 and h_p < l_p and len(high_domain)/len(loc1) < 0.6:
                jaccard=jaccard_high
                p_v=h_p
                tmp=1
            elif l_p < 0.01 and l_p < h_p and len(low_domain)/len(loc1) < 0.6:
                jaccard=jaccard_low
                p_v=l_p
                tmp=-1
            if jaccard > jaccard_max and p_v < 0.01:
                jaccard_max=jaccard
                expr_level=tmp
                explained_domain=domain_label_set[i]
                p_v_now=p_v
        if p_v_now < 0.01:
            sig_act.append(act_name[ge])
            jaccard_list.append(jaccard_max)
            expr_level_list.append(expr_level)
            explained_domain_list.append(explained_domain)
            p_v_list.append(p_v_now)
    
        pbar.update(1)

    if len(sig_act) >= 1:
        p_v_sorted=[]
        for i in range(len(p_v_list)):
	        p_v_sorted.append([p_v_list[i],i])
        p_v_sorted.sort(key=lambda x : x[0])

        tot=len(p_v_sorted)
        for i in range(tot-2,-1,-1):
            p_v_sorted[i][0]=min(p_v_sorted[i+1][0],p_v_sorted[i][0]*tot/(i+1))
    
        q_v_list=[1 for i in range(len(p_v_sorted))]
        for i in range(len(p_v_sorted)):
            q_v_list[p_v_sorted[i][1]]=p_v_sorted[i][0]

        if test_mode==0:
            df_res=pd.DataFrame(columns = ['','jaccard','expr_level','explained_domain','p_v','q_v'])
            df_res['']=sig_act
            df_res['jaccard']=jaccard_list
            df_res['expr_level']=expr_level_list
            df_res['explained_domain']=explained_domain_list
            df_res['p_v']=p_v_list
            df_res['q_v']=q_v_list
            df_res=df_res.set_index('')
            df_res.to_csv('/home/zkl/Desktop/9.3/'+dataset_name+'/step6/random_test/'+'random_'+str(t)+'_dbscan_'+use_act+'_res_p_v.csv')

def pre_svip():
    res = pd.read_excel('/home/zkl/Desktop/9.3/'+dataset_name+'/step5/svip_flitered.xlsx',index_col=0)
    res = res[res['qval']<0.05]
    svip = res.index.values

    df_clusters = pd.read_csv('/home/zkl/Desktop/9.3/'+dataset_name+'/step5/cluster/cluster.csv',header=0,index_col=0)
    edges_top = df_clusters.index

    svip_clustered=[]
    for i in svip:
        if i in edges_top:
            svip_clustered.append(i)
 
    df_res = df_clusters
    df_250 = df_res.iloc[:,0:1]
    df_250 = df_250[df_250.index.isin(svip_clustered)]
    print(df_250)

    a = df_250.values.tolist()
    b = [] # 250条代表边关于250个簇的标签
    for i in range(len(a)):
       b.append(a[i][0])
    labels_ = np.array(b)


    ### 代表元素在各spot上的活性
    df_act = pd.read_csv(r'/home/zkl/Desktop/9.3/'+dataset_name+'/step3/act_new.csv', sep='\t', index_col=0)
    df_act = df_act[df_act.index.isin(svip_clustered)]
    df_act = df_act.T
    print(df_act)
    act_name=df_act.columns
    #loc=pd.read_csv('/home/zkl/Desktop/9.3/data/'+dataset_name+'/loc.csv', sep=',', index_col=0)
    list_act = df_act.T.values.tolist()
    
    #loc1 = loc['x'].values
    #loc2 = loc['y'].values
    
    pbar=tqdm(total=len(act_name))
    high_domain_list=[]
    low_domain_list=[]
    high_random_list=[]
    low_random_list=[]
    #dis_mtx, neighbours=cal_dist(loc1,loc2,dis_threshold=1.5)
    global dis_mtx, neighbours

    for ge in range(len(act_name)):
        value_threshold=(max(list_act[ge])-min(list_act[ge]))/10
        cluster_label, value_label, cluster_num = my_DBScan(loc1, loc2, list_act[ge], value_threshold)
        
        inds=[]
        for i in range(cluster_num):
            inds.append([ind for ind, x in enumerate(cluster_label) if x == i])
        
        high_domain=[]
        low_domain=[]
        for i in range(1,cluster_num):
            if value_label[i]==1:
                high_domain.extend(inds[i])
            elif value_label[i]==-1:
                low_domain.extend(inds[i])
        high_domain_list.append(high_domain)
        low_domain_list.append(low_domain)
        
        random_times=1000
        high_random=random_select(len(loc1),len(high_domain),random_times)
        low_random=random_select(len(loc1),len(low_domain),random_times)
        high_random_list.append(high_random)
        low_random_list.append(low_random)
        
        pbar.update(1)
    
    return list_act, act_name, high_domain_list, low_domain_list,  high_random_list, low_random_list

### SVG
def pre_svg():
    res = pd.read_excel('/home/zkl/Desktop/9.3/'+dataset_name+'/step1/SpatialDE_res.xlsx')
    res = res[res['qval']<0.05]
    svg = res['g'].values

    df_count = pd.read_csv('/home/zkl/Desktop/9.3/data/'+dataset_name+'/count.csv', sep='\t', index_col=0)
    df_act = df_count[[x for x in svg]]
    print(df_act)
    act_name=df_act.columns
    #loc=pd.read_csv('/home/zkl/Desktop/9.3/data/'+dataset_name+'/loc.csv', sep=',', index_col=0)
    list_act = df_act.T.values.tolist()
    
    #loc1 = loc['x'].values
    #loc2 = loc['y'].values
    
    pbar=tqdm(total=len(act_name))
    high_domain_list=[]
    low_domain_list=[]
    high_random_list=[]
    low_random_list=[]
    #dis_mtx, neighbours=cal_dist(loc1,loc2,dis_threshold=1.5)
    global dis_mtx, neighbours

    for ge in range(len(act_name)):
        value_threshold=(max(list_act[ge])-min(list_act[ge]))/10
        cluster_label, value_label, cluster_num = my_DBScan(loc1, loc2, list_act[ge], value_threshold)
        
        inds=[]
        for i in range(cluster_num):
            inds.append([ind for ind, x in enumerate(cluster_label) if x == i])
        
        high_domain=[]
        low_domain=[]
        for i in range(1,cluster_num):
            if value_label[i]==1:
                high_domain.extend(inds[i])
            elif value_label[i]==-1:
                low_domain.extend(inds[i])
        high_domain_list.append(high_domain)
        low_domain_list.append(low_domain)
        
        random_times=1000
        high_random=random_select(len(loc1),len(high_domain),random_times)
        low_random=random_select(len(loc1),len(low_domain),random_times)
        high_random_list.append(high_random)
        low_random_list.append(low_random)
        
        pbar.update(1)
    
    return list_act, act_name, high_domain_list, low_domain_list, high_random_list, low_random_list

    

### 空间域识别结果(index要对应)
domains= pd.read_excel('/home/zkl/Desktop/9.3/data/'+dataset_name+'/domain.xlsx', header=0, index_col=0)
print(domains)
domain_labels = np.array(domains.iloc[:,2].tolist())
domain_label_set = list(set(domain_labels))


###location
loc=pd.read_csv('/home/zkl/Desktop/9.3/data/'+dataset_name+'/loc.csv', sep='\t', index_col=0)
#loc = loc.loc[domains.index,:]
loc1 = loc['x'].values
loc2 = loc['y'].values

dis_mtx, neighbours=cal_dist(loc1,loc2,dis_threshold=1.5)

colors=['silver','green','blue','yellow','orange','r','cyan','purple','lightpink','sandybrown','olive']

svip_list_act, svip_act_name, svip_high_domain_list, svip_low_domain_list, svip_high_random_list, svip_low_random_list = pre_svip()
svg_list_act, svg_act_name, svg_high_domain_list, svg_low_domain_list, svg_high_random_list, svg_low_random_list = pre_svg()
pbar=tqdm(total=(len(svip_act_name)+len(svg_act_name))*100)

for t in range(100):
    #random.shuffle(domain_labels)
    domains= pd.read_excel('/home/zkl/Desktop/9.3/data/'+dataset_name+'/random/GraphST_'+str(t)+'.xlsx', header=0, index_col=0)
    #print(domains)
    domain_labels = np.array(domains.iloc[:,0].tolist())
    domain_label_set = list(set(domain_labels))
    '''
    plt.figure(dpi=300,figsize=(6,4))
    for i in range(len(domain_label_set)):
        inds = [ind for ind, x in enumerate(domain_labels) if x == domain_label_set[i]]
        plt_x=[loc1[k] for k in inds]
        plt_y=[loc2[k] for k in inds]
        plt.scatter(plt_x,plt_y,color=colors[i%len(colors)],s=20)
    plt.savefig('/home/zkl/Desktop/9.3/'+dataset_name+'/step6/random_test/random_'+str(t)+'_domain.png')
    plt.cla()
    plt.clf()
    plt.close()
    '''
    random_cal('svip',svip_act_name,svip_high_domain_list,svip_low_domain_list,svip_high_random_list,svip_low_random_list,domain_labels,domain_label_set,loc1,loc2)
    random_cal('svg',svg_act_name,svg_high_domain_list,svg_low_domain_list,svg_high_random_list,svg_low_random_list,domain_labels,domain_label_set,loc1,loc2)
    
    print(t)
    

print(rerandom)

    
    

