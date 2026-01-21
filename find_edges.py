import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from pandas.core.frame import DataFrame
import scanpy as sc
import scipy as sp
import matplotlib.image as mpimg
from matplotlib import rcParams
import seaborn as sb
import anndata as ad
import h5py
from PIL import Image
import cv2
import sys
from scipy.stats.mstats import gmean
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

'''
count_csv: spot_id * gene_id
'''
def find_edges(svg_excel, count_csv, species, edges_qv_out, act_edge_name_out_excel, act_gene_count_out_csv, gene_network='./Network/newHnet-2015name.txt'):
    if species not in ['human','mouse']:
        print('species must be one of below:'+str(['human','mouse']))
        
    df_Human_Interactome = pd.read_csv(gene_network,sep='\t',header=None,index_col=False,names=['a','b','c','d'])
    
    edges_Human_Interactome=[]
    c=df_Human_Interactome['c'].tolist()
    d=df_Human_Interactome['d'].tolist()
    for i in range(len(c)):
        edges_Human_Interactome.append([c[i],d[i]])

    res=pd.read_excel(svg_excel,index_col=0)
    res_svg=res[res['qval']<0.05]
    genes_filtered = res.index
    
    if species == 'mouse':
        f = open('./mouse_human_homologs.txt','r')
        m2h = []
        for i in f:
            m2h.append(i[:-1].split('\t'))
        f.close()

        genes_human = []
        genes_qvalue = []
        for i in range(len(genes_filtered)):
            for j in range(len(m2h)):
                if m2h[j][0] == genes_filtered[i]:
                    if m2h[j][1] not in genes_human:
                        genes_human.append(m2h[j][1])
                        if type(res.loc[genes_filtered[i],'qval']) != np.float64:
                            genes_qvalue.append(res.loc[genes_filtered[i],'qval'].iloc[0])
                        else:
                            genes_qvalue.append(res.loc[genes_filtered[i],'qval'])
                    break
    elif species == 'human':
        genes_human = []
        genes_qvalue = []
        for i in range(len(genes_filtered)):
            if genes_filtered[i] not in genes_human:
                genes_human.append(genes_filtered[i])
                if type(res.loc[genes_filtered[i],'qval']) == pd.Series:
                    genes_qvalue.append(res.loc[genes_filtered[i],'qval'].iloc[0])
                else:
                    genes_qvalue.append(res.loc[genes_filtered[i],'qval'])
    
    
    edges_total = []
    for i in range(len(edges_Human_Interactome)):
        if((c[i] in genes_human) & (d[i] in genes_human)):
            edges_total.append(edges_Human_Interactome[i])
    
    '''
    if ['PTGDS','TBXAS1'] in edges_total:
        print('have!')
    '''
            
    self_loop = []
    for i in range(len(edges_total)):
        if(edges_total[i][0] == edges_total[i][1]):
            self_loop.append(i)
    self_loop.reverse()
    for i in self_loop:
        edges_total.pop(i)
    
    print('total edges:'+str(len(edges_total)))
    print('top 10 edges:'+str(edges_total[0:10]))
    
    
    gene1 = DataFrame(edges_total)[0].tolist()
    gene2 = DataFrame(edges_total)[1].tolist()
    
    res_dic = {}
    for i in range(len(genes_qvalue)):
        if genes_human[i] not in res_dic.keys():
            res_dic[genes_human[i]]=genes_qvalue[i]
    
    gene1_qval = []
    for i in range(len(gene1)):
        gene1_qval.append(float(res_dic[gene1[i]]))
    gene2_qval = []
    for i in range(len(gene2)):
        gene2_qval.append(float(res_dic[gene2[i]]))
    
    min_qval = []
    for i in range(len(gene1_qval)):
        min_qval.append(min(gene1_qval[i],gene2_qval[i]))
    avg_qval = []
    for i in range(len(gene1_qval)):
        avg_qval.append((gene1_qval[i]+gene2_qval[i])/2)
    
    edges_qv = pd.DataFrame(columns = ['edges','gene1_qval','gene2_qval','avg_qval','min_qval'])
    edges_qv['edges']=edges_total
    edges_qv['gene1_qval']=gene1_qval
    edges_qv['gene2_qval']=gene2_qval
    edges_qv['avg_qval']=avg_qval
    edges_qv['min_qval']=min_qval
    edges_qv = edges_qv.set_index('edges')
    
    edges_qv_sotred = edges_qv.sort_values(by=['min_qval','avg_qval'],ascending=[True,True])
    edges_qv_sotred.to_excel(edges_qv_out)
    
    edges_qv_sotred = pd.read_excel(edges_qv_out,index_col=0)
    df_edges = edges_qv_sotred[edges_qv_sotred['min_qval'] <= 0.05]
    gene1 = []
    gene2 = []
    edges_list = df_edges.index.tolist()
    for i in range(len(edges_list)):
        gene1.append(edges_list[i].strip('[').strip(']').split(',')[0].replace("'","").replace("'",""))
        gene2.append(edges_list[i].strip('[').strip(']').split(',')[1].replace("'","").replace("'","").replace(" ",""))
    edges_svg = pd.DataFrame(columns = [0,1])
    edges_svg[0]=gene1
    edges_svg[1]=gene2
    edges_svg.to_excel(act_edge_name_out_excel)
    print('total svg edges:'+str(len(edges_svg)))
    
    
    gene_tmp = gene1+gene2
    new_genes=[]
    for i in gene_tmp:
        if i not in new_genes:
            new_genes.append(i)
    print('total genes involved in svg edges:'+str(len(new_genes)))
    
    df_count = pd.read_csv(count_csv, sep='\t',index_col=0)
    df_count = df_count.T[df_count.sum(0) >= 3].T
    
    if species == 'mouse':
        df_value = df_count.values.T.tolist()
        df_tmp=[]
        re_col=[]
        for i in range(len(df_count.columns)):
            for j in range(len(m2h)):
                if m2h[j][0] == df_count.columns[i]:
                    re_col.append(m2h[j][1])
                    df_tmp.append(df_value[i])
                    break
        df_human = DataFrame(df_tmp,index=re_col,columns=df_count.index).T
    elif species == 'human':
        df_human = df_count
    
    act_genes = DataFrame(df_human[new_genes])
    act_genes.to_csv(act_gene_count_out_csv)
    

def find_CCI_edges(count_csv, species, LR_path, edges_qv_out, act_edge_name_out_excel, act_gene_count_out_csv, gene_network='./Network/newHnet-2015name.txt', min_cell=0, mean='algebra', test_mode=0):
    
    df_count = pd.read_csv(count_csv, sep='\t',index_col=0)
    df_count = df_count.T[df_count.sum(0) >= 3].T
    
    if species == 'mouse':
        df_value = df_count.values.T.tolist()
        df_tmp=[]
        re_col=[]
        for i in range(len(df_count.columns)):
            for j in range(len(m2h)):
                if m2h[j][0] == df_count.columns[i]:
                    re_col.append(m2h[j][1])
                    df_tmp.append(df_value[i])
                    break
        df_human = DataFrame(df_tmp,index=re_col,columns=df_count.index).T
    elif species == 'human':
        df_human = df_count
    
    geneInter = pd.read_csv(LR_path + 'interaction_input_CellChatDB.csv.gz', index_col=0, compression='gzip')
    comp = pd.read_csv(LR_path + 'complex_input_CellChatDB.csv', header=0, index_col=0)
    
    geneInter = geneInter.sort_values('annotation')
    ligand = geneInter.ligand.values
    receptor = geneInter.receptor.values
    geneInter.pop('ligand')
    geneInter.pop('receptor')
    print('total ligands:'+str(len(ligand)))
    print(ligand)
    print('total receptors:'+str(len(receptor)))
    print(receptor)

    t = []
    ligand_selected = []
    receptor_selected = []
    for i in range(len(ligand)):
        for n,flag in [[ligand,0], [receptor,1]]:
            l = n[i]
            if l in comp.index:
                n[i] = comp.loc[l].dropna().values[pd.Series \
                    (comp.loc[l].dropna().values).isin(df_human.columns)]
            else:
                n[i] = pd.Series(l).values[pd.Series(l).isin(df_human.columns)]
        if (len(ligand[i]) > 0) * (len(receptor[i]) > 0):
            if mean=='geometric':
                meanL = gmean(df_human.loc[:,ligand[i]].values, axis=1)
                meanR = gmean(df_human.loc[:, receptor[i]].values, axis=1)
            else:
                meanL = df_human.loc[:,ligand[i]].values.mean(axis=1)
                meanR = df_human.loc[:, receptor[i]].values.mean(axis=1)
            if (sum(meanL > 0) >= min_cell) * \
                    (sum(meanR > 0) >= min_cell):
                ligand_selected.append(ligand[i])
                receptor_selected.append(receptor[i])
    print('selected ligands:'+str(len(ligand_selected)))
    print(ligand_selected[0:10])
    print('selected receptors:'+str(len(receptor_selected)))
    print(receptor_selected[0:10])
    
    act_edges = pd.DataFrame(columns = [0,1])
    act_genes=[]
    cnt=0
    for i in range(len(ligand_selected)):
        for l in ligand_selected[i]:
            for r in receptor_selected[i]:
                if l != r:
                    if str([l,r]) in act_edges.index:
                        cnt+=1
                    act_edges.loc[str([l,r])] = [l,r]
                    if l not in act_genes:
                        act_genes.append(l)
                    if r not in act_genes:
                        act_genes.append(r)
    act_edges.index = range(len(act_edges))
    print('selected edges:'+str(len(act_edges)))
    print(act_edges)
    print('selected genes:'+str(len(act_genes)))
    act_genes = DataFrame(df_human[act_genes])
    print(act_genes.iloc[0:5,0:5])
    
    if test_mode == 0:
        act_edges.to_excel(act_edge_name_out_excel)
        act_genes.to_csv(act_gene_count_out_csv)

if __name__ == '__main__':
    species='human'
    
    dataset_name = sys.argv[1]
    df_count='/home/zkl/Desktop/9.3/data/'+dataset_name+'/count.csv'
    svg_excel='/home/zkl/Desktop/9.3/'+dataset_name+'/step1/SpatialDE_res.xlsx'
    edges_qv_out=111
    act_edge_name_out_excel=111
    act_gene_count_out_csv=111
    LR_path = '/home/zkl/Desktop/SpatialDM-main/spatialdm/datasets/LR_data/human-'
    find_CCI_edges(df_count, species, LR_path, edges_qv_out, act_edge_name_out_excel, act_gene_count_out_csv, test_mode=1)



    
    
    
    
    
    
    
    
    
