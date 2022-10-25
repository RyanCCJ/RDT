# packages
import gc
import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from statannot import add_stat_annotation

# main
def tool1_2(data):
    gene_path = data['gene_path']
    title_map_gene = data['title_map_gene']
    title_map_rna = data['title_map_rna']

    for ANALYZE_TYPE in data['ANALYZE_TYPE']:
        #print('===={}===='.format(ANALYZE_TYPE))
        
        if ANALYZE_TYPE == 'TARGET_SCORE':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['TARGET_SCORE_FILE']:
                    for gene in data['ANALYZE_GROUP']:
                        for FILTER in data['LEN_FILTER']:
                            gene_name = check_new(data, gene)
                            utr5_density, cds_density, utr3_density = main_filter_mRNA_target_score(FILE_NAME, rc, gene_path, gene, gene_name, FILTER)
                            main_plot_target_score(utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna[FILE_NAME], gene, title_map_gene[gene],  rc, FILTER, data['img_title'])
        
        elif ANALYZE_TYPE == 'PIRSCAN':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['PIR_FILE']:
                    for gene in data['ANALYZE_GROUP']:
                        for FILTER in data['LEN_FILTER']:
                            gene_name = check_new(data, gene)
                            utr5_density, cds_density, utr3_density = main_filter_mRNA_pirscan(FILE_NAME, rc, gene_path, gene, gene_name, FILTER)
                            main_plot_pirscan(utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna[FILE_NAME], gene, title_map_gene[gene], rc, FILTER, data['img_title'])         

        elif ANALYZE_TYPE == 'RNAUP':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['RNAUP_FILE']:
                    for gene in data['ANALYZE_GROUP']:
                        for FILTER in data['LEN_FILTER']:
                            gene_name = check_new(data, gene)
                            utr5_density, cds_density, utr3_density = main_filter_mRNA_rnaup(FILE_NAME, rc, gene_path, gene, gene_name, FILTER)
                            main_plot_rnaup(utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna[FILE_NAME], gene, title_map_gene[gene], rc, FILTER, data['img_title'])
                            
        elif ANALYZE_TYPE == 'G22':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['G22_FILE']:
                    for gene in data['ANALYZE_GROUP']:
                        for FILTER in data['LEN_FILTER']:
                            gene_name = check_new(data, gene)

                            #### all region
                            #all_density = main_filter_mRNA_g22_new(FILE_NAME, rc, gene_path, gene, gene_name, FILTER)
                            #main_plot_g22_new(all_density, FILE_NAME, title_map_rna[FILE_NAME], gene, title_map_gene[gene], rc, FILTER)
                            
                            #### 3 regions
                            all_density, utr5_density, cds_density, utr3_density = main_filter_mRNA_g22(FILE_NAME, rc, gene_path, gene, gene_name, FILTER)
                            main_plot_g22(all_density, utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna[FILE_NAME], gene, title_map_gene[gene], rc, FILTER, data['img_title'])

def add_two_mRNA_list(new, gene_path, gene_index, gene):
    if gene_index == 0:
        trans = gene.split('\n') 
        df = new[new['Gene name'].isin(trans)].reset_index(drop=True)
        return df
    else:
        if gene == 'all mRNAs':
            return new
        else:
            data_path = os.path.join(gene_path, gene)
            with open(os.path.abspath(__file__+'/../../'+data_path), 'r') as f:
                trans = f.readlines()
            trans = [n.replace('\n', '') for n in trans]  
            df = new[new['Gene name'].isin(trans)].reset_index(drop=True)
            return df

def check_new(data, gene):
    tran_name = 'new_transcript'
    if gene == 0:
        gene_name = data['NEW_ANALYZE_GROUP']
        data['title_map_gene'][gene] = tran_name
    else:
        gene_name = data['title_map_gene'][gene]
    return gene_name

# test
def KS_test(x,y):
    if x and y:
        less = stats.mstats.ks_2samp(x, y,alternative = 'less')[1]
        greater = stats.mstats.ks_2samp(x, y,alternative = 'greater')[1]
      #  if (x.all() == y.all()) or (sum(x) == 0 and sum(y) == 0): #樣本相同 或 兩個樣本皆為0
      #      two_sided = 1.0
      #  else:
        two_sided = stats.mstats.ks_2samp(x, y,alternative = 'two-sided')[1]
    else:
        less = greater = two_sided = 0
    #print(two_sided, less, greater)
    return [two_sided, round(less,6), round(greater,6)]

def T_test(x,y):
    if x and y:
        d, two_sided = stats.ttest_ind(x, y, equal_var=False)
        if d < 0:
            greater = 1 - two_sided/2 #"greater" is the alternative that x has a larger mean than y
            less = two_sided/2 #"less" is the alternative that x has a larger mean than y
        elif d >0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0
    #print(two_sided, greater, less)
    return [two_sided, round(greater,6), round(less,6)]

def U_test(x,y):
    if x and y:
        d, two_sided = stats.ranksums(x, y)
        if d < 0:
            greater = 1 - two_sided/2
            less = two_sided/2
        elif d >0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0

    #print(two_sided, greater,less )
    return [two_sided, round(greater,6), round(less,6)]

def color_decide(val):
    P_CUTOFF = 0.01
    if val < P_CUTOFF:
        return 'r'
    else:
        return 'black'

def reject_outliers(data, m = 1.5):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

# 1
def main_filter_mRNA_target_score(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    #group = mRNA_list(data, gene)
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene)
    #print('group'+str(gene)+' mRNA=',len(group))
    
    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    utr5_density = done_group['count5'] / done_group['len5']

    done_group = group[(group['lencds']>=FILTER)]
    cds_density = done_group['countcds'] / done_group['lencds']

    done_group = group[(group['len3']>=FILTER)]
    utr3_density = done_group['count3'] / done_group['len3']
    #print(len(utr5_density), len(cds_density), len(utr3_density))
    return utr5_density, cds_density, utr3_density

def main_plot_target_score(utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, img_title):
    # preprocess for plot
    total_utr5 = [0.000001 for i in utr5_density if i == 0] + [i for i in utr5_density if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density if i == 0] + [i for i in cds_density if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density if i == 0] + [i for i in utr3_density if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]

    single_plot_target_score(total_utr5, total_cds, total_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title)

def single_plot_target_score(all_utr5, all_cds, all_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER,utr5_density, cds_density, utr3_density, img_title):
    #res_path = os.path.abspath(__file__+'/../../../target_score_'+need_rc+'_output/density_fig/'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+'/../../../../../static/paper/density_fig/single/density_'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+rc+'.png')
    #print('ploting group:',FILE_NAME,rc,FILTER,str(gene))

    #### non-Log ####
    # U-test
    plt.figure(figsize=(8, 8), dpi=300)
    plt.title(img_title+'\n\n\n\n', fontsize=25)
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    UTR5 = "5'UTR"
    CDS = "CDS"
    UTR3 = "3'UTR"
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000),CDS:pd.Series(cds_density*1000),UTR3:pd.Series(utr3_density*1000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm")
    if rc == 'need_rc':
        ax.set_ylabel('read counts x $\mathregular{10^3}$ / nt',fontsize=22)
    else:
        ax.set_ylabel('sites x $\mathregular{10^3}$ / nt',fontsize=22)
    try:
        add_stat_annotation(ax,data=tmp2,
                        box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    except ValueError:
        print('=========== ValueError =========== ')
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png')) 
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    

# 2
def main_filter_mRNA_pirscan(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    #print('group'+str(gene)+' mRNA=',len(group))
    
    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    utr5_density = done_group['count5'] / done_group['len5']

    done_group = group[(group['lencds']>=FILTER)]
    cds_density = done_group['countcds'] / done_group['lencds']

    done_group = group[(group['len3']>=FILTER)]
    utr3_density = done_group['count3'] / done_group['len3']
    #print(len(utr5_density), len(cds_density), len(utr3_density))
    return utr5_density, cds_density, utr3_density

def main_plot_pirscan(utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, img_title):
    # preprocess for plot
    total_utr5 = [0.000001 for i in utr5_density if i == 0] + [i for i in utr5_density if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density if i == 0] + [i for i in cds_density if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density if i == 0] + [i for i in utr3_density if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]

    single_plot_pirscan(total_utr5, total_cds, total_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title)

def single_plot_pirscan(all_utr5, all_cds, all_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title):
    #res_path = os.path.abspath(__file__+'/../../../pirscan_'+need_rc+'_output/density_fig/'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+'/../../../../../static/paper/density_fig/single/density_'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+rc+'.png')
    #print('ploting group:',FILE_NAME,rc,FILTER,str(gene))

    #### non-Log ####
    # U-test
    plt.figure(figsize=(8, 8), dpi=300)
    plt.title(img_title+'\n\n\n\n', fontsize=25)
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    UTR5 = "5'UTR"
    CDS = "CDS"
    UTR3 = "3'UTR"
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000),CDS:pd.Series(cds_density*1000),UTR3:pd.Series(utr3_density*1000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm")
    if rc == 'need_rc':
        ax.set_ylabel('read counts x $\mathregular{10^3}$ / nt',fontsize=22)
    else:
        ax.set_ylabel('sites x $\mathregular{10^3}$ / nt',fontsize=22)
    try:
        add_stat_annotation(ax,data=tmp2,
                        box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    except ValueError:
        print('=========== ValueError =========== ')
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png')) 
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()


# 3
def main_filter_mRNA_rnaup(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    #print('group'+str(gene)+' mRNA=',len(group))
    
    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    utr5_density = done_group['count5'] / done_group['len5']

    done_group = group[(group['lencds']>=FILTER)]
    cds_density = done_group['countcds'] / done_group['lencds']

    done_group = group[(group['len3']>=FILTER)]
    utr3_density = done_group['count3'] / done_group['len3']
    #print(len(utr5_density), len(cds_density), len(utr3_density))
    return utr5_density, cds_density, utr3_density

def main_plot_rnaup(utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, img_title):
    # preprocess for plot
    #total_utr5 = [i for i in utr5_density if i != 0]
    #total_utr52 = [math.log10(i*1000000) for i in total_utr5]

    #total_cds = [i for i in cds_density if i != 0]
    #total_cds2 = [math.log10(i*1000000) for i in total_cds]

    #total_utr3 = [i for i in utr3_density if i != 0]
    #total_utr32 = [math.log10(i*1000000) for i in total_utr3]

    total_utr5 = [0.000001 for i in utr5_density if i == 0] + [i for i in utr5_density if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density if i == 0] + [i for i in cds_density if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density if i == 0] + [i for i in utr3_density if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]
    
    #single_plot_rnaup(total_utr52, total_cds2, total_utr32, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, total_utr5, total_cds, total_utr3)
    single_plot_rnaup(total_utr5, total_cds, total_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title)

def single_plot_rnaup(all_utr5, all_cds, all_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title):
    #res_path = os.path.abspath(__file__+'/../../../RNAup_'+need_rc+'_output/density_fig/'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+'/../../../../../static/paper/density_fig/single/density_non_zero_'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+rc+'.png')
    #print('ploting group:',FILE_NAME,rc,FILTER,str(gene))

    #### non-Log ####
    # U-test
    plt.figure(figsize=(8, 8), dpi=300)
    plt.title(img_title+'\n\n\n\n', fontsize=25)
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    
    UTR5 = "5'UTR"
    CDS = "CDS"
    UTR3 = "3'UTR"
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000),CDS:pd.Series(cds_density*1000),UTR3:pd.Series(utr3_density*1000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm")
    if rc == 'need_rc':
        ax.set_ylabel('read counts x $\mathregular{10^3}$ / nt',fontsize=22)
    else:
        ax.set_ylabel('sites x $\mathregular{10^3}$ / nt',fontsize=22)
    try:
        add_stat_annotation(ax, data=tmp2,
                        box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    except ValueError:
        print('=========== ValueError =========== ')
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png')) 
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

# 4

def main_filter_mRNA_g22(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    #print('group'+str(gene)+' mRNA=',len(group))

    # all mRNA
    group['all_rc'] = group['count5']+group['countcds']+group['count3']
    group['all_len'] = group['len5']+group['lencds']+group['len3']
    all_density = group['all_rc'] / group['all_len']
    
    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    utr5_density = done_group['count5'] / done_group['len5']

    done_group = group[(group['lencds']>=FILTER)]
    cds_density = done_group['countcds'] / done_group['lencds']

    done_group = group[(group['len3']>=FILTER)]
    utr3_density = done_group['count3'] / done_group['len3']
    #print(len(utr5_density), len(cds_density), len(utr3_density))

    return all_density, utr5_density, cds_density, utr3_density   

def main_filter_mRNA_g22_with_overlap(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+'/../../../../G22_near_site/usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool1_with_overlap.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    #print('group'+str(gene)+' mRNA=',len(group))
    
    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    utr5_density = done_group['count5'] / done_group['len5']

    done_group = group[(group['lencds']>=FILTER)]
    cds_density = done_group['countcds'] / done_group['lencds']

    done_group = group[(group['len3']>=FILTER)]
    utr3_density = done_group['count3'] / done_group['len3']
    #print(len(utr5_density), len(cds_density), len(utr3_density))
    return utr5_density, cds_density, utr3_density

def main_filter_mRNA_g22_new(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    all_group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    #print('group'+str(gene)+' mRNA=',len(all_group))
    
    # FILTER LENGTH
    '''
    done_group = group[(group['len5']>=FILTER)]
    utr5_density = done_group['count5'] / done_group['len5']

    done_group = group[(group['lencds']>=FILTER)]
    cds_density = done_group['countcds'] / done_group['lencds']

    done_group = group[(group['len3']>=FILTER)]
    utr3_density = done_group['count3'] / done_group['len3']
    '''
    all_group['all_rc'] = all_group['count5']+all_group['countcds']+all_group['count3']
    all_group['all_len'] = all_group['len5']+all_group['lencds']+all_group['len3']
    
    all_density = all_group['all_rc'] / all_group['all_len']
    #print(len(all_density))
    return all_density

def main_plot_g22_new(all_density, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, img_title):
    # preprocess for plot
    total_all = [0.000001 for i in all_density if i == 0] + [i for i in all_density if i != 0]
    total_all = [math.log10(i*1000000) for i in total_all]
    
    single_plot_g22_new(total_all, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, all_density, img_title) 

def main_plot_g22_with_overlap(utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, img_title):
    # preprocess for plot
    total_utr5 = [0.000001 for i in utr5_density if i == 0] + [i for i in utr5_density if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density if i == 0] + [i for i in cds_density if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density if i == 0] + [i for i in utr3_density if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]

    single_plot_g22_with_overlap(total_utr5, total_cds, total_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title)

def main_plot_g22(all_density, utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, img_title):

    # preprocess for plot
    total_all = [0.000001 for i in all_density if i == 0] + [i for i in all_density if i != 0]
    total_all = [math.log10(i*1000000) for i in total_all]

    total_utr5 = [0.000001 for i in utr5_density if i == 0] + [i for i in utr5_density if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density if i == 0] + [i for i in cds_density if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density if i == 0] + [i for i in utr3_density if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]

    single_plot_g22(total_all, total_utr5, total_cds, total_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, all_density, utr5_density, cds_density, utr3_density, img_title, img_title)

def single_plot_g22_new(all_region, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, all_density, img_title):
    res_path = os.path.abspath(__file__+'/../../../../../static/paper/density_fig/single/density_ALL_'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    #print('ploting group:',FILE_NAME,rc,FILTER,str(gene))

    #### Log ####
    plt.figure(figsize=(8, 6), dpi=300)
    #all_region_75 = np.percentile(all_region, 75)
    #all_region_25 = np.percentile(all_region, 25)
    #all_region_max = all_region_75 + 1.5*(all_region_75 - all_region_25)
    #all_region_min = all_region_25 - 1.5*(all_region - all_region_25)
    plt.title(img_title+'\n\n\n\n', fontsize=25) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    tmp2 = pd.DataFrame({'ALL':pd.Series(all_region)})
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('log10(read counts*M/nt)',fontsize=22) #log10(read counts*M/nt)
    else:
        ax.set_ylabel('log10(sites counts*M/nt)',fontsize=22)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

    #### non-Log ####
    plt.figure(figsize=(8, 6), dpi=300)
    plt.title(img_title+'\n\n\n\n', fontsize=25) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    tmp2 = pd.DataFrame({'ALL':pd.Series(all_density*1000000)})
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm", showmeans=True)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('read counts*M/nt',fontsize=22)
    else:
        ax.set_ylabel('sites counts*M/nt',fontsize=22)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

def single_plot_g22_with_overlap(all_utr5, all_cds, all_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title):
    #res_path = os.path.abspath(__file__+'/../../../22G_output/density_fig/'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+'/../../../../../static/paper/density_fig/'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'_with_overlap.png')
    #print('ploting group:',FILE_NAME,rc,FILTER,str(gene))

    #### Log ####
    # U-test
    plt.figure(figsize=(14, 8), dpi=300)
    utr5_75 = np.percentile(all_utr5, 75)
    utr3_75 = np.percentile(all_utr3,75 )
    cds_75 = np.percentile(all_cds, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr3_25 = np.percentile(all_utr3,25 )
    cds_25 = np.percentile(all_cds, 25)
    utr5_max = utr5_75 + 1.5*(utr5_75 - utr5_25)
    cds_max = cds_75 + 1.5*(cds_75 - cds_25)
    utr3_max = utr3_75 + 1.5*(utr3_75 - utr3_25)
    utr5_min = utr5_25 - 1.5*(utr5_75 - utr5_25)
    cds_min = cds_25 - 1.5*(cds_75 - cds_25)
    utr3_min = utr3_25 - 1.5*(utr3_75 - utr3_25)
    plt.title(img_title+'\n\n\n\n', fontsize=25) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(all_utr5)),3)
    median_5 = round(np.median(all_utr5),3)
    mean_cds = round(np.mean(np.array(all_cds)),3)
    median_cds = round(np.median(all_cds),3)
    mean_3 = round(np.mean(np.array(all_utr3)),3)
    median_3 = round(np.median(all_utr3),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),CDS:pd.Series(all_cds),UTR3:pd.Series(all_utr3)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('log10(read counts*M/nt)',fontsize=22)
    else:
        ax.set_ylabel('log10(sites counts*M/nt)',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_U.png')) 
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    plt.figure(figsize=(14, 8), dpi=300)
    utr5_75 = np.percentile(all_utr5, 75)
    utr3_75 = np.percentile(all_utr3,75 )
    cds_75 = np.percentile(all_cds, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr3_25 = np.percentile(all_utr3,25 )
    cds_25 = np.percentile(all_cds, 25)
    utr5_max = utr5_75 + 1.5*(utr5_75 - utr5_25)
    cds_max = cds_75 + 1.5*(cds_75 - cds_25)
    utr3_max = utr3_75 + 1.5*(utr3_75 - utr3_25)
    utr5_min = utr5_25 - 1.5*(utr5_75 - utr5_25)
    cds_min = cds_25 - 1.5*(cds_75 - cds_25)
    utr3_min = utr3_25 - 1.5*(utr3_75 - utr3_25)
    plt.title(img_title+'\n\n\n\n', fontsize=20) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(all_utr5)),3)
    median_5 = round(np.median(all_utr5),3)
    mean_cds = round(np.mean(np.array(all_cds)),3)
    median_cds = round(np.median(all_cds),3)
    mean_3 = round(np.mean(np.array(all_utr3)),3)
    median_3 = round(np.median(all_utr3),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),CDS:pd.Series(all_cds),UTR3:pd.Series(all_utr3)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('log10(read counts*M/nt)',fontsize=22)
    else:
        ax.set_ylabel('log10(sites counts*M/nt)',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_T.png'))  
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    #### non-Log ####
    # U-test
    plt.figure(figsize=(14, 8), dpi=300)
    plt.title(img_title+'\n\n\n\n', fontsize=25) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(utr5_density*1000000)),3)
    median_5 = round(np.median(utr5_density*1000000),3)
    mean_cds = round(np.mean(np.array(cds_density*1000000)),3)
    median_cds = round(np.median(cds_density*1000000),3)
    mean_3 = round(np.mean(np.array(utr3_density*1000000)),3)
    median_3 = round(np.median(utr3_density*1000000),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000000),CDS:pd.Series(cds_density*1000000),UTR3:pd.Series(utr3_density*1000000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('read counts*M/nt',fontsize=22)
    else:
        ax.set_ylabel('sites counts*M/nt',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png')) 
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    plt.figure(figsize=(14, 8), dpi=300)
    plt.title(img_title+'\n\n\n\n', fontsize=25) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(utr5_density*1000000)),3)
    median_5 = round(np.median(utr5_density*1000000),3)
    mean_cds = round(np.mean(np.array(cds_density*1000000)),3)
    median_cds = round(np.median(cds_density*1000000),3)
    mean_3 = round(np.mean(np.array(utr3_density*1000000)),3)
    median_3 = round(np.median(utr3_density*1000000),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000000),CDS:pd.Series(cds_density*1000000),UTR3:pd.Series(utr3_density*1000000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('read counts*M/nt',fontsize=22)
    else:
        ax.set_ylabel('sites counts*M/nt',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_T.png'))
    
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

def single_plot_g22(all_region, all_utr5, all_cds, all_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, all_density, utr5_density, cds_density, utr3_density, img_title):
    
    res_path = os.path.abspath(__file__+'/../../../../../static/paper/density_fig/single/density_'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    #print('ploting group:',FILE_NAME,rc,FILTER,str(gene))
    
    #### non-Log ####
    UTR5 = "5'UTR"
    CDS = "CDS"
    UTR3 = "3'UTR"
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000),CDS:pd.Series(cds_density*1000),UTR3:pd.Series(utr3_density*1000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    
    # U-test
    plt.figure(figsize=(8, 8), dpi=300)
    plt.title(img_title+'\n\n\n\n', fontsize=25)
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm")
    if rc == 'need_rc':
        ax.set_ylabel('read counts x $\mathregular{10^3}$ / nt',fontsize=22)
    else:
        ax.set_ylabel('sites x $\mathregular{10^3}$ / nt',fontsize=22)
    try:
        add_stat_annotation(ax,data=tmp2,
                        box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                        test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    except ValueError:
        print('=========== ValueError =========== ')
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png')) 
    plt.clf()
    plt.close()
    #del tmp2
    gc.collect()
    

# per site density
def main_filter_mRNA_g22_per_site(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool1_per_site.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    #print('group'+str(gene)+' mRNA=',len(group))
    
    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    done_group = group[(group['site5']>0)]
    utr5_density = done_group['count5'] / done_group['site5']
    
    done_group = group[(group['lencds']>=FILTER)]
    done_group = group[(group['sitecds']>0)]
    cds_density = done_group['countcds'] / done_group['sitecds']
    
    done_group = group[(group['len3']>=FILTER)]
    done_group = group[(group['site3']>0)]
    utr3_density = done_group['count3'] / done_group['site3']
    #print(len(utr5_density), len(cds_density), len(utr3_density))
    return utr5_density, cds_density, utr3_density

def main_filter_mRNA_g22_per_rc(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool1_per_rc.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    #print('group'+str(gene)+' mRNA=',len(group))
    
    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    done_group = group[(group['site5']>0)]
    utr5_density = done_group['count5'] / done_group['rc5']
    
    done_group = group[(group['lencds']>=FILTER)]
    done_group = group[(group['sitecds']>0)]
    cds_density = done_group['countcds'] / done_group['rccds']
    
    done_group = group[(group['len3']>=FILTER)]
    done_group = group[(group['site3']>0)]
    utr3_density = done_group['count3'] / done_group['rc3']
    #print(len(utr5_density), len(cds_density), len(utr3_density))
    return utr5_density, cds_density, utr3_density  

def main_plot_g22_per_site(utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, img_title):
    # preprocess for plot
    total_utr5 = [0.000001 for i in utr5_density if i == 0] + [i for i in utr5_density if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density if i == 0] + [i for i in cds_density if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density if i == 0] + [i for i in utr3_density if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]

    single_plot_g22_per_site(total_utr5, total_cds, total_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title)  

def main_plot_g22_per_rc(utr5_density, cds_density, utr3_density, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, img_title):
    # preprocess for plot
    total_utr5 = [0.000001 for i in utr5_density if i == 0] + [i for i in utr5_density if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density if i == 0] + [i for i in cds_density if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density if i == 0] + [i for i in utr3_density if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]

    single_plot_g22_per_rc(total_utr5, total_cds, total_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title)  

def single_plot_g22_per_site(all_utr5, all_cds, all_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title):
    #res_path = os.path.abspath(__file__+'/../../../22G_output/density_fig/'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+'/../../../../../static/paper/density_fig/per_site/'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'_per_site.png')
    #print('ploting group:',FILE_NAME,rc,FILTER,str(gene))

    #### Log ####
    # U-test
    plt.figure(figsize=(14, 8), dpi=300)
    utr5_75 = np.percentile(all_utr5, 75)
    utr3_75 = np.percentile(all_utr3,75 )
    cds_75 = np.percentile(all_cds, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr3_25 = np.percentile(all_utr3,25 )
    cds_25 = np.percentile(all_cds, 25)
    utr5_max = utr5_75 + 1.5*(utr5_75 - utr5_25)
    cds_max = cds_75 + 1.5*(cds_75 - cds_25)
    utr3_max = utr3_75 + 1.5*(utr3_75 - utr3_25)
    utr5_min = utr5_25 - 1.5*(utr5_75 - utr5_25)
    cds_min = cds_25 - 1.5*(cds_75 - cds_25)
    utr3_min = utr3_25 - 1.5*(utr3_75 - utr3_25)
    plt.title(title_map_rna+'_'+title_map_gene+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3))+'\n\n\n\n', fontsize=20) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(all_utr5)),3)
    median_5 = round(np.median(all_utr5),3)
    mean_cds = round(np.mean(np.array(all_cds)),3)
    median_cds = round(np.median(all_cds),3)
    mean_3 = round(np.mean(np.array(all_utr3)),3)
    median_3 = round(np.median(all_utr3),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),CDS:pd.Series(all_cds),UTR3:pd.Series(all_utr3)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('log10(read counts*M/sites)',fontsize=22)
    else:
        ax.set_ylabel('log10(sites counts*M/sites)',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_U.png')) 
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    plt.figure(figsize=(14, 8), dpi=300)
    utr5_75 = np.percentile(all_utr5, 75)
    utr3_75 = np.percentile(all_utr3,75 )
    cds_75 = np.percentile(all_cds, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr3_25 = np.percentile(all_utr3,25 )
    cds_25 = np.percentile(all_cds, 25)
    utr5_max = utr5_75 + 1.5*(utr5_75 - utr5_25)
    cds_max = cds_75 + 1.5*(cds_75 - cds_25)
    utr3_max = utr3_75 + 1.5*(utr3_75 - utr3_25)
    utr5_min = utr5_25 - 1.5*(utr5_75 - utr5_25)
    cds_min = cds_25 - 1.5*(cds_75 - cds_25)
    utr3_min = utr3_25 - 1.5*(utr3_75 - utr3_25)
    plt.title(title_map_rna+'_'+title_map_gene+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3))+'\n\n\n\n', fontsize=20) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(all_utr5)),3)
    median_5 = round(np.median(all_utr5),3)
    mean_cds = round(np.mean(np.array(all_cds)),3)
    median_cds = round(np.median(all_cds),3)
    mean_3 = round(np.mean(np.array(all_utr3)),3)
    median_3 = round(np.median(all_utr3),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),CDS:pd.Series(all_cds),UTR3:pd.Series(all_utr3)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('log10(read counts*M/sites)',fontsize=22)
    else:
        ax.set_ylabel('log10(sites counts*M/sites)',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_T.png'))  
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    #### non-Log ####
    # U-test
    plt.figure(figsize=(14, 8), dpi=300)
    plt.title(title_map_rna+'_'+title_map_gene+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3))+'\n\n\n\n', fontsize=20) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(utr5_density*1000000)),3)
    median_5 = round(np.median(utr5_density*1000000),3)
    mean_cds = round(np.mean(np.array(cds_density*1000000)),3)
    median_cds = round(np.median(cds_density*1000000),3)
    mean_3 = round(np.mean(np.array(utr3_density*1000000)),3)
    median_3 = round(np.median(utr3_density*1000000),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000000),CDS:pd.Series(cds_density*1000000),UTR3:pd.Series(utr3_density*1000000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('read counts*M/sites',fontsize=22)
    else:
        ax.set_ylabel('sites counts*M/sites',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png')) 
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    plt.figure(figsize=(14, 8), dpi=300)
    plt.title(title_map_rna+'_'+title_map_gene+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3))+'\n\n\n\n', fontsize=20) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(utr5_density*1000000)),3)
    median_5 = round(np.median(utr5_density*1000000),3)
    mean_cds = round(np.mean(np.array(cds_density*1000000)),3)
    median_cds = round(np.median(cds_density*1000000),3)
    mean_3 = round(np.mean(np.array(utr3_density*1000000)),3)
    median_3 = round(np.median(utr3_density*1000000),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000000),CDS:pd.Series(cds_density*1000000),UTR3:pd.Series(utr3_density*1000000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('read counts*M/sites',fontsize=22)
    else:
        ax.set_ylabel('sites counts*M/sites',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_T.png'))
    
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

def single_plot_g22_per_rc(all_utr5, all_cds, all_utr3, FILE_NAME, title_map_rna, gene, title_map_gene, rc, FILTER, utr5_density, cds_density, utr3_density, img_title):
    #res_path = os.path.abspath(__file__+'/../../../22G_output/density_fig/'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+'/../../../../../static/paper/density_fig/per_site/'+FILE_NAME+'_G'+str(gene)+'_L'+str(FILTER)+'_per_rc.png')
    #print('ploting group:',FILE_NAME,rc,FILTER,str(gene))

    #### Log ####
    # U-test
    plt.figure(figsize=(14, 8), dpi=300)
    utr5_75 = np.percentile(all_utr5, 75)
    utr3_75 = np.percentile(all_utr3,75 )
    cds_75 = np.percentile(all_cds, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr3_25 = np.percentile(all_utr3,25 )
    cds_25 = np.percentile(all_cds, 25)
    utr5_max = utr5_75 + 1.5*(utr5_75 - utr5_25)
    cds_max = cds_75 + 1.5*(cds_75 - cds_25)
    utr3_max = utr3_75 + 1.5*(utr3_75 - utr3_25)
    utr5_min = utr5_25 - 1.5*(utr5_75 - utr5_25)
    cds_min = cds_25 - 1.5*(cds_75 - cds_25)
    utr3_min = utr3_25 - 1.5*(utr3_75 - utr3_25)
    plt.title(title_map_rna+'_'+title_map_gene+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3))+'\n\n\n\n', fontsize=20) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(all_utr5)),3)
    median_5 = round(np.median(all_utr5),3)
    mean_cds = round(np.mean(np.array(all_cds)),3)
    median_cds = round(np.median(all_cds),3)
    mean_3 = round(np.mean(np.array(all_utr3)),3)
    median_3 = round(np.median(all_utr3),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),CDS:pd.Series(all_cds),UTR3:pd.Series(all_utr3)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('log10(read counts*M/sites)',fontsize=22)
    else:
        ax.set_ylabel('log10(sites counts*M/sites)',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_U.png')) 
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    plt.figure(figsize=(14, 8), dpi=300)
    utr5_75 = np.percentile(all_utr5, 75)
    utr3_75 = np.percentile(all_utr3,75 )
    cds_75 = np.percentile(all_cds, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr3_25 = np.percentile(all_utr3,25 )
    cds_25 = np.percentile(all_cds, 25)
    utr5_max = utr5_75 + 1.5*(utr5_75 - utr5_25)
    cds_max = cds_75 + 1.5*(cds_75 - cds_25)
    utr3_max = utr3_75 + 1.5*(utr3_75 - utr3_25)
    utr5_min = utr5_25 - 1.5*(utr5_75 - utr5_25)
    cds_min = cds_25 - 1.5*(cds_75 - cds_25)
    utr3_min = utr3_25 - 1.5*(utr3_75 - utr3_25)
    plt.title(title_map_rna+'_'+title_map_gene+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3))+'\n\n\n\n', fontsize=20) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(all_utr5)),3)
    median_5 = round(np.median(all_utr5),3)
    mean_cds = round(np.mean(np.array(all_cds)),3)
    median_cds = round(np.median(all_cds),3)
    mean_3 = round(np.mean(np.array(all_utr3)),3)
    median_3 = round(np.median(all_utr3),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),CDS:pd.Series(all_cds),UTR3:pd.Series(all_utr3)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('log10(read counts*M/sites)',fontsize=22)
    else:
        ax.set_ylabel('log10(sites counts*M/sites)',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_T.png'))  
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    #### non-Log ####
    # U-test
    plt.figure(figsize=(14, 8), dpi=300)
    plt.title(title_map_rna+'_'+title_map_gene+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3))+'\n\n\n\n', fontsize=20) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(utr5_density*1000000)),3)
    median_5 = round(np.median(utr5_density*1000000),3)
    mean_cds = round(np.mean(np.array(cds_density*1000000)),3)
    median_cds = round(np.median(cds_density*1000000),3)
    mean_3 = round(np.mean(np.array(utr3_density*1000000)),3)
    median_3 = round(np.median(utr3_density*1000000),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000000),CDS:pd.Series(cds_density*1000000),UTR3:pd.Series(utr3_density*1000000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('read counts*M/sites',fontsize=22)
    else:
        ax.set_ylabel('sites counts*M/sites',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png')) 
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    plt.figure(figsize=(14, 8), dpi=300)
    plt.title(title_map_rna+'_'+title_map_gene+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3))+'\n\n\n\n', fontsize=20) #Algorithm_name+'_'+
    plt.xticks(fontsize = 22) 
    plt.yticks(fontsize = 22)
    mean_5 = round(np.mean(np.array(utr5_density*1000000)),3)
    median_5 = round(np.median(utr5_density*1000000),3)
    mean_cds = round(np.mean(np.array(cds_density*1000000)),3)
    median_cds = round(np.median(cds_density*1000000),3)
    mean_3 = round(np.mean(np.array(utr3_density*1000000)),3)
    median_3 = round(np.median(utr3_density*1000000),3)
    UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
    CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
    UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)
    
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density*1000000),CDS:pd.Series(cds_density*1000000),UTR3:pd.Series(utr3_density*1000000)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    ax.set_xlabel('region',fontsize=22)
    if rc == 'need_rc':
        ax.set_ylabel('read counts*M/sites',fontsize=22)
    else:
        ax.set_ylabel('sites counts*M/sites',fontsize=22)
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=0, fontsize=14)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_T.png'))
    
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
