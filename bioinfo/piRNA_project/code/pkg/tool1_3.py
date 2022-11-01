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
def tool1_3(data):
    gene_path = data['gene_path']
    title_map_gene = data['title_map_gene']
    title_map_rna = data['title_map_rna']

    for ANALYZE_TYPE in data['ANALYZE_TYPE']:
        #print('===={}===='.format(ANALYZE_TYPE))
        
        if ANALYZE_TYPE == 'TARGET_SCORE':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['TARGET_SCORE_FILE']:
                    for FILTER in data['LEN_FILTER']:#[100,150]:
                        for gene1, gene2 in zip(data['GENE_LIST1'], data['GENE_LIST2']):
                            gene_name1, gene_name2 = check_new(data, gene1, gene2)
                            utr5_density1, cds_density1, utr3_density1 = main_filter_mRNA_target_score(FILE_NAME, rc, gene_path, gene1, gene_name1, FILTER)
                            utr5_density2, cds_density2, utr3_density2 = main_filter_mRNA_target_score(FILE_NAME, rc, gene_path, gene2, gene_name2, FILTER)
                            main_plot_target_score(utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, FILE_NAME, title_map_rna[FILE_NAME], gene1, title_map_gene[gene1], gene2, title_map_gene[gene2], rc, FILTER, data['img_title'])

        elif ANALYZE_TYPE == 'PIRSCAN':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['PIR_FILE']:
                    for FILTER in data['LEN_FILTER']:#[25,50,75,100,125,150]:
                        for gene1, gene2 in zip(data['GENE_LIST1'], data['GENE_LIST2']):
                            gene_name1, gene_name2 = check_new(data, gene1, gene2)
                            utr5_density1, cds_density1, utr3_density1 = main_filter_mRNA_pirscan(FILE_NAME, rc, gene_path, gene1, gene_name1, FILTER)
                            utr5_density2, cds_density2, utr3_density2 = main_filter_mRNA_pirscan(FILE_NAME, rc, gene_path, gene2, gene_name2, FILTER)
                            main_plot_pirscan(utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, FILE_NAME, title_map_rna[FILE_NAME], gene1, title_map_gene[gene1], gene2, title_map_gene[gene2], rc, FILTER)
                            
        elif ANALYZE_TYPE == 'RNAUP':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['RNAUP_FILE']:
                    for FILTER in data['LEN_FILTER']:#[25,50,75,100,125,150]:
                        for gene1, gene2 in zip(data['GENE_LIST1'], data['GENE_LIST2']):
                            gene_name1, gene_name2 = check_new(data, gene1, gene2)
                            utr5_density1, cds_density1, utr3_density1 = main_filter_mRNA_rnaup(FILE_NAME, rc, gene_path, gene1, gene_name1, FILTER)
                            utr5_density2, cds_density2, utr3_density2 = main_filter_mRNA_rnaup(FILE_NAME, rc, gene_path, gene2, gene_name2, FILTER)
                            main_plot_rnaup(utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, FILE_NAME, title_map_rna[FILE_NAME], gene1, title_map_gene[gene1], gene2, title_map_gene[gene2], rc, FILTER)
                            
        elif ANALYZE_TYPE == 'G22':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['G22_FILE']:
                    for FILTER in data['LEN_FILTER']:#[100,150]:
                        for gene1, gene2 in zip(data['GENE_LIST1'], data['GENE_LIST2']):
                            gene_name1, gene_name2 = check_new(data, gene1, gene2)

                            #### all region
                            all_density1 = main_filter_mRNA_g22_new(FILE_NAME, rc, gene_path, gene1, gene_name1, FILTER)
                            all_density2 = main_filter_mRNA_g22_new(FILE_NAME, rc, gene_path, gene2, gene_name2, FILTER)
                            main_plot_g22_new(all_density1, all_density2, FILE_NAME, title_map_rna[FILE_NAME], gene1, title_map_gene[gene1], gene2, title_map_gene[gene2], rc, FILTER)

                            #### 3 regions
                            all_density1, utr5_density1, cds_density1, utr3_density1 = main_filter_mRNA_g22(FILE_NAME, rc, gene_path, gene1, gene_name1, FILTER)
                            all_density2, utr5_density2, cds_density2, utr3_density2 = main_filter_mRNA_g22(FILE_NAME, rc, gene_path, gene2, gene_name2, FILTER)
                            main_plot_g22(all_density1, all_density2, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, FILE_NAME, title_map_rna[FILE_NAME], gene1, title_map_gene[gene1], gene2, title_map_gene[gene2], rc, FILTER)
 

def check_new(data, gene1, gene2):
    tran_name = 'new_transcript'
    if gene1 == 0:
        gene_name1 = data['NEW_GENE_LIST1']
        data['title_map_gene'][gene1] = tran_name
    else:
        gene_name1 = data['title_map_gene'][gene1]
    if gene2 == 0:
        gene_name2 = data['NEW_GENE_LIST2']
        data['title_map_gene'][gene2] = tran_name
    else:
        gene_name2 = data['title_map_gene'][gene2]
    return gene_name1, gene_name2

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

def main_filter_mRNA_target_score(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    #group = mRNA_list(data, gene)
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    print('group'+str(gene)+' mRNA=',len(group))

    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    utr5_density = done_group['count5'] / done_group['len5']

    done_group = group[(group['lencds']>=FILTER)]
    cds_density = done_group['countcds'] / done_group['lencds']

    done_group = group[(group['len3']>=FILTER)]
    utr3_density = done_group['count3'] / done_group['len3']
    print(len(utr5_density), len(cds_density), len(utr3_density))
    return utr5_density, cds_density, utr3_density

def main_plot_target_score(utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc, FILTER, img_title):
    # preprocess for plot
    total_utr5 = [0.000001 for i in utr5_density1 if i == 0] + [i for i in utr5_density1 if i != 0]
    total_utr5 = [math.log10(i*1000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density1 if i == 0] + [i for i in cds_density1 if i != 0]
    total_cds = [math.log10(i*1000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density1 if i == 0] + [i for i in utr3_density1 if i != 0]
    total_utr3 = [math.log10(i*1000) for i in total_utr3]

    total_utr52 = [0.000001 for i in utr5_density2 if i == 0] + [i for i in utr5_density2 if i != 0]
    total_utr52 = [math.log10(i*1000) for i in total_utr52]

    total_cds2 = [0.000001 for i in cds_density2 if i == 0] + [i for i in cds_density2 if i != 0]
    total_cds2 = [math.log10(i*1000) for i in total_cds2]

    total_utr32 = [0.000001 for i in utr3_density2 if i == 0] + [i for i in utr3_density2 if i != 0]
    total_utr32 = [math.log10(i*1000) for i in total_utr32]
    
    single_plot_target_score(total_utr5, total_cds, total_utr3, total_utr52, total_cds2, total_utr32, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2,
                             rc, FILTER, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, img_title)
    
def single_plot_target_score(all_utr5, all_cds, all_utr3, all_utr52, all_cds2, all_utr32, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2,
                             rc, FILTER, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, img_title):
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/density_fig/COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/density_fig/compare_list/density_COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(FILTER)+rc+'.png')
    #print('ploting group:',FILE_NAME,rc,FILTER,str(gene1),str(gene2))
    '''
    #### Log ####
    # U-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(all_utr5),3)
    median_5 = round(np.median(all_utr5),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(all_utr52),3)
    median_52 = round(np.median(all_utr52),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(all_cds),3)
    median_cds = round(np.median(all_cds),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(all_cds2),3)
    median_cds2 = round(np.median(all_cds2),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(all_utr3),3)
    median_3 = round(np.median(all_utr3),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(all_utr32),3)
    median_32 = round(np.median(all_utr32),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(all_utr5, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr5_max = utr5_25 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(all_utr52, 75)
    utr5_25 = np.percentile(all_utr52, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(all_cds, 75)
    cds_25 = np.percentile(all_cds, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(all_cds2, 75)
    cds_25 = np.percentile(all_cds2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(all_utr3, 75)
    utr3_25 = np.percentile(all_utr3, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(all_utr32, 75)
    utr3_25 = np.percentile(all_utr32, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    
    plt.ylim(0, max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),UTR52:pd.Series(all_utr52)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('log10(read counts*M/nt)',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('log10(sites counts*M/nt)',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(all_cds),CDS2:pd.Series(all_cds2)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(all_utr3),UTR32:pd.Series(all_utr32)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    #tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        #title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        #title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_U.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

    # T-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(all_utr5),3)
    median_5 = round(np.median(all_utr5),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(all_utr52),3)
    median_52 = round(np.median(all_utr52),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(all_cds),3)
    median_cds = round(np.median(all_cds),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(all_cds2),3)
    median_cds2 = round(np.median(all_cds2),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(all_utr3),3)
    median_3 = round(np.median(all_utr3),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(all_utr32),3)
    median_32 = round(np.median(all_utr32),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(all_utr5, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr5_max = utr5_25 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(all_utr52, 75)
    utr5_25 = np.percentile(all_utr52, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(all_cds, 75)
    cds_25 = np.percentile(all_cds, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(all_cds2, 75)
    cds_25 = np.percentile(all_cds2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(all_utr3, 75)
    utr3_25 = np.percentile(all_utr3, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(all_utr32, 75)
    utr3_25 = np.percentile(all_utr32, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    
    plt.ylim(0, max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),UTR52:pd.Series(all_utr52)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('log10(read counts*M/nt)',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('log10(sites counts*M/nt)',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(all_cds),CDS2:pd.Series(all_cds2)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(all_utr3),UTR32:pd.Series(all_utr32)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_T.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    '''
    #### non-Log ####
    # U-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.tick_params(axis='y', labelsize=16)
    
    UTR5 = title_map_gene1
    UTR52 = title_map_gene2
    CDS = title_map_gene1
    CDS2 = title_map_gene2
    UTR3 = title_map_gene1
    UTR32 = title_map_gene2

    utr5_75 = np.percentile(utr5_density1, 75)
    utr5_25 = np.percentile(utr5_density1, 25)
    utr5_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(utr5_density2, 75)
    utr5_25 = np.percentile(utr5_density2, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(cds_density1, 75)
    cds_25 = np.percentile(cds_density1, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(cds_density2, 75)
    cds_25 = np.percentile(cds_density2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(utr3_density1, 75)
    utr3_25 = np.percentile(utr3_density1, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(utr3_density2, 75)
    utr3_25 = np.percentile(utr3_density2, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    
    plt.ylim(0, 1000*max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density1*1000),UTR52:pd.Series(utr5_density2*1000)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=0)
    if rc == 'need_rc':
        ax1.set_ylabel('reads x $\mathregular{10^3}$ / nt',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('sites x $\mathregular{10^3}$ / nt',fontsize=22)

    
    tmp2 = pd.DataFrame({CDS:pd.Series(cds_density1*1000),CDS2:pd.Series(cds_density2*1000)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=0)
    
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(utr3_density1*1000),UTR32:pd.Series(utr3_density2*1000)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=0)
    
    ax2.set_title(img_title+'\n\n', fontsize=25)
    #tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
    #                    title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
    #                    title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)

    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png'))
    plt.clf()
    plt.close()
    #del tmp2
    gc.collect()
    '''
    # T-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(utr5_density1*1000000),3)
    median_5 = round(np.median(utr5_density1*1000000),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(utr5_density2*1000000),3)
    median_52 = round(np.median(utr5_density2*1000000),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(cds_density1*1000000),3)
    median_cds = round(np.median(cds_density1*1000000),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(cds_density2*1000000),3)
    median_cds2 = round(np.median(cds_density2*1000000),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(utr3_density1*1000000),3)
    median_3 = round(np.median(utr3_density1*1000000),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(utr3_density2*1000000),3)
    median_32 = round(np.median(utr3_density2*1000000),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(utr5_density1, 75)
    utr5_25 = np.percentile(utr5_density1, 25)
    utr5_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(utr5_density2, 75)
    utr5_25 = np.percentile(utr5_density2, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(cds_density1, 75)
    cds_25 = np.percentile(cds_density1, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(cds_density2, 75)
    cds_25 = np.percentile(cds_density2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(utr3_density1, 75)
    utr3_25 = np.percentile(utr3_density1, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(utr3_density2, 75)
    utr3_25 = np.percentile(utr3_density2, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    plt.ylim(0, 1000000*max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density1*1000000),UTR52:pd.Series(utr5_density2*1000000)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('read counts*M/nt',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('sites counts*M/nt',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(cds_density1*1000000),CDS2:pd.Series(cds_density2*1000000)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(utr3_density1*1000000),UTR32:pd.Series(utr3_density2*1000000)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_T.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    '''

def main_filter_mRNA_pirscan(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+ '/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    print('group'+str(gene)+' mRNA=',len(group))
    
    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    utr5_density = done_group['count5'] / done_group['len5']

    done_group = group[(group['lencds']>=FILTER)]
    cds_density = done_group['countcds'] / done_group['lencds']

    done_group = group[(group['len3']>=FILTER)]
    utr3_density = done_group['count3'] / done_group['len3']
    print(len(utr5_density), len(cds_density), len(utr3_density))
    return utr5_density, cds_density, utr3_density

def main_plot_pirscan(utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc, FILTER):
    # preprocess for plot
    total_utr5 = [0.000001 for i in utr5_density1 if i == 0] + [i for i in utr5_density1 if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density1 if i == 0] + [i for i in cds_density1 if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density1 if i == 0] + [i for i in utr3_density1 if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]

    total_utr52 = [0.000001 for i in utr5_density2 if i == 0] + [i for i in utr5_density2 if i != 0]
    total_utr52 = [math.log10(i*1000000) for i in total_utr52]

    total_cds2 = [0.000001 for i in cds_density2 if i == 0] + [i for i in cds_density2 if i != 0]
    total_cds2 = [math.log10(i*1000000) for i in total_cds2]

    total_utr32 = [0.000001 for i in utr3_density2 if i == 0] + [i for i in utr3_density2 if i != 0]
    total_utr32 = [math.log10(i*1000000) for i in total_utr32]
    
    single_plot_pirscan(total_utr5, total_cds, total_utr3, total_utr52, total_cds2, total_utr32,FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2,
                        rc, FILTER, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2)

def single_plot_pirscan(all_utr5, all_cds, all_utr3, all_utr52, all_cds2, all_utr32,FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc,
                        FILTER, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2):
    #res_path = os.path.abspath(__file__+ '/../../pirscan_'+rc+'_output/density_fig/COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/density_fig/compare_list/density_COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(FILTER)+rc+'.png')
    
    print('ploting group:',FILE_NAME,rc,FILTER,str(gene1),str(gene2))
    
    #### Log ####
    # U-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(all_utr5),3)
    median_5 = round(np.median(all_utr5),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(all_utr52),3)
    median_52 = round(np.median(all_utr52),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(all_cds),3)
    median_cds = round(np.median(all_cds),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(all_cds2),3)
    median_cds2 = round(np.median(all_cds2),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(all_utr3),3)
    median_3 = round(np.median(all_utr3),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(all_utr32),3)
    median_32 = round(np.median(all_utr32),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(all_utr5, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr5_max = utr5_25 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(all_utr52, 75)
    utr5_25 = np.percentile(all_utr52, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(all_cds, 75)
    cds_25 = np.percentile(all_cds, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(all_cds2, 75)
    cds_25 = np.percentile(all_cds2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(all_utr3, 75)
    utr3_25 = np.percentile(all_utr3, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(all_utr32, 75)
    utr3_25 = np.percentile(all_utr32, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    
    plt.ylim(0, max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),UTR52:pd.Series(all_utr52)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('log10(read counts*M/nt)',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('log10(sites counts*M/nt)',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(all_cds),CDS2:pd.Series(all_cds2)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(all_utr3),UTR32:pd.Series(all_utr32)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_U.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(all_utr5),3)
    median_5 = round(np.median(all_utr5),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(all_utr52),3)
    median_52 = round(np.median(all_utr52),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(all_cds),3)
    median_cds = round(np.median(all_cds),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(all_cds2),3)
    median_cds2 = round(np.median(all_cds2),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(all_utr3),3)
    median_3 = round(np.median(all_utr3),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(all_utr32),3)
    median_32 = round(np.median(all_utr32),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(all_utr5, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr5_max = utr5_25 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(all_utr52, 75)
    utr5_25 = np.percentile(all_utr52, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(all_cds, 75)
    cds_25 = np.percentile(all_cds, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(all_cds2, 75)
    cds_25 = np.percentile(all_cds2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(all_utr3, 75)
    utr3_25 = np.percentile(all_utr3, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(all_utr32, 75)
    utr3_25 = np.percentile(all_utr32, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    
    plt.ylim(0, max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),UTR52:pd.Series(all_utr52)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('log10(read counts*M/nt)',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('log10(sites counts*M/nt)',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(all_cds),CDS2:pd.Series(all_cds2)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(all_utr3),UTR32:pd.Series(all_utr32)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_T.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    #### non-Log ####
    # U-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(utr5_density1*1000000),3)
    median_5 = round(np.median(utr5_density1*1000000),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(utr5_density2*1000000),3)
    median_52 = round(np.median(utr5_density2*1000000),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(cds_density1*1000000),3)
    median_cds = round(np.median(cds_density1*1000000),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(cds_density2*1000000),3)
    median_cds2 = round(np.median(cds_density2*1000000),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(utr3_density1*1000000),3)
    median_3 = round(np.median(utr3_density1*1000000),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(utr3_density2*1000000),3)
    median_32 = round(np.median(utr3_density2*1000000),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(utr5_density1, 75)
    utr5_25 = np.percentile(utr5_density1, 25)
    utr5_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(utr5_density2, 75)
    utr5_25 = np.percentile(utr5_density2, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(cds_density1, 75)
    cds_25 = np.percentile(cds_density1, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(cds_density2, 75)
    cds_25 = np.percentile(cds_density2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(utr3_density1, 75)
    utr3_25 = np.percentile(utr3_density1, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(utr3_density2, 75)
    utr3_25 = np.percentile(utr3_density2, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    plt.ylim(0, 1000000*max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density1*1000000),UTR52:pd.Series(utr5_density2*1000000)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('read counts*M/nt',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('sites counts*M/nt',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(cds_density1*1000000),CDS2:pd.Series(cds_density2*1000000)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(utr3_density1*1000000),UTR32:pd.Series(utr3_density2*1000000)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(utr5_density1*1000000),3)
    median_5 = round(np.median(utr5_density1*1000000),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(utr5_density2*1000000),3)
    median_52 = round(np.median(utr5_density2*1000000),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(cds_density1*1000000),3)
    median_cds = round(np.median(cds_density1*1000000),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(cds_density2*1000000),3)
    median_cds2 = round(np.median(cds_density2*1000000),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(utr3_density1*1000000),3)
    median_3 = round(np.median(utr3_density1*1000000),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(utr3_density2*1000000),3)
    median_32 = round(np.median(utr3_density2*1000000),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(utr5_density1, 75)
    utr5_25 = np.percentile(utr5_density1, 25)
    utr5_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(utr5_density2, 75)
    utr5_25 = np.percentile(utr5_density2, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(cds_density1, 75)
    cds_25 = np.percentile(cds_density1, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(cds_density2, 75)
    cds_25 = np.percentile(cds_density2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(utr3_density1, 75)
    utr3_25 = np.percentile(utr3_density1, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(utr3_density2, 75)
    utr3_25 = np.percentile(utr3_density2, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    plt.ylim(0, 1000000*max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density1*1000000),UTR52:pd.Series(utr5_density2*1000000)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('read counts*M/nt',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('sites counts*M/nt',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(cds_density1*1000000),CDS2:pd.Series(cds_density2*1000000)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(utr3_density1*1000000),UTR32:pd.Series(utr3_density2*1000000)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_T.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
 
def main_filter_mRNA_rnaup(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+ '/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    print('group'+str(gene)+' mRNA=',len(group))
    
    # FILTER LENGTH
    done_group = group[(group['len5']>=FILTER)]
    utr5_density = done_group['count5'] / done_group['len5']

    done_group = group[(group['lencds']>=FILTER)]
    cds_density = done_group['countcds'] / done_group['lencds']

    done_group = group[(group['len3']>=FILTER)]
    utr3_density = done_group['count3'] / done_group['len3']
    print(len(utr5_density), len(cds_density), len(utr3_density))
    return utr5_density, cds_density, utr3_density

def single_plot_rnaup(all_utr5, all_cds, all_utr3, all_utr52, all_cds2, all_utr32,FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc,
                      FILTER, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2):
    #res_path = os.path.abspath(__file__+ '/../../../RNAup_'+rc+'_output/density_fig/COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/density_fig/compare_list/density_COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(FILTER)+rc+'.png')
    
    print('ploting group:',FILE_NAME,rc,FILTER,str(gene1),str(gene2))
    
    #### Log ####
    # U-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(all_utr5),3)
    median_5 = round(np.median(all_utr5),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(all_utr52),3)
    median_52 = round(np.median(all_utr52),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(all_cds),3)
    median_cds = round(np.median(all_cds),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(all_cds2),3)
    median_cds2 = round(np.median(all_cds2),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(all_utr3),3)
    median_3 = round(np.median(all_utr3),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(all_utr32),3)
    median_32 = round(np.median(all_utr32),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(all_utr5, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr5_max = utr5_25 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(all_utr52, 75)
    utr5_25 = np.percentile(all_utr52, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(all_cds, 75)
    cds_25 = np.percentile(all_cds, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(all_cds2, 75)
    cds_25 = np.percentile(all_cds2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(all_utr3, 75)
    utr3_25 = np.percentile(all_utr3, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(all_utr32, 75)
    utr3_25 = np.percentile(all_utr32, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    
    plt.ylim(0, max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),UTR52:pd.Series(all_utr52)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('log10(read counts*M/nt)',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('log10(sites counts*M/nt)',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(all_cds),CDS2:pd.Series(all_cds2)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(all_utr3),UTR32:pd.Series(all_utr32)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_U.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(all_utr5),3)
    median_5 = round(np.median(all_utr5),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(all_utr52),3)
    median_52 = round(np.median(all_utr52),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(all_cds),3)
    median_cds = round(np.median(all_cds),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(all_cds2),3)
    median_cds2 = round(np.median(all_cds2),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(all_utr3),3)
    median_3 = round(np.median(all_utr3),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(all_utr32),3)
    median_32 = round(np.median(all_utr32),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(all_utr5, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr5_max = utr5_25 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(all_utr52, 75)
    utr5_25 = np.percentile(all_utr52, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(all_cds, 75)
    cds_25 = np.percentile(all_cds, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(all_cds2, 75)
    cds_25 = np.percentile(all_cds2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(all_utr3, 75)
    utr3_25 = np.percentile(all_utr3, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(all_utr32, 75)
    utr3_25 = np.percentile(all_utr32, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    
    plt.ylim(0, max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),UTR52:pd.Series(all_utr52)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('log10(read counts*M/nt)',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('log10(sites counts*M/nt)',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(all_cds),CDS2:pd.Series(all_cds2)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(all_utr3),UTR32:pd.Series(all_utr32)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_T.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    #### non-Log ####
    # U-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(utr5_density1*1000000),3)
    median_5 = round(np.median(utr5_density1*1000000),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(utr5_density2*1000000),3)
    median_52 = round(np.median(utr5_density2*1000000),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(cds_density1*1000000),3)
    median_cds = round(np.median(cds_density1*1000000),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(cds_density2*1000000),3)
    median_cds2 = round(np.median(cds_density2*1000000),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(utr3_density1*1000000),3)
    median_3 = round(np.median(utr3_density1*1000000),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(utr3_density2*1000000),3)
    median_32 = round(np.median(utr3_density2*1000000),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(utr5_density1, 75)
    utr5_25 = np.percentile(utr5_density1, 25)
    utr5_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(utr5_density2, 75)
    utr5_25 = np.percentile(utr5_density2, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(cds_density1, 75)
    cds_25 = np.percentile(cds_density1, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(cds_density2, 75)
    cds_25 = np.percentile(cds_density2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(utr3_density1, 75)
    utr3_25 = np.percentile(utr3_density1, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(utr3_density2, 75)
    utr3_25 = np.percentile(utr3_density2, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    plt.ylim(0, 1000000*max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density1*1000000),UTR52:pd.Series(utr5_density2*1000000)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('read counts*M/nt',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('sites counts*M/nt',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(cds_density1*1000000),CDS2:pd.Series(cds_density2*1000000)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(utr3_density1*1000000),UTR32:pd.Series(utr3_density2*1000000)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()
    
    # T-test
    fig ,(ax1,ax2,ax3) = plt.subplots(1,3, sharey=True, figsize=(14,8),dpi=200)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax1.tick_params(axis='y', labelsize=16)
    
    mean_5 = round(np.mean(utr5_density1*1000000),3)
    median_5 = round(np.median(utr5_density1*1000000),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(utr5_density2*1000000),3)
    median_52 = round(np.median(utr5_density2*1000000),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(cds_density1*1000000),3)
    median_cds = round(np.median(cds_density1*1000000),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(cds_density2*1000000),3)
    median_cds2 = round(np.median(cds_density2*1000000),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(utr3_density1*1000000),3)
    median_3 = round(np.median(utr3_density1*1000000),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(utr3_density2*1000000),3)
    median_32 = round(np.median(utr3_density2*1000000),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    utr5_75 = np.percentile(utr5_density1, 75)
    utr5_25 = np.percentile(utr5_density1, 25)
    utr5_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(utr5_density2, 75)
    utr5_25 = np.percentile(utr5_density2, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(cds_density1, 75)
    cds_25 = np.percentile(cds_density1, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(cds_density2, 75)
    cds_25 = np.percentile(cds_density2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(utr3_density1, 75)
    utr3_25 = np.percentile(utr3_density1, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(utr3_density2, 75)
    utr3_25 = np.percentile(utr3_density2, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    plt.ylim(0, 1000000*max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density1*1000000),UTR52:pd.Series(utr5_density2*1000000)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax1.set_ylabel('read counts*M/nt',fontsize=22) #log10(read counts*M/nt)
    else:
        ax1.set_ylabel('sites counts*M/nt',fontsize=22)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(cds_density1*1000000),CDS2:pd.Series(cds_density2*1000000)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(utr3_density1*1000000),UTR32:pd.Series(utr3_density2*1000000)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)   
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_T.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

def main_plot_rnaup(utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc, FILTER):
    # preprocess for plot
    total_utr5 = [0.000001 for i in utr5_density1 if i == 0] + [i for i in utr5_density1 if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density1 if i == 0] + [i for i in cds_density1 if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density1 if i == 0] + [i for i in utr3_density1 if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]

    total_utr52 = [0.000001 for i in utr5_density2 if i == 0] + [i for i in utr5_density2 if i != 0]
    total_utr52 = [math.log10(i*1000000) for i in total_utr52]

    total_cds2 = [0.000001 for i in cds_density2 if i == 0] + [i for i in cds_density2 if i != 0]
    total_cds2 = [math.log10(i*1000000) for i in total_cds2]

    total_utr32 = [0.000001 for i in utr3_density2 if i == 0] + [i for i in utr3_density2 if i != 0]
    total_utr32 = [math.log10(i*1000000) for i in total_utr32]
    
    single_plot_rnaup(total_utr5, total_cds, total_utr3, total_utr52, total_cds2, total_utr32,FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2,
                      rc, FILTER, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2)

def main_filter_mRNA_g22(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+ '/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    print('group'+str(gene)+' mRNA=',len(group))
    
    # all_mRNA
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

def main_filter_mRNA_g22_new(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+ '/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    all_group = add_two_mRNA_list(data, gene_path, gene, title_map_gene,)
    print('group'+str(gene)+' mRNA=',len(all_group))
    
    # FILTER LENGTH
    all_group['all_rc'] = all_group['count5']+all_group['countcds']+all_group['count3']
    all_group['all_len'] = all_group['len5']+all_group['lencds']+all_group['len3']
    
    all_density = all_group['all_rc'] / all_group['all_len']
    print(len(all_density))
    
    return all_density

def single_plot_g22_new(all_region1, all_region2, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc, FILTER, all_density1, all_density2):
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/density_fig/compare_list/density_ALL_COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(FILTER)+'.png')
    print('ploting group:',FILE_NAME,rc,FILTER,str(gene1),str(gene2))
    
    #### Log ####
    # U-test
    all_region_75 = np.percentile(all_region1, 75)
    all_region_25 = np.percentile(all_region1, 25)
    all_region1_max = all_region_25 + 1.6*(all_region_75 - all_region_25)
    all_region_75 = np.percentile(all_region2, 75)
    all_region_25 = np.percentile(all_region2, 25)
    all_region2_max = all_region_75 + 1.6*(all_region_75 - all_region_25)

    plt.figure(figsize=(8, 6),dpi=200)
    plt.title(title_map_rna+'\n'+title_map_gene1+' N='+str(len(all_region1))+'\n'+title_map_gene2+' N='+str(len(all_region2)), y=1.1)
    plt.yticks(fontsize = 20)
    plt.ylim(0, max([all_region1_max, all_region2_max]))
    tmp2 = pd.DataFrame({title_map_gene1:pd.Series(all_region1),title_map_gene2:pd.Series(all_region2)})
    my_pal = {title_map_gene1:'#FFA07A', title_map_gene2:'#98F898'}
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3,showmeans=True, palette=my_pal)
    #ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(title_map_gene1, title_map_gene2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        plt.ylabel('log10(read counts*M/nt)',fontsize=20)
    else:
        plt.ylabel('log10(sites counts*M/nt)',fontsize=20)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_U.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

    # T-test
    plt.figure(figsize=(8, 6),dpi=200)
    plt.title(title_map_rna+'\n'+title_map_gene1+' N='+str(len(all_region1))+'\n'+title_map_gene2+' N='+str(len(all_region2)), y=1.1)
    plt.yticks(fontsize = 20)
    plt.ylim(0, max([all_region1_max, all_region2_max]))
    tmp2 = pd.DataFrame({title_map_gene1:pd.Series(all_region1),title_map_gene2:pd.Series(all_region2)})
    my_pal = {title_map_gene1:'#FFA07A', title_map_gene2:'#98F898'}
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3,showmeans=True, palette=my_pal)
    #ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(title_map_gene1, title_map_gene2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        plt.ylabel('log10(read counts*M/nt)',fontsize=20)
    else:
        plt.ylabel('log10(sites counts*M/nt)',fontsize=20)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_T.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

    #### non-Log ####
    # U-test
    all_density_75 = np.percentile(all_density1, 75)
    all_density_25 = np.percentile(all_density1, 25)
    all_density1_max = all_density_25 + 1.6*(all_density_75 - all_density_25)
    all_density_75 = np.percentile(all_density2, 75)
    all_density_25 = np.percentile(all_density2, 25)
    all_density2_max = all_density_75 + 1.6*(all_density_75 - all_density_25)

    plt.figure(figsize=(8, 6),dpi=200)
    plt.title(title_map_rna+'\n'+title_map_gene1+' N='+str(len(all_density1))+'\n'+title_map_gene2+' N='+str(len(all_density2)), y=1.1)
    plt.yticks(fontsize = 20)
    plt.ylim(0, 1000000*max([all_density1_max, all_density2_max]))
    tmp2 = pd.DataFrame({title_map_gene1:pd.Series(all_density1*1000000),title_map_gene2:pd.Series(all_density2*1000000)})
    my_pal = {title_map_gene1:'#FFA07A', title_map_gene2:'#98F898'}
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3,showmeans=True, palette=my_pal)
    #ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(title_map_gene1, title_map_gene2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        plt.ylabel('read counts*M/nt',fontsize=20)
    else:
        plt.ylabel('sites counts*M/nt',fontsize=20)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

    # T-test
    plt.figure(figsize=(8, 6),dpi=200)
    plt.title(title_map_rna+'\n'+title_map_gene1+' N='+str(len(all_density1))+'\n'+title_map_gene2+' N='+str(len(all_density2)), y=1.1)
    plt.yticks(fontsize = 20)
    plt.ylim(0, 1000000*max([all_density1_max, all_density2_max]))
    tmp2 = pd.DataFrame({title_map_gene1:pd.Series(all_density1*1000000),title_map_gene2:pd.Series(all_density2*1000000)})
    my_pal = {title_map_gene1:'#FFA07A', title_map_gene2:'#98F898'}
    ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3,showmeans=True, palette=my_pal)
    #ax = sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette="coolwarm",showmeans=True)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax,data=tmp2,
                    box_pairs=[(title_map_gene1, title_map_gene2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        plt.ylabel('read counts*M/nt',fontsize=20)
    else:
        plt.ylabel('sites counts*M/nt',fontsize=20)
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_T.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

def single_plot_g22(all_region1, all_region2, all_utr5, all_cds, all_utr3, all_utr52, all_cds2, all_utr32, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc, FILTER,
                    all_density1, all_density2, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2):
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/density_fig/compare_list/density_COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(FILTER)+'.png')
    print('ploting group:',FILE_NAME,rc,FILTER,str(gene1),str(gene2))
    
    #### Log ####
    mean_all = round(np.mean(all_region1),3)
    median_all = round(np.median(all_region1),3)
    ALL = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_all, median_all)
    mean_all2 = round(np.mean(all_region2),3)
    median_all2 = round(np.median(all_region2),3)
    ALL2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_all2, median_all2)

    mean_5 = round(np.mean(all_utr5),3)
    median_5 = round(np.median(all_utr5),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(all_utr52),3)
    median_52 = round(np.median(all_utr52),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(all_cds),3)
    median_cds = round(np.median(all_cds),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(all_cds2),3)
    median_cds2 = round(np.median(all_cds2),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(all_utr3),3)
    median_3 = round(np.median(all_utr3),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(all_utr32),3)
    median_32 = round(np.median(all_utr32),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    
    all_75 = np.percentile(all_region1, 75)
    all_25 = np.percentile(all_region1, 25)
    all_max = all_25 + 1.6*(all_75 - all_25)
    all_75 = np.percentile(all_region2, 75)
    all_25 = np.percentile(all_region2, 25)
    all2_max = all_75 + 1.6*(all_75 - all_25)

    utr5_75 = np.percentile(all_utr5, 75)
    utr5_25 = np.percentile(all_utr5, 25)
    utr5_max = utr5_25 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(all_utr52, 75)
    utr5_25 = np.percentile(all_utr52, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(all_cds, 75)
    cds_25 = np.percentile(all_cds, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(all_cds2, 75)
    cds_25 = np.percentile(all_cds2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(all_utr3, 75)
    utr3_25 = np.percentile(all_utr3, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(all_utr32, 75)
    utr3_25 = np.percentile(all_utr32, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)

    # U-test
    fig ,(ax0, ax1,ax2,ax3) = plt.subplots(1,4, sharey=True, figsize=(14,8),dpi=200)

    ax0.set_title(title_map_rna+'_ALL\n'+title_map_gene1+' N='+str(len(all_region1))+'\n'+
                 title_map_gene2+' N='+str(len(all_region2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax0.tick_params(axis='y', labelsize=16)
    plt.ylim(0, max([all_max, all2_max, utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({ALL:pd.Series(all_region1),ALL2:pd.Series(all_region2)})
    my_pal = {ALL:'#FFA07A', ALL2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax0)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax0,data=tmp2,
                    box_pairs=[(ALL, ALL2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax0.set_ylabel('log10(read counts*M/nt)',fontsize=20) #log10(read counts*M/nt)
    else:
        ax0.set_ylabel('log10(sites counts*M/nt)',fontsize=20)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),UTR52:pd.Series(all_utr52)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(all_cds),CDS2:pd.Series(all_cds2)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(all_utr3),UTR32:pd.Series(all_utr32)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_U.png'))
    plt.clf()
    plt.close()
    #del tmp2
    gc.collect()
    
    # T-test
    fig ,(ax0, ax1,ax2,ax3) = plt.subplots(1,4, sharey=True, figsize=(14,8),dpi=200)

    ax0.set_title(title_map_rna+'_ALL\n'+title_map_gene1+' N='+str(len(all_region1))+'\n'+
                 title_map_gene2+' N='+str(len(all_region2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax0.tick_params(axis='y', labelsize=16)
    plt.ylim(0, max([all_max, all2_max, utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({ALL:pd.Series(all_region1),ALL2:pd.Series(all_region2)})
    my_pal = {ALL:'#FFA07A', ALL2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax0)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax0,data=tmp2,
                    box_pairs=[(ALL, ALL2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax0.set_ylabel('log10(read counts*M/nt)',fontsize=20) #log10(read counts*M/nt)
    else:
        ax0.set_ylabel('log10(sites counts*M/nt)',fontsize=20)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    tmp2 = pd.DataFrame({UTR5:pd.Series(all_utr5),UTR52:pd.Series(all_utr52)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(all_cds),CDS2:pd.Series(all_cds2)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(all_utr3),UTR32:pd.Series(all_utr32)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)  
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    
    
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_log_T.png'))
    plt.clf()
    plt.close()
    #del tmp2
    gc.collect()
    
    #### non-Log ####
    mean_all = round(np.mean(all_density1*1000000),3)
    median_all = round(np.median(all_density1*1000000),3)
    ALL = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_all, median_all)
    mean_all2 = round(np.mean(all_density2*1000000),3)
    median_all2 = round(np.median(all_density2*1000000),3)
    ALL2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_all2, median_all2)

    mean_5 = round(np.mean(utr5_density1*1000000),3)
    median_5 = round(np.median(utr5_density1*1000000),3)
    UTR5 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_5, median_5)
    mean_52 = round(np.mean(utr5_density2*1000000),3)
    median_52 = round(np.median(utr5_density2*1000000),3)
    UTR52 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_52, median_52)
    
    mean_cds = round(np.mean(cds_density1*1000000),3)
    median_cds = round(np.median(cds_density1*1000000),3)
    CDS = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_cds, median_cds)
    mean_cds2 = round(np.mean(cds_density2*1000000),3)
    median_cds2 = round(np.median(cds_density2*1000000),3)
    CDS2 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_cds2, median_cds2)
    
    mean_3 = round(np.mean(utr3_density1*1000000),3)
    median_3 = round(np.median(utr3_density1*1000000),3)
    UTR3 = '{}\navg:{}\nmedian:{}'.format(title_map_gene1, mean_3, median_3)
    mean_32 = round(np.mean(utr3_density2*1000000),3)
    median_32 = round(np.median(utr3_density2*1000000),3)
    UTR32 = '{}\navg:{}\nmedian:{}'.format(title_map_gene2, mean_32, median_32)
    
    all_75 = np.percentile(all_density1, 75)
    all_25 = np.percentile(all_density1, 25)
    all_max = all_25 + 1.6*(all_75 - all_25)
    all_75 = np.percentile(all_density2, 75)
    all_25 = np.percentile(all_density2, 25)
    all2_max = all_75 + 1.6*(all_75 - all_25)

    utr5_75 = np.percentile(utr5_density1, 75)
    utr5_25 = np.percentile(utr5_density1, 25)
    utr5_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    utr5_75 = np.percentile(utr5_density2, 75)
    utr5_25 = np.percentile(utr5_density2, 25)
    utr52_max = utr5_75 + 1.6*(utr5_75 - utr5_25)
    
    cds_75 = np.percentile(cds_density1, 75)
    cds_25 = np.percentile(cds_density1, 25)
    cds_max = cds_75 + 1.6*(cds_75 - cds_25)
    cds_75 = np.percentile(cds_density2, 75)
    cds_25 = np.percentile(cds_density2, 25)
    cds2_max = cds_75 + 1.6*(cds_75 - cds_25)
    
    utr3_75 = np.percentile(utr3_density1, 75)
    utr3_25 = np.percentile(utr3_density1, 25)
    utr3_max = utr3_75 + 1.6*(utr3_75 - utr3_25)
    utr3_75 = np.percentile(utr3_density2, 75)
    utr3_25 = np.percentile(utr3_density2, 25)
    utr32_max = utr3_75 + 1.6*(utr3_75 - utr3_25)

    # U-test
    fig ,(ax0, ax1,ax2,ax3) = plt.subplots(1,4, sharey=True, figsize=(14,8),dpi=200)

    ax0.set_title(title_map_rna+'_ALL\n'+title_map_gene1+' N='+str(len(all_density1))+'\n'+
                 title_map_gene2+' N='+str(len(all_density2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax0.tick_params(axis='y', labelsize=16)
    plt.ylim(0, 1000000*max([all_max, all2_max, utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({ALL:pd.Series(all_density1*1000000),ALL2:pd.Series(all_density2*1000000)})
    my_pal = {ALL:'#FFA07A', ALL2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax0)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax0,data=tmp2,
                    box_pairs=[(ALL, ALL2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax0.set_ylabel('read counts*M/nt',fontsize=20) #log10(read counts*M/nt)
    else:
        ax0.set_ylabel('sites counts*M/nt',fontsize=20)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density1*1000000),UTR52:pd.Series(utr5_density2*1000000)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)

    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(cds_density1*1000000),CDS2:pd.Series(cds_density2*1000000)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(utr3_density1*1000000),UTR32:pd.Series(utr3_density2*1000000)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3) 
    
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_U.png'))
    plt.clf()
    plt.close()
    #del tmp2
    gc.collect()
    
    # T-test
    fig ,(ax0,ax1,ax2,ax3) = plt.subplots(1,4, sharey=True, figsize=(14,8),dpi=200)

    ax0.set_title(title_map_rna+'_ALL\n'+title_map_gene1+' N='+str(len(all_density1))+'\n'+
                 title_map_gene2+' N='+str(len(all_density2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    ax0.tick_params(axis='y', labelsize=16)
    plt.ylim(0, 1000000*max([all_max, all2_max, utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({ALL:pd.Series(all_density1*1000000),ALL2:pd.Series(all_density2*1000000)})
    my_pal = {ALL:'#FFA07A', ALL2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax0)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax0,data=tmp2,
                    box_pairs=[(ALL, ALL2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    if rc == 'need_rc':
        ax0.set_ylabel('read counts*M/nt',fontsize=20) #log10(read counts*M/nt)
    else:
        ax0.set_ylabel('sites counts*M/nt',fontsize=20)

    ax1.set_title(title_map_rna+'_UTR5\n'+title_map_gene1+' N='+str(len(all_utr5))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr52))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    plt.ylim(0, 1000000*max([utr5_max, utr52_max, cds_max, cds2_max, utr3_max, utr32_max]))
    tmp2 = pd.DataFrame({UTR5:pd.Series(utr5_density1*1000000),UTR52:pd.Series(utr5_density2*1000000)})
    my_pal = {UTR5:'#FFA07A', UTR52:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax1)
    #ax = sns.violinplot(data=tmp2, scale="count", split=True, inner="quartile", palette="coolwarm")
    add_stat_annotation(ax1,data=tmp2,
                    box_pairs=[(UTR5, UTR52)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    
    ax2.set_title(title_map_rna+'_CDS\n'+title_map_gene1+' N='+str(len(all_cds))+'\n'+
                 title_map_gene2+' N='+str(len(all_cds2))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({CDS:pd.Series(cds_density1*1000000),CDS2:pd.Series(cds_density2*1000000)})
    my_pal = {CDS:'#FFA07A', CDS2:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax2)    
    add_stat_annotation(ax2,data=tmp2,
                    box_pairs=[(CDS, CDS2)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    
    ax3.set_title(title_map_rna+'_UTR3\n'+title_map_gene1+' N='+str(len(all_utr3))+'\n'+
                 title_map_gene2+' N='+str(len(all_utr32))+'\n\n\n')#+title_map_gene[str(gene)]+'\n #UTR5='+str(len(all_utr5))+' #CDS='+str(len(all_cds))+' #UTR3='+str(len(all_utr3)), fontsize=20) #Algorithm_name+'_'+
    
    tmp2 = pd.DataFrame({UTR3:pd.Series(utr3_density1*1000000),UTR32:pd.Series(utr3_density2*1000000)})
    my_pal = {UTR3:'#FFA07A', UTR32:'#98F898'}
    sns.boxplot(data=tmp2, showfliers = False, width=0.3, palette=my_pal,showmeans=True, ax=ax3)
    add_stat_annotation(ax3,data=tmp2,
                    box_pairs=[(UTR3, UTR32)],
                    test='t-test_welch', text_format='full', loc='outside', verbose=2)
    tmp2 = pd.DataFrame({title_map_gene1+' UTR5':pd.Series(all_utr5),title_map_gene2+' UTR5':pd.Series(all_utr52),
                        title_map_gene1+' CDS':pd.Series(all_cds),title_map_gene2+' CDS':pd.Series(all_cds2),
                        title_map_gene1+' UTR3':pd.Series(all_utr3),title_map_gene2+' UTR3':pd.Series(all_utr32)})
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
    
    plt.tight_layout()
    plt.savefig(res_path.replace('.png', '_non-log_T.png'))
    plt.clf()
    plt.close()
    del tmp2
    gc.collect()

def main_plot_g22_new(all_density1, all_density2, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc, FILTER):
    # preprocess for plot
    total_all1 = [0.000001 for i in all_density1 if i == 0] + [i for i in all_density1 if i != 0]
    total_all1 = [math.log10(i*1000000) for i in total_all1]

    total_all2 = [0.000001 for i in all_density2 if i == 0] + [i for i in all_density2 if i != 0]
    total_all2 = [math.log10(i*1000000) for i in total_all2]
    
    single_plot_g22_new(total_all1, total_all2, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc, FILTER, all_density1, all_density2)
    
def main_plot_g22(all_density1, all_density2, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2, FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2, rc, FILTER):
    # all_mRNA
    total_all1 = [0.000001 for i in all_density1 if i == 0] + [i for i in all_density1 if i != 0]
    total_all1 = [math.log10(i*1000000) for i in total_all1]

    total_all2 = [0.000001 for i in all_density2 if i == 0] + [i for i in all_density2 if i != 0]
    total_all2 = [math.log10(i*1000000) for i in total_all2]
    
    # preprocess for plot
    total_utr5 = [0.000001 for i in utr5_density1 if i == 0] + [i for i in utr5_density1 if i != 0]
    total_utr5 = [math.log10(i*1000000) for i in total_utr5]

    total_cds = [0.000001 for i in cds_density1 if i == 0] + [i for i in cds_density1 if i != 0]
    total_cds = [math.log10(i*1000000) for i in total_cds]

    total_utr3 = [0.000001 for i in utr3_density1 if i == 0] + [i for i in utr3_density1 if i != 0]
    total_utr3 = [math.log10(i*1000000) for i in total_utr3]

    total_utr52 = [0.000001 for i in utr5_density2 if i == 0] + [i for i in utr5_density2 if i != 0]
    total_utr52 = [math.log10(i*1000000) for i in total_utr52]

    total_cds2 = [0.000001 for i in cds_density2 if i == 0] + [i for i in cds_density2 if i != 0]
    total_cds2 = [math.log10(i*1000000) for i in total_cds2]

    total_utr32 = [0.000001 for i in utr3_density2 if i == 0] + [i for i in utr3_density2 if i != 0]
    total_utr32 = [math.log10(i*1000000) for i in total_utr32]
    
    single_plot_g22(total_all1, total_all2, total_utr5, total_cds, total_utr3, total_utr52, total_cds2, total_utr32,FILE_NAME, title_map_rna, gene1, title_map_gene1, gene2, title_map_gene2,
                    rc, FILTER, all_density1, all_density2, utr5_density1, cds_density1, utr3_density1, utr5_density2, cds_density2, utr3_density2)    

