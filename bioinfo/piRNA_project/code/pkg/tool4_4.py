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
import warnings
warnings.filterwarnings("ignore")

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

def tool4_4(data):
    all_mRNA = pd.read_csv(os.path.abspath(__file__+ '/../../../../data/mRNA_WS275_IDtoName.csv'))
    gene_path = data['gene_path']
    title_map_gene = data['title_map_gene']
    title_map_rna = data['title_map_rna']

    for ANALYZE_TYPE in data['ANALYZE_TYPE']:
        #print('===={}===='.format(ANALYZE_TYPE))
        if ANALYZE_TYPE == 'TARGET_SCORE':
            for FILE_NAME1,FILE_NAME2 in zip(data['COMPARE_DATA1'], data['COMPARE_DATA2']):
                for rc in data['READ_COUNT_TYPE']:
                    for gene in data['ANALYZE_GROUP']:
                        for FILTER in data['FOLD_LEN_FILTER']:
                            gene_name = check_new(data, gene)
                            utr5_density1, cds_density1, utr3_density1 = main_filter_mRNA(FILE_NAME1, rc, gene_path, gene, gene_name, FILTER)
                            utr5_density2, cds_density2, utr3_density2 = main_filter_mRNA(FILE_NAME2, rc, gene_path, gene, gene_name, FILTER)
                            plot_ratio_pseuodo_count(FILE_NAME1, title_map_rna[FILE_NAME1], FILE_NAME2, title_map_rna[FILE_NAME2], utr5_density1 ,utr5_density2, cds_density1, cds_density2, utr3_density1, utr3_density2, gene, title_map_gene[gene], FILTER, rc)
                            #plot_scatter(FILE_NAME1, title_map_rna[FILE_NAME1], FILE_NAME2, title_map_rna[FILE_NAME2], utr5_density1 ,utr5_density2, cds_density1, cds_density2, utr3_density1, utr3_density2, gene, title_map_gene[gene], FILTER, rc)

        elif ANALYZE_TYPE == 'G22':
            #two_set_group_dict = {}
            for FILE_NAME1,FILE_NAME2 in zip(data['COMPARE_DATA1'], data['COMPARE_DATA2']):
                for rc in data['READ_COUNT_TYPE']:
                    for gene in data['ANALYZE_GROUP']:
                        for FILTER in data['FOLD_LEN_FILTER']:
                            gene_name = check_new(data, gene)

                            #### all region
                            #all_density1 = main_filter_mRNA_g22_new(FILE_NAME1, rc, gene_path, gene, gene_name, FILTER)
                            #all_density2 = main_filter_mRNA_g22_new(FILE_NAME2, rc, gene_path, gene, gene_name, FILTER)
                            #plot_scatter_g22_new(FILE_NAME1, title_map_rna[FILE_NAME1], FILE_NAME2, title_map_rna[FILE_NAME2], all_density1, all_density2, gene, title_map_gene[gene], FILTER, rc)

                            #### 3 regions
                            all_density1, utr5_density1, cds_density1, utr3_density1 = main_filter_mRNA_g22(FILE_NAME1, rc, gene_path, gene, gene_name, FILTER)
                            all_density2, utr5_density2, cds_density2, utr3_density2 = main_filter_mRNA_g22(FILE_NAME2, rc, gene_path, gene, gene_name, FILTER)
                            plot_ratio_pseuodo_count_g22(FILE_NAME1, title_map_rna[FILE_NAME1], FILE_NAME2, title_map_rna[FILE_NAME2], all_density1, all_density2, utr5_density1 ,utr5_density2, cds_density1, cds_density2, utr3_density1, utr3_density2, gene, title_map_gene[gene], FILTER, rc, data['img_title'])
                            #plot_scatter_g22(FILE_NAME1, title_map_rna[FILE_NAME1], FILE_NAME2, title_map_rna[FILE_NAME2], all_density1, all_density2, utr5_density1 ,utr5_density2, cds_density1, cds_density2, utr3_density1, utr3_density2, gene, title_map_gene[gene], FILTER, rc)
                            #two_set_group_dict.setdefault(str(gene)+'_'+FILE_NAME1, [utr5_density1, cds_density1, utr3_density1])
                            #two_set_group_dict.setdefault(str(gene)+'_'+FILE_NAME2, [utr5_density2, cds_density2, utr3_density2])
                                    
                            # two group, two dataset FOLD CHANGE compare
                            #GROUP_NAME1 = ANALYZE_GROUP[0]
                            #GROUP_NAME2 = ANALYZE_GROUP[1]
                            #plot_ratio_pseuodo_count_g22_2data_3region(GROUP_NAME1, GROUP_NAME2, FILE_NAME1, FILE_NAME2, two_set_group_dict, FILTER, rc)                 
                            #plot_scatter_2data_g22_3region(GROUP_NAME1, GROUP_NAME2, FILE_NAME1, FILE_NAME2, two_set_group_dict, FILTER, rc)


def main_filter_mRNA(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+ '/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
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

def plot_ratio_pseuodo_count(FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, d51, d52, dc1, dc2, d31, d32, gene, title_map_gene, FILTER, rc):
    #### with 0
    ALPHA1 = min(find_alpha(d51), find_alpha(dc1), find_alpha(d31),
               find_alpha(d52), find_alpha(dc2), find_alpha(d32))#非零的所有值的最小值*0.0001
                
    #res_path = os.path.abspath(__file__+ '/../../../usetRNA_target_score_'+rc+'_output/fold_change_fig/FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/fold_change_fig/compare/FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'_G'+str(gene)+'_L'+str(FILTER)+'_'+rc+'.png') 
    
    z_dict = {1:['_log_U.png', 'Mann-Whitney'],
             2:['_log_T.png', 't-test_welch'],
             3:['_non-log_U.png', 'Mann-Whitney'],
             4:['_non-log_T.png', 't-test_welch']} 
    for z in range(1,5):
        
        x1 = [i+ALPHA1 for i in d51 ] #control
        x2 = [i+ALPHA1 for i in dc1 ] #control
        x3 = [i+ALPHA1 for i in d31 ] #control

        y1 = [i+ALPHA1 for i in d52 ] #mutant
        y2 = [i+ALPHA1 for i in dc2 ] #mutant
        y3 = [i+ALPHA1 for i in d32 ] #mutant
        
        if z == 1 or z == 2:
            plot_value1 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)] # 5'
            plot_value3 = [math.log(mut/con, 2) for mut,con in zip(y2,x2)] # cds
            plot_value5 = [math.log(mut/con, 2) for mut,con in zip(y3,x3)] # 3'
            plot_max = max(max(plot_value1),max(plot_value3),max(plot_value5))
            plot_min = min(min(plot_value1),min(plot_value3),min(plot_value5))
        elif z == 3 or z == 4:
            plot_value1 = [mut/con for mut,con in zip(y1,x1)]
            plot_value3 = [mut/con for mut,con in zip(y2,x2)]
            plot_value5 = [mut/con for mut,con in zip(y3,x3)]
            plot_max = max(max(plot_value1),max(plot_value3),max(plot_value5))
            plot_min = min(min(plot_value1),min(plot_value3),min(plot_value5))
            
        #### without 0
        d51_no_b = []
        d52_no_b = []

        dc1_no_b = []
        dc2_no_b = []

        d31_no_b = []
        d32_no_b = []
        for i,j in zip(d51 ,d52):
            if i != 0 and j != 0:
                d51_no_b.append(i)
                d52_no_b.append(j)

        for i,j in zip(dc1 ,dc2):
            if i != 0 and j != 0:
                dc1_no_b.append(i)
                dc2_no_b.append(j)

        for i,j in zip(d31 ,d32):
            if i != 0 and j != 0:
                d31_no_b.append(i)
                d32_no_b.append(j)
        ALPHA2 = min(find_alpha(d51_no_b), find_alpha(dc1_no_b), find_alpha(d31_no_b),
                   find_alpha(d52_no_b), find_alpha(dc2_no_b), find_alpha(d32_no_b))#非零的所有值的最小值*0.0001

        x1 = [i+ALPHA2 for i in d51_no_b ] #control
        x2 = [i+ALPHA2 for i in dc1_no_b ] #control
        x3 = [i+ALPHA2 for i in d31_no_b ] #control

        y1 = [i+ALPHA2 for i in d52_no_b ] #mutant
        y2 = [i+ALPHA2 for i in dc2_no_b ] #mutant
        y3 = [i+ALPHA2 for i in d32_no_b ] #mutant
    
        if z == 1 or z == 2:
            plot_value2 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)] # 5'
            plot_value4 = [math.log(mut/con, 2) for mut,con in zip(y2,x2)] # cds
            plot_value6 = [math.log(mut/con, 2) for mut,con in zip(y3,x3)] # 3'
            plot_max = max(max(plot_value2),max(plot_value4),max(plot_value6))
            plot_min = min(min(plot_value2),min(plot_value4),min(plot_value6))
        elif z == 3 or z == 4:
            plot_value2 = [mut/con for mut,con in zip(y1,x1)]
            plot_value4 = [mut/con for mut,con in zip(y2,x2)]
            plot_value6 = [mut/con for mut,con in zip(y3,x3)]
            plot_max = max(max(plot_value2),max(plot_value4),max(plot_value6))
            plot_min = min(min(plot_value2),min(plot_value4),min(plot_value6))

        #### plot
        f, ax = plt.subplots(figsize=(15, 10))
        tmp2 = pd.DataFrame({'UTR5':pd.Series(plot_value1), 'UTR5\nWITHOUT 0':pd.Series(plot_value2),
                             'CDS':pd.Series(plot_value3), 'CDS\nWITHOUT 0':pd.Series(plot_value4),
                             'UTR3':pd.Series(plot_value5), 'UTR3\nWITHOUT 0':pd.Series(plot_value6)})

        my_pal = {'UTR5':'#FFA07A', 'UTR5\nWITHOUT 0':'#FF7F50',
                  'CDS':'#98F898', 'CDS\nWITHOUT 0':'#90EE90',
                  'UTR3':'#BA55D3','UTR3\nWITHOUT 0':'#8A2BE2'}

        #### plot with blue points ####
        mean_5 = round(np.mean(plot_value1),3)
        median_5 = round(np.median(plot_value1),3)
        mean_cds = round(np.mean(plot_value3),3)
        median_cds = round(np.median(plot_value3),3)
        mean_3 = round(np.mean(plot_value5),3)
        median_3 = round(np.median(plot_value5),3)
        UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
        CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
        UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)

        tmp2 = pd.DataFrame({UTR5:pd.Series(plot_value1), CDS:pd.Series(plot_value3), UTR3:pd.Series(plot_value5)})
        my_pal = {UTR5:'#FFA07A',CDS:'#98F898',UTR3:'#BA55D3'}
        box1 = sns.boxplot(data=tmp2,showfliers = False, width=0.3, palette=my_pal, showmeans=True)#.set(xlabel='region',ylabel='log2 ('+s2+' / '+s1+')')
        add_stat_annotation(box1,data=tmp2,
                        box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                        test=z_dict[z][1], text_format='full', loc='outside', verbose=0,fontsize=16)
        # plot point
        plt.xticks(fontsize = 22) 
        plt.yticks(fontsize = 22)
        box1.axes.set_title(title_map_rna1 +' v.s. '+title_map_rna2+ ' '+title_map_gene+
                     '\n #UTR5='+str(len(plot_value1))+'  #CDS='+str(len(plot_value3))+'  #UTR3='+str(len(plot_value5))+
                      '\n#WITH 0 α='+str(round(ALPHA1,7))+'\n\n\n\n\n\n', fontsize=20)
        box1.set_xticklabels(box1.get_xmajorticklabels(), fontsize = 16)
        if z == 1 or z == 2:
            box1.set_ylabel('log2 ('+title_map_rna2+'+α / '+title_map_rna1+'+α)',fontsize=20)
        elif z == 3 or z == 4:
            box1.set_ylabel(title_map_rna2+'+α / '+title_map_rna1+'+α',fontsize=20)
        ax.axhline(y=0, c='black' ,linestyle="--")

        #tmp2.to_csv(res_path.replace('.png', z_dict[z][1]), index=False)
        plt.tight_layout()
        plt.savefig(res_path.replace('.png', z_dict[z][0]))

        #### plot without blue points ####

        f, ax = plt.subplots(figsize=(15, 10))

        mean_5 = round(np.mean(plot_value2),3)
        median_5 = round(np.median(plot_value2),3)
        mean_cds = round(np.mean(plot_value4),3)
        median_cds = round(np.median(plot_value4),3)
        mean_3 = round(np.mean(plot_value6),3)
        median_3 = round(np.median(plot_value6),3)
        UTR5 = 'UTR5\navg:{}\nmedian:{}'.format(mean_5, median_5)
        CDS = 'CDS\navg:{}\nmedian:{}'.format(mean_cds, median_cds)
        UTR3 = 'UTR3\navg:{}\nmedian:{}'.format(mean_3, median_3)

        tmp2 = pd.DataFrame({UTR5:pd.Series(plot_value2), CDS:pd.Series(plot_value4), UTR3:pd.Series(plot_value6)})
        my_pal = {UTR5:'#FFA07A',CDS:'#98F898',UTR3:'#BA55D3'}
        box1 = sns.boxplot(data=tmp2,showfliers = False, width=0.3, palette=my_pal, showmeans=True)#.set(xlabel='region',ylabel='log2 ('+s2+' / '+s1+')')
        add_stat_annotation(box1,data=tmp2,
                        box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                        test=z_dict[z][1], text_format='full', loc='outside', verbose=0,fontsize=16)
        # plot point
        plt.xticks(fontsize = 22) 
        plt.yticks(fontsize = 22)
        box1.axes.set_title(title_map_rna1 +' v.s. '+title_map_rna2+ ' '+title_map_gene+
                     '\n #UTR5='+str(len(plot_value1))+'  #CDS='+str(len(plot_value3))+'  #UTR3='+str(len(plot_value5))+
                      '\n#WITHOUT 0 α='+str(round(ALPHA1,7))+'\n\n\n\n\n\n', fontsize=20)
        box1.set_xticklabels(box1.get_xmajorticklabels(), fontsize = 16)
        if z == 1 or z == 2:
            box1.set_ylabel('log2 ('+title_map_rna2+'+α / '+title_map_rna1+'+α)',fontsize=20)
        elif z == 3 or z == 4:
            box1.set_ylabel(title_map_rna2+'+α / '+title_map_rna1+'+α',fontsize=20)
            
        ax.axhline(y=0, c='black' ,linestyle="--")

        #tmp2.to_csv(res_path.replace('.png', 'WITHOUT_B.csv'), index=False)
        plt.tight_layout()
        plt.savefig(res_path.replace('.png', '_WITHOUT0'+z_dict[z][0]))

def plot_scatter(FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, d51, d52, dc1, dc2, d31, d32, gene, title_map_gene, FILTER, rc):
    #res_path = os.path.abspath(__file__+ '/../../../usetRNA_target_score_'+rc+'_output/fold_change_fig/Sca_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/fold_change_fig/compare/Sca_'+FILE_NAME1+'_'+FILE_NAME2+'_G'+str(gene)+'_L'+str(FILTER)+'_'+rc+'.png')
    x1_f = []
    x2_f = []
    x3_f = []
    y1_f = []
    y2_f = []
    y3_f = []
    blue_p_x1 = 0
    blue_p_y1 = 0
    blue_p_x2 = 0
    blue_p_y2 = 0
    blue_p_x3 = 0
    blue_p_y3 = 0
    blue_both1 = 0
    blue_both2 = 0
    blue_both3 = 0
    for i,j in zip(d51 ,d52):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y1 += 1
                else:
                    blue_both1 += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x1 += 1
            x1_f.append(math.log(i, 2))
            y1_f.append(math.log(j, 2))
    for i,j in zip(dc1, dc2):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y2 += 1
                else:
                    blue_both2 += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x2 += 1
            x2_f.append(math.log(i, 2))
            y2_f.append(math.log(j, 2))
    for i,j in zip(d31, d32):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y3 += 1
                else:
                    blue_both3 += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x3 += 1
            x3_f.append(math.log(i, 2))
            y3_f.append(math.log(j, 2))
            
    x1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in d51]
    x2 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in dc1]
    x3 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in d31]
    
    y1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in d52]
    y2 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in dc2]
    y3 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in d32]


    #f, ax = plt.subplots(figsize=(8, 8))
    fig , ((ax1,ax2,ax3)) = plt.subplots(1,3,figsize=(23,10))
    fig.suptitle(title_map_rna1 +' v.s. '+title_map_rna2+ ' '+title_map_gene+
                 '\n #UTR5='+str(len(x1))+' #CDS='+str(len(x2))+' #UTR3='+str(len(x3)), fontsize=14,y=1)
    if rc == 'need_rc':
        ylab = 'read count'
    else:
        ylab = 'site'
    fig.text(0.5, 0.05, 'log2('+ylab+') '+title_map_rna1, ha='center',size=14)
    fig.text(0.07, 0.5, 'log2('+ylab+') '+title_map_rna2, va='center', rotation='vertical',size=14)
    
    ax1.scatter(x1, y1, s=1, c='#AE0000',alpha=0.5)
    ax1.scatter(x1_f, y1_f, s=1, c='#0066FF')
    
    ax1.set(xlim=ax1.get_xlim(), ylim=ax1.get_ylim())
    diag_x = list(ax1.get_xlim())
    diag_y = list(ax1.get_ylim())
    ax1.plot(diag_x, diag_y, c='black')
    ax1.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax1.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    #print(len(x1),len(y1),len(x1_f),len(y1_f))
    ax1.axes.set_title('UTR5\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x1)-len(x1_f),
                                                                                           title_map_rna1,blue_p_x1-blue_both1,
                                                                                           title_map_rna2,blue_p_y1,
                                                                                            blue_both1))
    
    ax2.scatter(x2, y2, s=1, c='#AE0000',alpha=0.5)
    ax2.scatter(x2_f, y2_f, s=1, c='#0066FF')
    ax2.set(xlim=ax2.get_xlim(), ylim=ax2.get_ylim())
    diag_x = list(ax2.get_xlim())
    diag_y = list(ax2.get_ylim())
    ax2.plot(diag_x, diag_y, c='black')
    ax2.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax2.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    #print(len(x2),len(y2),len(x2_f),len(y2_f))
    ax2.axes.set_title('CDS\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x2)-len(x2_f),
                                                                                           title_map_rna1,blue_p_x2-blue_both2,
                                                                                           title_map_rna2,blue_p_y2,
                                                                                            blue_both2))
    
    ax3.scatter(x3, y3, s=1, c='#AE0000',alpha=0.5)
    ax3.scatter(x3_f, y3_f, s=1, c='#0066FF')
    plt.legend(['valid value','replace 0 by 0.00001'], bbox_to_anchor=(1.05, 1), loc='upper left')
    ax3.set(xlim=ax3.get_xlim(), ylim=ax3.get_ylim())
    diag_x = list(ax3.get_xlim())
    diag_y = list(ax3.get_ylim())
    ax3.plot(diag_x, diag_y, c='black')
    ax3.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax3.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    #print(len(x3),len(y3),len(x3_f),len(y3_f))
    ax3.axes.set_title('UTR3\nRED POINT={}\nBLUE BLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x3)-len(x3_f),
                                                                                           title_map_rna1,blue_p_x3-blue_both3,
                                                                                           title_map_rna2,blue_p_y3,
                                                                                            blue_both3))
    
    plt.savefig(res_path)
    tmp2 = pd.DataFrame({'UTR5_1':d51, 'UTR5_2':d52,
                         'CDS_1':dc1, 'CDS_2':dc2,
                         'UTR3_1':d31, 'UTR3_2':d32,
    })
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)

def find_alpha(l):
    l2 = [i for i in l if i != 0]
    return min(l2) #非零當中最小值的萬分之一倍

def main_filter_mRNA_g22(FILE_NAME, rc, gene_path, gene, title_map_gene, FILTER):
    #mRNA list
    data = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    group = add_two_mRNA_list(data, gene_path, gene, title_map_gene)
    #print('group'+str(gene)+' mRNA=',len(group))
    
    # FILTER LENGTH
    group['all_rc'] = group['count5']+group['countcds']+group['count3']
    group['all_len'] = group['len5']+group['lencds']+group['len3']
    all_group = group[(group['all_len']>=FILTER)]
    all_density = all_group['all_rc'] / all_group['all_len']

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
    data = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    #data['Gene ID'] = data['Gene ID'].apply(lambda x:x.split('=')[1])
    all_group = add_two_mRNA_list(data, gene_path, gene, title_map_gene)
    #print('group'+str(gene)+' mRNA=',len(all_group))
    
    # FILTER LENGTH
    all_group['all_rc'] = all_group['count5']+all_group['countcds']+all_group['count3']
    all_group['all_len'] = all_group['len5']+all_group['lencds']+all_group['len3']
    all_group = all_group[(all_group['all_len']>=FILTER)]
    all_density = all_group['all_rc'] / all_group['all_len']
    
    #print(len(all_density))
    
    return all_density

def plot_ratio_pseuodo_count_g22(FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, dall1, dall2, d51, d52, dc1, dc2, d31, d32, gene, title_map_gene, FILTER, rc, img_title):

    #### with 0 ####
    ALPHA1 = min(find_alpha(d51), find_alpha(dc1), find_alpha(d31),
                find_alpha(d52), find_alpha(dc2), find_alpha(d32))#非零的所有值的最小值*0.0001
                
    #res_path = os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/fold_change_fig/FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/fold_change_fig/compare/FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    z_dict = {1:['_log_U.png', 'Mann-Whitney'],} 
    for z in [1]:
        
        x1 = [i+ALPHA1 for i in d51 ] #control
        x2 = [i+ALPHA1 for i in dc1 ] #control
        x3 = [i+ALPHA1 for i in d31 ] #control

        y1 = [i+ALPHA1 for i in d52 ] #mutant
        y2 = [i+ALPHA1 for i in dc2 ] #mutant
        y3 = [i+ALPHA1 for i in d32 ] #mutant
        
        if z == 1 or z == 2:
            plot_value1 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)] # 5'
            plot_value3 = [math.log(mut/con, 2) for mut,con in zip(y2,x2)] # cds
            plot_value5 = [math.log(mut/con, 2) for mut,con in zip(y3,x3)] # 3'
            plot_max = max(max(plot_value1),max(plot_value3),max(plot_value5))
            plot_min = min(min(plot_value1),min(plot_value3),min(plot_value5))
        elif z == 3 or z == 4:
            plot_value1 = [mut/con for mut,con in zip(y1,x1)]
            plot_value3 = [mut/con for mut,con in zip(y2,x2)]
            plot_value5 = [mut/con for mut,con in zip(y3,x3)]
            plot_max = max(max(plot_value1),max(plot_value3),max(plot_value5))
            plot_min = min(min(plot_value1),min(plot_value3),min(plot_value5))
            
        #### without 0
        d51_no_b = []
        d52_no_b = []

        dc1_no_b = []
        dc2_no_b = []

        d31_no_b = []
        d32_no_b = []

        for i,j in zip(d51 ,d52):
            if i != 0 and j != 0:
                d51_no_b.append(i)
                d52_no_b.append(j)

        for i,j in zip(dc1 ,dc2):
            if i != 0 and j != 0:
                dc1_no_b.append(i)
                dc2_no_b.append(j)

        for i,j in zip(d31 ,d32):
            if i != 0 and j != 0:
                d31_no_b.append(i)
                d32_no_b.append(j)
        ALPHA2 = min(find_alpha(d51_no_b), find_alpha(dc1_no_b), find_alpha(d31_no_b),
                    find_alpha(d52_no_b), find_alpha(dc2_no_b), find_alpha(d32_no_b))#非零的所有值的最小值*0.0001

        x1 = [i+ALPHA2 for i in d51_no_b ] #control
        x2 = [i+ALPHA2 for i in dc1_no_b ] #control
        x3 = [i+ALPHA2 for i in d31_no_b ] #control

        y1 = [i+ALPHA2 for i in d52_no_b ] #mutant
        y2 = [i+ALPHA2 for i in dc2_no_b ] #mutant
        y3 = [i+ALPHA2 for i in d32_no_b ] #mutant
    
        if z == 1 or z == 2:
            plot_value2 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)] # 5'
            plot_value4 = [math.log(mut/con, 2) for mut,con in zip(y2,x2)] # cds
            plot_value6 = [math.log(mut/con, 2) for mut,con in zip(y3,x3)] # 3'
            plot_max = max(max(plot_value2),max(plot_value4),max(plot_value6))
            plot_min = min(min(plot_value2),min(plot_value4),min(plot_value6))
        elif z == 3 or z == 4:
            plot_value2 = [mut/con for mut,con in zip(y1,x1)]
            plot_value4 = [mut/con for mut,con in zip(y2,x2)]
            plot_value6 = [mut/con for mut,con in zip(y3,x3)]
            plot_max = max(max(plot_value2),max(plot_value4),max(plot_value6))
            plot_min = min(min(plot_value2),min(plot_value4),min(plot_value6))

        #### plot
        f, ax = plt.subplots(figsize=(8, 8))
        tmp2 = pd.DataFrame({'UTR5':pd.Series(plot_value1), 'UTR5\nWITHOUT 0':pd.Series(plot_value2),
                             'CDS':pd.Series(plot_value3), 'CDS\nWITHOUT 0':pd.Series(plot_value4),
                             'UTR3':pd.Series(plot_value5), 'UTR3\nWITHOUT 0':pd.Series(plot_value6)})

        #### plot with blue points ####
        UTR5 = "5'UTR"
        CDS = "CDS"
        UTR3 = "3'UTR"

        tmp2 = pd.DataFrame({UTR5:pd.Series(plot_value1), CDS:pd.Series(plot_value3), UTR3:pd.Series(plot_value5)})
        box1 = sns.boxplot(data=tmp2,showfliers = False, width=0.3, palette="coolwarm")#.set(xlabel='region',ylabel='log2 ('+s2+' / '+s1+')')
        add_stat_annotation(box1,data=tmp2,
                        box_pairs=[(UTR5,CDS),(UTR5,UTR3),(CDS,UTR3)],
                        test=z_dict[z][1], text_format='full', loc='outside', verbose=0, fontsize=16)
        # plot point
        plt.xticks(fontsize = 22) 
        plt.yticks(fontsize = 22)
        box1.axes.set_title(img_title+'\n\n\n\n', fontsize=25)
        box1.set_xticklabels(box1.get_xmajorticklabels(), fontsize = 16)
        if z == 1 or z == 2:
            box1.set_ylabel('$\mathregular{log_2}$ ('+title_map_rna2+' / '+title_map_rna1+')', fontsize=20)
        elif z == 3 or z == 4:
            box1.set_ylabel(title_map_rna2+' / '+title_map_rna1, fontsize=20)
        ax.axhline(y=0, c='black' ,linestyle="--")

        #tmp2.to_csv(res_path.replace('.png', z_dict[z][1]), index=False)
        plt.tight_layout()
        plt.savefig(res_path.replace('.png', z_dict[z][0]))

        

def plot_scatter_g22(FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, dall1, dall2, d51, d52, dc1, dc2, d31, d32, gene, title_map_gene, FILTER, rc):
    #res_path = os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/fold_change_fig/Sca_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/fold_change_fig/compare/Sca_'+FILE_NAME1+'_'+FILE_NAME2+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    x0_f = []
    x1_f = []
    x2_f = []
    x3_f = []
    y0_f = []
    y1_f = []
    y2_f = []
    y3_f = []
    blue_p_x0 = 0
    blue_p_y0 = 0
    blue_p_x1 = 0
    blue_p_y1 = 0
    blue_p_x2 = 0
    blue_p_y2 = 0
    blue_p_x3 = 0
    blue_p_y3 = 0
    blue_both0 = 0
    blue_both1 = 0
    blue_both2 = 0
    blue_both3 = 0
    for i,j in zip(dall1 ,dall2):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y0 += 1
                else:
                    blue_both0 += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x0 += 1
            x0_f.append(math.log(i, 2))
            y0_f.append(math.log(j, 2))
    for i,j in zip(d51 ,d52):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y1 += 1
                else:
                    blue_both1 += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x1 += 1
            x1_f.append(math.log(i, 2))
            y1_f.append(math.log(j, 2))
    for i,j in zip(dc1, dc2):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y2 += 1
                else:
                    blue_both2 += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x2 += 1
            x2_f.append(math.log(i, 2))
            y2_f.append(math.log(j, 2))
    for i,j in zip(d31, d32):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y3 += 1
                else:
                    blue_both3 += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x3 += 1
            x3_f.append(math.log(i, 2))
            y3_f.append(math.log(j, 2))
            
    x0 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in dall1]
    x1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in d51]
    x2 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in dc1]
    x3 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in d31]
    
    y0 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in dall2]
    y1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in d52]
    y2 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in dc2]
    y3 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in d32]


    #f, ax = plt.subplots(figsize=(8, 8))
    fig , ((ax0,ax1,ax2,ax3)) = plt.subplots(1,4,figsize=(32,10))
    fig.suptitle(title_map_rna1 +' v.s. '+title_map_rna2+ ' '+title_map_gene+
                 '\n #UTR5='+str(len(x1))+' #CDS='+str(len(x2))+' #UTR3='+str(len(x3)), fontsize=14,y=1)
    if rc == 'need_rc':
        ylab = 'read count'
    else:
        ylab = 'site'
    fig.text(0.5, 0.05, 'log2('+ylab+') '+title_map_rna1, ha='center',size=14)
    fig.text(0.07, 0.5, 'log2('+ylab+') '+title_map_rna2, va='center', rotation='vertical',size=14)
    
    ax0.scatter(x0, y0, s=1, c='#AE0000',alpha=0.5)
    ax0.scatter(x0_f, y0_f, s=1, c='#0066FF')
    ax0.set(xlim=ax0.get_xlim(), ylim=ax0.get_ylim())
    diag_x = list(ax0.get_xlim())
    diag_y = list(ax0.get_ylim())
    ax0.plot(diag_x, diag_y, c='black')
    ax0.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax0.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    #print(len(x0),len(y0),len(x0_f),len(y0_f))
    ax0.axes.set_title('ALL\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x0)-len(x0_f),
                                                                                           title_map_rna1,blue_p_x0-blue_both0,
                                                                                           title_map_rna2,blue_p_y0,
                                                                                           blue_both0))

    ax1.scatter(x1, y1, s=1, c='#AE0000',alpha=0.5)
    ax1.scatter(x1_f, y1_f, s=1, c='#0066FF')
    ax1.set(xlim=ax1.get_xlim(), ylim=ax1.get_ylim())
    diag_x = list(ax1.get_xlim())
    diag_y = list(ax1.get_ylim())
    ax1.plot(diag_x, diag_y, c='black')
    ax1.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax1.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    #print(len(x1),len(y1),len(x1_f),len(y1_f))
    ax1.axes.set_title('UTR5\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x1)-len(x1_f),
                                                                                           title_map_rna1,blue_p_x1-blue_both1,
                                                                                           title_map_rna2,blue_p_y1,
                                                                                            blue_both1))
    
    ax2.scatter(x2, y2, s=1, c='#AE0000',alpha=0.5)
    ax2.scatter(x2_f, y2_f, s=1, c='#0066FF')
    ax2.set(xlim=ax2.get_xlim(), ylim=ax2.get_ylim())
    diag_x = list(ax2.get_xlim())
    diag_y = list(ax2.get_ylim())
    ax2.plot(diag_x, diag_y, c='black')
    ax2.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax2.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    #print(len(x2),len(y2),len(x2_f),len(y2_f))
    ax2.axes.set_title('CDS\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x2)-len(x2_f),
                                                                                           title_map_rna1,blue_p_x2-blue_both2,
                                                                                           title_map_rna2,blue_p_y2,
                                                                                            blue_both2))
    
    ax3.scatter(x3, y3, s=1, c='#AE0000',alpha=0.5)
    ax3.scatter(x3_f, y3_f, s=1, c='#0066FF')
    plt.legend(['valid value','replace 0 by 0.00001'], bbox_to_anchor=(1.05, 1), loc='upper left')
    ax3.set(xlim=ax3.get_xlim(), ylim=ax3.get_ylim())
    diag_x = list(ax3.get_xlim())
    diag_y = list(ax3.get_ylim())
    ax3.plot(diag_x, diag_y, c='black')
    ax3.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax3.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    #print(len(x3),len(y3),len(x3_f),len(y3_f))
    ax3.axes.set_title('UTR3\nRED POINT={}\nBLUE BLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x3)-len(x3_f),
                                                                                           title_map_rna1,blue_p_x3-blue_both3,
                                                                                           title_map_rna2,blue_p_y3,
                                                                                            blue_both3))
    plt.savefig(res_path)
    tmp2 = pd.DataFrame({'UTR5_1':d51, 'UTR5_2':d52,
                         'CDS_1':dc1, 'CDS_2':dc2,
                         'UTR3_1':d31, 'UTR3_2':d32,
    })
    #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)

def plot_scatter_g22_new(FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, dall1, dall2, gene, title_map_gene, FILTER, rc):
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/fold_change_fig/compare/Sca_ALL_'+FILE_NAME1+'_'+FILE_NAME2+'_G'+str(gene)+'_L'+str(FILTER)+'.png')
    # scatter
    
    x1_f = []
    y1_f = []
    blue_p_x = 0
    blue_p_y = 0
    blue_both = 0
    
    for i,j in zip(dall1 ,dall2):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y += 1
                else:
                    blue_both += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x += 1
            x1_f.append(math.log(i, 2))
            y1_f.append(math.log(j, 2))
            
    x1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in dall1]
    y1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in dall2]
    
    fig , ((ax1,ax2,ax3)) = plt.subplots(1,3,figsize=(22,9))
    fig.suptitle(title_map_rna1 +' v.s. '+title_map_rna2+ ' '+title_map_gene+
                 '\n #FULL GENE='+str(len(x1)),fontsize=15,y=1)
    if rc == 'need_rc':
        ylab = 'read count'
    else:
        ylab = 'site'
    ax1.set_xlabel('log2('+ylab+') '+title_map_rna1)
    ax1.set_ylabel('log2('+ylab+') '+title_map_rna2)
    
    ax1.scatter(x1, y1, s=1, c='#AE0000',alpha=0.5)
    ax1.scatter(x1_f, y1_f, s=1, c='#0066FF')
    ax1.legend(['valid value','replace 0 by 0.00001'], bbox_to_anchor=(1.05, 1), loc='upper left')
    #ax1.set(xlim=ax1.get_xlim(), ylim=ax1.get_ylim())
    diag_x = list(ax1.get_xlim())
    diag_y = list(ax1.get_ylim())
    ax1.plot(diag_x, diag_y, c='black')
    ax1.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax1.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    
    #print(len(x1),len(y1),len(x1_f),len(y1_f), blue_p_x, blue_p_y)
    ax1.set_title('FULL GENE\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x1)-len(x1_f)
                                                                                             ,title_map_rna1,blue_p_x-blue_both
                                                                                             ,title_map_rna2,blue_p_y
                                                                                             ,blue_both))
    
    # fold change
    #### Log ####
    ALPHA1 = min(find_alpha(dall1), find_alpha(dall2))#非零的所有值的最小值*0.0001
    x1 = [i+ALPHA1 for i in dall1 ] #control
    y1 = [i+ALPHA1 for i in dall2 ] #mutant
    plot_value1 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)] # all'
    #fold = pd.DataFrame()
    #fold['Gene name'] = all_group2['Gene name']
    #fold['fold'] = plot_value1
    #print(plot_value1)
    plot_max = max(plot_value1)
    plot_min = min(plot_value1)    
    
    #tmp2 = pd.DataFrame({'FULL GENE':pd.Series(plot_value1)})
    #sns.boxplot(data=tmp2,showfliers = False, width=0.3, color="#FFE66F", showmeans=True, ax=ax2)#.set(xlabel='region',ylabel='log2 ('+s2+' / '+s1+')')
    #add_stat_annotation(box1,data=tmp2,
    #                box_pairs=[('UTR5','CDS'),('UTR5','UTR3'),('CDS','UTR3')],
    #                test='Mann-Whitney', text_format='full', loc='inside', verbose=0,fontsize=16)
    # plot point
    #sns.stripplot(color='grey',data=tmp2,s=2,alpha=0.3)
    #ax2.set_title(title_map_rna1 +' v.s. '+title_map_rna2+ ' '+title_map_gene+
                 #'\n #ALL='+str(len(plot_value1)), fontsize=14)
    #ax2.set_title('FULL GENE \n α='+str(round(ALPHA,7)))
    #ax2.set_ylabel('log2 ('+title_map_rna2+'+α / '+title_map_rna1+'+α)',fontsize=20)
    #ax2.axhline(y=0, c='black' ,linestyle="--")
    
    # fold change without blue point
    dall1_no_b = []
    dall2_no_b = []
    #out = []
    for i,j in zip(dall1 ,dall2):
        if i != 0 and j != 0:
            dall1_no_b.append(i)
            dall2_no_b.append(j)
            #out.append(math.log(j/i, 2))
        #else:
           #out.append('NULL')
    #fold['fold_without0'] = out
    #ALPHA2 = min(find_alpha(dall1_no_b), find_alpha(dall2_no_b))#非零的所有值的最小值*0.0001
    
    x1 = [i for i in dall1_no_b ] #control

    y1 = [i for i in dall2_no_b ] #mutant

    plot_value2 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)] # all'
    #print(plot_value1)
    plot_max = max(max(plot_value1), max(plot_value2))
    plot_min = min(min(plot_value1), min(plot_value2)) 
    
    tmp2 = pd.DataFrame({'FULL GENE':pd.Series(plot_value1), 'FULL GENE WITHOUT 0':pd.Series(plot_value2)})
    sns.boxplot(data=tmp2,showfliers = False, width=0.3, color="#FFE66F", showmeans=True, ax=ax2)#.set(xlabel='region',ylabel='log2 ('+s2+' / '+s1+')')
    #add_stat_annotation(box1,data=tmp2,
    #                box_pairs=[('UTR5','CDS'),('UTR5','UTR3'),('CDS','UTR3')],
    #                test='Mann-Whitney', text_format='full', loc='inside', verbose=0,fontsize=16)
    # plot point
    #sns.stripplot(color='grey',data=tmp2,s=2,alpha=0.3)
    #ax2.set_title(title_map_rna1 +' v.s. '+title_map_rna2+ ' '+title_map_gene+
                 #'\n #ALL='+str(len(plot_value1)), fontsize=14)
    ax2.set_title('FULL GENE α='+str(round(ALPHA1,7)))
    ax2.set_ylabel('log2 ('+title_map_rna2+'+α / '+title_map_rna1+'+α)')
    ax2.axhline(y=0, c='black' ,linestyle="--")

    #### non-Log ####
    ALPHA1 = min(find_alpha(dall1), find_alpha(dall2))#非零的所有值的最小值*0.0001
    x1 = [i+ALPHA1 for i in dall1 ] #control
    y1 = [i+ALPHA1 for i in dall2 ] #mutant
    plot_value1 = [mut/con for mut,con in zip(y1,x1)] # all'  
    
    dall1_no_b = []
    dall2_no_b = []
    for i,j in zip(dall1 ,dall2):
        if i != 0 and j != 0:
            dall1_no_b.append(i)
            dall2_no_b.append(j)
    
    x1 = [i for i in dall1_no_b ] #control
    y1 = [i for i in dall2_no_b ] #mutant
    plot_value2 = [mut/con for mut,con in zip(y1,x1)] # all'
    
    plot_max = max(max(plot_value1), max(plot_value2))
    plot_min = min(min(plot_value1), min(plot_value2)) 
    
    tmp2 = pd.DataFrame({'FULL GENE':pd.Series(plot_value1), 'FULL GENE WITHOUT 0':pd.Series(plot_value2)})
    sns.boxplot(data=tmp2,showfliers = False, width=0.3, color="#FFE66F", showmeans=True, ax=ax3)
    ax3.set_title('FULL GENE α='+str(round(ALPHA1,7)))
    ax3.set_ylabel(title_map_rna2+'+α / '+title_map_rna1+'+α')
    ax3.axhline(y=0, c='black' ,linestyle="--")

    plt.tight_layout()
    plt.savefig(res_path)
    
'''
def plot_scatter_2data_g22_new(GROUP_NAME1, GROUP_NAME2, FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, two_set_group_dict, FILTER, rc):
    res_path = os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/fold_change_fig/2d_Sca_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(GROUP_NAME1)+'_'+str(GROUP_NAME2)+'_L'+str(FILTER)+'.png')
    #res_path = os.path.abspath(__file__+ '/../../../try/Sca_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    g1d1 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME1]
    g1d2 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME2]
    g2d1 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME1]
    g2d2 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME2]
    x1_f = []
    x2_f = []
    y1_f = []
    y2_f = []
    blue_p_x1 = 0
    blue_p_y1 = 0
    blue_p_x2 = 0
    blue_p_y2 = 0
    blue_both1 = 0
    blue_both2 = 0
    for i,j in zip(g1d1 ,g1d2):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y1 += 1
                else:
                    blue_both1 += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x1 += 1
            x1_f.append(math.log(i, 2))
            y1_f.append(math.log(j, 2))
            
    for i,j in zip(g2d1, g2d2):
        if i == 0 or j == 0:
            if i == 0:
                i = 0.00001
                if j != 0:
                    blue_p_y2 += 1
                else:
                    blue_both2 += 1
            if j == 0:
                j = 0.00001
                if i != 0:
                    blue_p_x2 += 1
            x2_f.append(math.log(i, 2))
            y2_f.append(math.log(j, 2))
    x1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in g1d1]
    x2 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in g2d1]
    
    y1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in g1d2]
    y2 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in g2d2]


    #f, ax = plt.subplots(figsize=(8, 8))
    fig , ((ax1,ax2)) = plt.subplots(1,2,figsize=(24,13))
    fig.suptitle(title_map_rna1 +' v.s. '+title_map_rna2+
                 ', '+title_map_gene[str(GROUP_NAME1)]+' v.s. '+title_map_gene[str(GROUP_NAME2)]+
                 '\n #{}='.format(title_map_gene[str(GROUP_NAME1)])+str(len(x1))+
                 '\n #{}='.format(title_map_gene[str(GROUP_NAME2)])+str(len(x2)),fontsize=18,y=1)
    if rc == 'need_rc':
        ylab = 'read count'
    else:
        ylab = 'site'
    fig.text(0.5, 0.05, 'log2('+ylab+') '+title_map_rna1, ha='center',size=20)
    fig.text(0.07, 0.5, 'log2('+ylab+') '+title_map_rna2, va='center', rotation='vertical',size=20)
    
    ax1.scatter(x1, y1, s=1, c='#AE0000',alpha=0.5)
    ax1.scatter(x1_f, y1_f, s=1, c='#0066FF')
    
    ax1.set(xlim=ax1.get_xlim(), ylim=ax1.get_ylim())
    diag_x = list(ax1.get_xlim())
    diag_y = list(ax1.get_ylim())
    ax1.plot(diag_x, diag_y, c='black')
    ax1.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax1.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    #print(len(x1),len(y1),len(x1_f),len(y1_f),blue_p_x1,blue_p_y1)
    ax1.axes.set_title(title_map_gene[str(GROUP_NAME1)]+'\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x1)-len(x1_f)
                                                                                             ,title_map_rna1,blue_p_x1-blue_both1
                                                                                             ,title_map_rna2,blue_p_y1
                                                                                              ,blue_both1))
    
    ax2.scatter(x2, y2, s=1, c='#AE0000',alpha=0.5)
    ax2.scatter(x2_f, y2_f, s=1, c='#0066FF')
    plt.legend(['valid value','replace 0 by 0.00001'], bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.set(xlim=ax2.get_xlim(), ylim=ax2.get_ylim())
    diag_x = list(ax2.get_xlim())
    diag_y = list(ax2.get_ylim())
    ax2.plot(diag_x, diag_y, c='black')
    ax2.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax2.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    
    #print(len(x2),len(y2),len(x2_f),len(y2_f),blue_p_x2,blue_p_y2)
    ax2.axes.set_title(title_map_gene[str(GROUP_NAME2)]+'\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x2)-len(x2_f)
                                                                                             ,title_map_rna1,blue_p_x2-blue_both2
                                                                                             ,title_map_rna2,blue_p_y2
                                                                                              ,blue_both2))
    
    plt.savefig(res_path)
    #plt.show()
def plot_ratio_pseuodo_count_g22_2data_new(GROUP_NAME1, GROUP_NAME2, FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, two_set_group_dict, FILTER, rc):
    
    g1d1 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME1]
    g1d2 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME2]
    g2d1 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME1]
    g2d2 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME2]
    
    #### fold change with blue point
    ALPHA = min(find_alpha(g1d1), find_alpha(g1d2), find_alpha(g2d1), find_alpha(g2d2))#非零的所有值的最小值*0.0001
                
    res_path = os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/fold_change_fig/2d_with0_FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(GROUP_NAME1)+'_'+str(GROUP_NAME2)+'_L'+str(FILTER)+'.png')
    #res_path = os.path.abspath(__file__+ '/../../../try/FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    x1 = [i+ALPHA for i in g1d1 ] #control1
    x2 = [i+ALPHA for i in g2d1 ] #control2

    y1 = [i+ALPHA for i in g1d2 ] #mutant1
    y2 = [i+ALPHA for i in g2d2 ] #mutant2

    plot_value1 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)]
    plot_value2 = [math.log(mut/con, 2) for mut,con in zip(y2,x2)]
    
    plot_max = max(max(plot_value1),max(plot_value2))
    plot_min = min(min(plot_value1),min(plot_value2))
    f, ax = plt.subplots(figsize=(12, 10))
    tmp2 = pd.DataFrame({title_map_gene[str(GROUP_NAME1)]:pd.Series(plot_value1), title_map_gene[str(GROUP_NAME2)]:pd.Series(plot_value2)})
    box1 = sns.boxplot(data=tmp2,showfliers = False, width=0.3, color="#FFE66F", showmeans=True)#.set(xlabel='region',ylabel='log2 ('+s2+' / '+s1+')')
    add_stat_annotation(box1,data=tmp2,
                    box_pairs=[(title_map_gene[str(GROUP_NAME1)],title_map_gene[str(GROUP_NAME2)])],
                    test='Mann-Whitney', text_format='full', loc='inside', verbose=0,fontsize=16)
    # plot point
    plt.xticks(fontsize = 20) 
    plt.yticks(fontsize = 20)
    #sns.stripplot(color='grey',data=tmp2,s=2,alpha=0.3)
    box1.axes.set_title(title_map_rna1 +' v.s. '+title_map_rna2+
                 '\n'+title_map_gene[str(GROUP_NAME1)]+' v.s. '+title_map_gene[str(GROUP_NAME2)]+
                 '\n #{}='.format(title_map_gene[str(GROUP_NAME1)])+str(len(plot_value1))+
                 '\n #{}='.format(title_map_gene[str(GROUP_NAME2)])+str(len(plot_value2))+
                        '\n # α={}'.format(round(ALPHA,7), fontsize=14))
    
    box1.set_xticklabels(box1.get_xmajorticklabels(), fontsize = 20)
    box1.set_ylabel('log2 ('+title_map_rna2+'+α / '+title_map_rna1+'+α)',fontsize=20)
    ax.axhline(y=0, c='black' ,linestyle="--")
    
    plt.tight_layout()
    plt.savefig(res_path) 
    
    #plt.show()
    
    #### fold change without blue point
    g1d1_no_b = []
    g1d2_no_b = []
    g2d1_no_b = []
    g2d2_no_b = []
    
    for a,b in zip(g1d1 ,g1d2):
        if a != 0 and b != 0:
            g1d1_no_b.append(a)
            g1d2_no_b.append(b)
            
    for c,d in zip(g2d1,g2d2):
        if c != 0 and d != 0:
            g2d1_no_b.append(c)
            g2d2_no_b.append(d)         
    #ALPHA = min(find_alpha(g1d1_no_b), find_alpha(g1d2_no_b), find_alpha(g2d1_no_b), find_alpha(g2d2_no_b))#非零的所有值的最小值*0.0001
      
    res_path = os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/fold_change_fig/2d_without0_FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(GROUP_NAME1)+'_'+str(GROUP_NAME2)+'_L'+str(FILTER)+'.png')
    #res_path = os.path.abspath(__file__+ '/../../../try/FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    x1 = [i for i in g1d1_no_b ] #control1
    x2 = [i for i in g2d1_no_b ] #control2

    y1 = [i for i in g1d2_no_b ] #mutant1
    y2 = [i for i in g2d2_no_b ] #mutant2

    plot_value1 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)]
    plot_value2 = [math.log(mut/con, 2) for mut,con in zip(y2,x2)]
    
    plot_max = max(max(plot_value1),max(plot_value2))
    plot_min = min(min(plot_value1),min(plot_value2))
    f, ax = plt.subplots(figsize=(12, 10))
    tmp2 = pd.DataFrame({title_map_gene[str(GROUP_NAME1)]:pd.Series(plot_value1), title_map_gene[str(GROUP_NAME2)]:pd.Series(plot_value2)})
    box1 = sns.boxplot(data=tmp2,showfliers = False, width=0.3, color="#FFE66F", showmeans=True)#.set(xlabel='region',ylabel='log2 ('+s2+' / '+s1+')')
    add_stat_annotation(box1,data=tmp2,
                    box_pairs=[(title_map_gene[str(GROUP_NAME1)],title_map_gene[str(GROUP_NAME2)])],
                    test='Mann-Whitney', text_format='full', loc='inside', verbose=0,fontsize=16)
    # plot point
    plt.xticks(fontsize = 20) 
    plt.yticks(fontsize = 20)
    #sns.stripplot(color='grey',data=tmp2,s=2,alpha=0.3)
    box1.axes.set_title(title_map_rna1 +' v.s. '+title_map_rna2+
                 '\n'+title_map_gene[str(GROUP_NAME1)]+' v.s. '+title_map_gene[str(GROUP_NAME2)]+
                 '\n #{}='.format(title_map_gene[str(GROUP_NAME1)])+str(len(plot_value1))+
                 '\n #{}='.format(title_map_gene[str(GROUP_NAME2)])+str(len(plot_value2)), fontsize=14)
    
    box1.set_xticklabels(box1.get_xmajorticklabels(), fontsize = 20)
    box1.set_ylabel('log2 ('+title_map_rna2+'+α / '+title_map_rna1+'+α)',fontsize=20)
    ax.axhline(y=0, c='black' ,linestyle="--")
    
    plt.tight_layout()
    plt.savefig(res_path) 
    #plt.show()
    
    plt.clf()
    plt.close()
    del tmp2
    gc.collect() 
    
def plot_ratio_pseuodo_count_g22_2data_3region(GROUP_NAME1, GROUP_NAME2, FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, two_set_group_dict, FILTER, rc):
    g1d1_1 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME1][0]
    g1d2_1 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME2][0]
    g2d1_1 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME1][0]
    g2d2_1 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME2][0] # utr5
    
    g1d1_2 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME1][1]
    g1d2_2 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME2][1]
    g2d1_2 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME1][1]
    g2d2_2 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME2][1] # cds
    
    g1d1_3 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME1][2]
    g1d2_3 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME2][2]
    g2d1_3 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME1][2]
    g2d2_3 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME2][2] # utr3
    
    
    #### with 0
    
    ALPHA = min(find_alpha(g1d1_1), find_alpha(g1d2_1), find_alpha(g2d1_1), find_alpha(g2d2_1),
               find_alpha(g1d1_2), find_alpha(g1d2_2), find_alpha(g2d1_2), find_alpha(g2d2_2),
               find_alpha(g1d1_3), find_alpha(g1d2_3), find_alpha(g2d1_3), find_alpha(g2d2_3))#非零的所有值的最小值*0.0001
                
    res_path = os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/fold_change_fig/2d_with0_FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(GROUP_NAME1)+'_'+str(GROUP_NAME2)+'_L'+str(FILTER)+'.png')
    #res_path = os.path.abspath(__file__+ '/../../../try/2d_FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    x1 = [i+ALPHA for i in g1d1_1 ] #control1
    x2 = [i+ALPHA for i in g2d1_1 ] #control2
    y1 = [i+ALPHA for i in g1d2_1 ] #mutant1
    y2 = [i+ALPHA for i in g2d2_1 ] #mutant2
    
    x3 = [i+ALPHA for i in g1d1_2 ] #control1
    x4 = [i+ALPHA for i in g2d1_2 ] #control2
    y3 = [i+ALPHA for i in g1d2_2 ] #mutant1
    y4 = [i+ALPHA for i in g2d2_2 ] #mutant2
    
    x5 = [i+ALPHA for i in g1d1_3 ] #control1
    x6 = [i+ALPHA for i in g2d1_3 ] #control2
    y5 = [i+ALPHA for i in g1d2_3 ] #mutant1
    y6 = [i+ALPHA for i in g2d2_3 ] #mutant2

    plot_value1 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)]
    plot_value2 = [math.log(mut/con, 2) for mut,con in zip(y2,x2)]
    
    plot_value3 = [math.log(mut/con, 2) for mut,con in zip(y3,x3)]
    plot_value4 = [math.log(mut/con, 2) for mut,con in zip(y4,x4)]
    
    plot_value5 = [math.log(mut/con, 2) for mut,con in zip(y5,x5)]
    plot_value6 = [math.log(mut/con, 2) for mut,con in zip(y6,x6)]
    plot_max = max(max(plot_value1),max(plot_value2),max(plot_value3),max(plot_value4),max(plot_value5),max(plot_value6))
    plot_min = min(min(plot_value1),min(plot_value2),min(plot_value3),min(plot_value4),min(plot_value5),min(plot_value6))
    
    f, ax = plt.subplots(figsize=(20, 13))
    tmp2 = pd.DataFrame({title_map_gene[str(GROUP_NAME1)]+'_UTR5':pd.Series(plot_value1), title_map_gene[str(GROUP_NAME2)]+'_UTR5':pd.Series(plot_value2),
                        title_map_gene[str(GROUP_NAME1)]+'_CDS':pd.Series(plot_value3), title_map_gene[str(GROUP_NAME2)]+'_CDS':pd.Series(plot_value4),
                        title_map_gene[str(GROUP_NAME1)]+'_UTR3':pd.Series(plot_value5), title_map_gene[str(GROUP_NAME2)]+'_UTR3':pd.Series(plot_value6)})
    
    my_pal = {title_map_gene[str(GROUP_NAME1)]+'_UTR5':'#FFA07A',
              title_map_gene[str(GROUP_NAME2)]+'_UTR5':'#FF7F50',
              title_map_gene[str(GROUP_NAME1)]+'_CDS':'#98F898', 
              title_map_gene[str(GROUP_NAME2)]+'_CDS':'#90EE90',
              title_map_gene[str(GROUP_NAME1)]+'_UTR3':'#BA55D3',
              title_map_gene[str(GROUP_NAME2)]+'_UTR3':'#8A2BE2'}
    
    box1 = sns.boxplot(data=tmp2,showfliers = False, width=0.3, palette=my_pal, showmeans=True)#.set(xlabel='region',ylabel='log2 ('+s2+' / '+s1+')')
    add_stat_annotation(box1,data=tmp2,
                    box_pairs=[(title_map_gene[str(GROUP_NAME1)]+'_UTR5',title_map_gene[str(GROUP_NAME2)]+'_UTR5'),
                               (title_map_gene[str(GROUP_NAME1)]+'_CDS',title_map_gene[str(GROUP_NAME2)]+'_CDS'),
                               (title_map_gene[str(GROUP_NAME1)]+'_UTR3',title_map_gene[str(GROUP_NAME2)]+'_UTR3')],
                    test='Mann-Whitney', text_format='full', loc='inside', verbose=0,fontsize=16)
    # plot point
    plt.xticks(fontsize = 20) 
    plt.yticks(fontsize = 20)
    #sns.stripplot(color='grey',data=tmp2,s=2,alpha=0.3)
    box1.axes.set_title(title_map_rna1 +' v.s. '+title_map_rna2+
                 '\n'+title_map_gene[str(GROUP_NAME1)]+' v.s. '+title_map_gene[str(GROUP_NAME2)]+
                 '\n #{} UTR5={}, CDS={}, UTR3={}'.format(title_map_gene[str(GROUP_NAME1)],str(len(plot_value1)),str(len(plot_value3)),str(len(plot_value5)))+
                 '\n #{} UTR5={}, CDS={}, UTR3={}'.format(title_map_gene[str(GROUP_NAME2)],str(len(plot_value2)),str(len(plot_value4)),str(len(plot_value6)))+
                    '\n # α={}'.format(round(ALPHA,7)), fontsize=14)
    
    box1.set_xticklabels(box1.get_xmajorticklabels(),fontsize=13)
    box1.set_ylabel('log2 ('+title_map_rna2+'+α / '+title_map_rna1+'+α)',fontsize=20)
    ax.axhline(y=0, c='black' ,linestyle="--")
    
    plt.savefig(res_path)
    #plt.show()
    
    #### without 0
    
    
    #### fold change without blue point
    g1d1_1_no_b = []
    g1d2_1_no_b = []
    g2d1_1_no_b = []
    g2d2_1_no_b = []
    
    g1d1_2_no_b = []
    g1d2_2_no_b = []
    g2d1_2_no_b = []
    g2d2_2_no_b = []
    
    g1d1_3_no_b = []
    g1d2_3_no_b = []
    g2d1_3_no_b = []
    g2d2_3_no_b = []
    
    for a,b in zip(g1d1_1 ,g1d2_1):
        if a != 0 and b != 0:
            g1d1_1_no_b.append(a)
            g1d2_1_no_b.append(b)
            
    for c,d in zip(g2d1_1,g2d2_1):
        if c != 0 and d != 0:
            g2d1_1_no_b.append(c)
            g2d2_1_no_b.append(d)    
            
    for a,b in zip(g1d1_2 ,g1d2_2):
        if a != 0 and b != 0:
            g1d1_2_no_b.append(a)
            g1d2_2_no_b.append(b)
            
    for c,d in zip(g2d1_2,g2d2_2):
        if c != 0 and d != 0:
            g2d1_2_no_b.append(c)
            g2d2_2_no_b.append(d)  
            
    for a,b in zip(g1d1_3 ,g1d2_3):
        if a != 0 and b != 0:
            g1d1_3_no_b.append(a)
            g1d2_3_no_b.append(b)
            
    for c,d in zip(g2d1_3,g2d2_3):
        if c != 0 and d != 0:
            g2d1_3_no_b.append(c)
            g2d2_3_no_b.append(d)  
            
    ALPHA = min(find_alpha(g1d1_1_no_b), find_alpha(g1d2_1_no_b), find_alpha(g2d1_1_no_b), find_alpha(g2d2_1_no_b),
               find_alpha(g1d1_2_no_b), find_alpha(g1d2_2_no_b), find_alpha(g2d1_2_no_b), find_alpha(g2d2_2_no_b),
               find_alpha(g1d1_3_no_b), find_alpha(g1d2_3_no_b), find_alpha(g2d1_3_no_b), find_alpha(g2d2_3_no_b))#非零的所有值的最小值*0.0001

    res_path = os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/fold_change_fig/2d_without0_FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(GROUP_NAME1)+'_'+str(GROUP_NAME2)+'_L'+str(FILTER)+'.png')
    #res_path = os.path.abspath(__file__+ '/../../../try/2d_FOLD_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    x1 = [i for i in g1d1_1_no_b ] #control1
    x2 = [i for i in g2d1_1_no_b ] #control2
    y1 = [i for i in g1d2_1_no_b ] #mutant1
    y2 = [i for i in g2d2_1_no_b ] #mutant2
    
    x3 = [i for i in g1d1_2_no_b ] #control1
    x4 = [i for i in g2d1_2_no_b ] #control2
    y3 = [i for i in g1d2_2_no_b ] #mutant1
    y4 = [i for i in g2d2_2_no_b ] #mutant2
    
    x5 = [i for i in g1d1_3_no_b ] #control1
    x6 = [i for i in g2d1_3_no_b ] #control2
    y5 = [i for i in g1d2_3_no_b ] #mutant1
    y6 = [i for i in g2d2_3_no_b ] #mutant2

    plot_value1 = [math.log(mut/con, 2) for mut,con in zip(y1,x1)]
    plot_value2 = [math.log(mut/con, 2) for mut,con in zip(y2,x2)]
    
    plot_value3 = [math.log(mut/con, 2) for mut,con in zip(y3,x3)]
    plot_value4 = [math.log(mut/con, 2) for mut,con in zip(y4,x4)]
    
    plot_value5 = [math.log(mut/con, 2) for mut,con in zip(y5,x5)]
    plot_value6 = [math.log(mut/con, 2) for mut,con in zip(y6,x6)]
    plot_max = max(max(plot_value1),max(plot_value2),max(plot_value3),max(plot_value4),max(plot_value5),max(plot_value6))
    plot_min = min(min(plot_value1),min(plot_value2),min(plot_value3),min(plot_value4),min(plot_value5),min(plot_value6))
    
    f, ax = plt.subplots(figsize=(20, 13))
    tmp2 = pd.DataFrame({title_map_gene[str(GROUP_NAME1)]+'_UTR5':pd.Series(plot_value1), title_map_gene[str(GROUP_NAME2)]+'_UTR5':pd.Series(plot_value2),
                        title_map_gene[str(GROUP_NAME1)]+'_CDS':pd.Series(plot_value3), title_map_gene[str(GROUP_NAME2)]+'_CDS':pd.Series(plot_value4),
                        title_map_gene[str(GROUP_NAME1)]+'_UTR3':pd.Series(plot_value5), title_map_gene[str(GROUP_NAME2)]+'_UTR3':pd.Series(plot_value6)})
    
    my_pal = {title_map_gene[str(GROUP_NAME1)]+'_UTR5':'#FFA07A',
              title_map_gene[str(GROUP_NAME2)]+'_UTR5':'#FF7F50',
              title_map_gene[str(GROUP_NAME1)]+'_CDS':'#98F898', 
              title_map_gene[str(GROUP_NAME2)]+'_CDS':'#90EE90',
              title_map_gene[str(GROUP_NAME1)]+'_UTR3':'#BA55D3',
              title_map_gene[str(GROUP_NAME2)]+'_UTR3':'#8A2BE2'}
    
    box1 = sns.boxplot(data=tmp2,showfliers = False, width=0.3, palette=my_pal, showmeans=True)#.set(xlabel='region',ylabel='log2 ('+s2+' / '+s1+')')
    add_stat_annotation(box1,data=tmp2,
                    box_pairs=[(title_map_gene[str(GROUP_NAME1)]+'_UTR5',title_map_gene[str(GROUP_NAME2)]+'_UTR5'),
                               (title_map_gene[str(GROUP_NAME1)]+'_CDS',title_map_gene[str(GROUP_NAME2)]+'_CDS'),
                               (title_map_gene[str(GROUP_NAME1)]+'_UTR3',title_map_gene[str(GROUP_NAME2)]+'_UTR3')],
                    test='Mann-Whitney', text_format='full', loc='inside', verbose=0,fontsize=16)
    # plot point
    plt.xticks(fontsize = 20) 
    plt.yticks(fontsize = 20)
    #sns.stripplot(color='grey',data=tmp2,s=2,alpha=0.3)
    box1.axes.set_title(title_map_rna1 +' v.s. '+title_map_rna2+
                 '\n'+title_map_gene[str(GROUP_NAME1)]+' v.s. '+title_map_gene[str(GROUP_NAME2)]+
                 '\n #{} UTR5={}, CDS={}, UTR3={}'.format(title_map_gene[str(GROUP_NAME1)],str(len(plot_value1)),str(len(plot_value3)),str(len(plot_value5)))+
                 '\n #{} UTR5={}, CDS={}, UTR3={}'.format(title_map_gene[str(GROUP_NAME2)],str(len(plot_value2)),str(len(plot_value4)),str(len(plot_value6)))
                   , fontsize=14)
    
    box1.set_xticklabels(box1.get_xmajorticklabels(),fontsize=13)
    box1.set_ylabel('log2 ('+title_map_rna2+'+α / '+title_map_rna1+'+α)',fontsize=20)
    ax.axhline(y=0, c='black' ,linestyle="--")
    
    plt.savefig(res_path)
    #plt.show()
    
def plot_scatter_2data_g22_3region(GROUP_NAME1, GROUP_NAME2, FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, two_set_group_dict, FILTER, rc):
    res_path = os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/fold_change_fig/2d_Sca_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(GROUP_NAME1)+'_'+str(GROUP_NAME2)+'_L'+str(FILTER)+'.png')
    #res_path = os.path.abspath(__file__+ '/../../../try/3_2d_Sca_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(FILTER)+'.png')
    g1d1_1 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME1][0]
    g1d2_1 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME2][0]
    g2d1_1 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME1][0]
    g2d2_1 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME2][0] # utr5
    
    g1d1_2 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME1][1]
    g1d2_2 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME2][1]
    g2d1_2 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME1][1]
    g2d2_2 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME2][1] # cds
    
    g1d1_3 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME1][2]
    g1d2_3 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME2][2]
    g2d1_3 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME1][2]
    g2d2_3 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME2][2] # utr3
    
    #f, ax = plt.subplots(figsize=(8, 8))
    fig , ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(15,25))
    
    x_list = []
    y_list = []
    xf_list = []
    yf_list = []
    blue_p_y_list = []
    blue_p_x_list = []
    blue_both_list = []
    for i in range(3):
        g1d1 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME1][i]
        g1d2 = two_set_group_dict[str(GROUP_NAME1)+'_'+FILE_NAME2][i]
        g2d1 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME1][i]
        g2d2 = two_set_group_dict[str(GROUP_NAME2)+'_'+FILE_NAME2][i]
    
        x1_f = []
        x2_f = []
        y1_f = []
        y2_f = []
        blue_p_x1 = 0
        blue_p_y1 = 0
        blue_p_x2 = 0
        blue_p_y2 = 0
        blue_both1 = 0
        blue_both2 = 0
        for i,j in zip(g1d1 ,g1d2):
            if i == 0 or j == 0:
                if i == 0:
                    i = 0.00001
                    if j != 0:
                        blue_p_y1 += 1
                    else:
                        blue_both1 += 1
                if j == 0:
                    j = 0.00001
                    if i != 0:
                        blue_p_x1 += 1
                x1_f.append(math.log(i, 2))
                y1_f.append(math.log(j, 2))

        for i,j in zip(g2d1, g2d2):
            if i == 0 or j == 0:
                if i == 0:
                    i = 0.00001
                    if j != 0:
                        blue_p_y2 += 1
                    else:
                        blue_both2 += 1
                if j == 0:
                    j = 0.00001
                    if i != 0:
                        blue_p_x2 += 1
                x2_f.append(math.log(i, 2))
                y2_f.append(math.log(j, 2))
        x1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in g1d1]
        x2 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in g2d1]

        y1 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in g1d2]
        y2 = [math.log(i,2) if i!= 0 else math.log(0.00001,2) for i in g2d2]
        
        x_list.append([x1,x2])
        y_list.append([y1,y2])
        xf_list.append([x1_f,x2_f])
        yf_list.append([y1_f,y2_f])
        blue_p_y_list.append([blue_p_y1,blue_p_y2])
        blue_p_x_list.append([blue_p_x1,blue_p_x2])
        blue_both_list.append([blue_both1,blue_both2])
    
    fig.suptitle(title_map_rna1 +' v.s. '+title_map_rna2+
                 '\n'+title_map_gene[str(GROUP_NAME1)]+' v.s. '+title_map_gene[str(GROUP_NAME2)]+
                 '\n{} #UTR5={} #CDS={} #UTR3={}'.format(title_map_gene[str(GROUP_NAME1)],str(len(x_list[0][0])),str(len(x_list[1][0])),str(len(x_list[2][0])))+
                 '\n{} #UTR5={} #CDS={} #UTR3={}'.format(title_map_gene[str(GROUP_NAME2)],str(len(x_list[0][1])),str(len(x_list[1][1])),str(len(x_list[2][1]))),fontsize=17,y=0.97)
    if rc == 'need_rc':
        ylab = 'read count'
    else:
        ylab = 'site'
    fig.text(0.5, 0.05, 'log2('+ylab+') '+title_map_rna1, ha='center',size=25)
    fig.text(0.07, 0.5, 'log2('+ylab+') '+title_map_rna2, va='center', rotation='vertical',size=25)
    
    #### utr5
    ax1.scatter(x_list[0][0], y_list[0][0], s=1, c='#AE0000',alpha=0.5)
    ax1.scatter(xf_list[0][0], yf_list[0][0], s=1, c='#0066FF')
    
    ax1.set(xlim=ax1.get_xlim(), ylim=ax1.get_ylim())
    diag_x = list(ax1.get_xlim())
    diag_y = list(ax1.get_ylim())
    ax1.plot(diag_x, diag_y, c='black')
    ax1.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax1.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    ax1.axes.set_title(title_map_gene[str(GROUP_NAME1)]+'\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x_list[0][0])-len(xf_list[0][0])
                                                                                             ,title_map_rna1,blue_p_x_list[0][0]-blue_both_list[0][0]
                                                                                             ,title_map_rna2,blue_p_y_list[0][0]
                                                                                              ,blue_both_list[0][0]))
    
    ax2.scatter(x_list[0][1], y_list[0][1], s=1, c='#AE0000',alpha=0.5)
    ax2.scatter(xf_list[0][1], yf_list[0][1], s=1, c='#0066FF')
    plt.legend(['valid value','replace 0 by 0.00001'], bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.set(xlim=ax2.get_xlim(), ylim=ax2.get_ylim())
    diag_x = list(ax2.get_xlim())
    diag_y = list(ax2.get_ylim())
    ax2.plot(diag_x, diag_y, c='black')
    ax2.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax2.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    ax2.axes.set_title(title_map_gene[str(GROUP_NAME2)]+'\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x_list[0][1])-len(xf_list[0][1])
                                                                                             ,title_map_rna1,blue_p_x_list[0][1]-blue_both_list[0][1]
                                                                                             ,title_map_rna2,blue_p_y_list[0][1]
                                                                                              ,blue_both_list[0][1]))
    
    #### utr3
    ax5.scatter(x_list[2][0], y_list[2][0], s=1, c='#AE0000',alpha=0.5)
    ax5.scatter(xf_list[2][0], yf_list[2][0], s=1, c='#0066FF')
    
    ax5.set(xlim=ax5.get_xlim(), ylim=ax5.get_ylim())
    diag_x = list(ax5.get_xlim())
    diag_y = list(ax5.get_ylim())
    ax5.plot(diag_x, diag_y, c='black')
    ax5.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax5.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    ax5.axes.set_title(title_map_gene[str(GROUP_NAME1)]+'\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x_list[2][0])-len(xf_list[2][0])
                                                                                             ,title_map_rna1,blue_p_x_list[2][0]-blue_both_list[2][0]
                                                                                             ,title_map_rna2,blue_p_y_list[2][0]
                                                                                              ,blue_both_list[2][0]))
    
    ax6.scatter(x_list[2][1], y_list[2][1], s=1, c='#AE0000',alpha=0.5)
    ax6.scatter(xf_list[2][1], yf_list[2][1], s=1, c='#0066FF')
    ax6.set(xlim=ax6.get_xlim(), ylim=ax6.get_ylim())
    diag_x = list(ax6.get_xlim())
    diag_y = list(ax6.get_ylim())
    ax6.plot(diag_x, diag_y, c='black')
    ax6.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax6.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    ax6.axes.set_title(title_map_gene[str(GROUP_NAME2)]+'\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x_list[2][1])-len(xf_list[2][1])
                                                                                             ,title_map_rna1,blue_p_x_list[2][1]-blue_both_list[2][1]
                                                                                             ,title_map_rna2,blue_p_y_list[2][1]
                                                                                              ,blue_both_list[2][1]))
    #### cds
    ax3.scatter(x_list[1][0], y_list[1][0], s=1, c='#AE0000',alpha=0.5)
    ax3.scatter(xf_list[1][0], yf_list[1][0], s=1, c='#0066FF')
    
    ax3.set(xlim=ax3.get_xlim(), ylim=ax3.get_ylim())
    diag_x = list(ax3.get_xlim())
    diag_y = list(ax3.get_ylim())
    ax3.plot(diag_x, diag_y, c='black')
    ax3.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax3.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    ax3.axes.set_title(title_map_gene[str(GROUP_NAME1)]+'\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x_list[1][0])-len(xf_list[1][0])
                                                                                             ,title_map_rna1,blue_p_x_list[1][0]-blue_both_list[1][0]
                                                                                             ,title_map_rna2,blue_p_y_list[1][0]
                                                                                              ,blue_both_list[1][0]))
    
    ax4.scatter(x_list[1][1], y_list[1][1], s=1, c='#AE0000',alpha=0.5)
    ax4.scatter(xf_list[1][1], yf_list[1][1], s=1, c='#0066FF')
    ax4.set(xlim=ax4.get_xlim(), ylim=ax4.get_ylim())
    diag_x = list(ax4.get_xlim())
    diag_y = list(ax4.get_ylim())
    ax4.plot(diag_x, diag_y, c='black')
    ax4.plot([diag_x[0]+1, diag_x[1]+1], [diag_y[0], diag_y[1]], c='black' ,linestyle="--")
    ax4.plot([diag_x[0], diag_x[1]], [diag_y[0]+1, diag_y[1]+1], c='black' ,linestyle="--")
    ax4.axes.set_title(title_map_gene[str(GROUP_NAME2)]+'\nRED POINT={}\nBLUE POINT IN {}={}, IN {}={}, BOTH={}'.format(len(x_list[1][1])-len(xf_list[1][1])
                                                                                             ,title_map_rna1,blue_p_x_list[1][1]-blue_both_list[1][1]
                                                                                             ,title_map_rna2,blue_p_y_list[1][1]
                                                                                              ,blue_both_list[1][1]))
    
    plt.savefig(res_path)
    #plt.show()

'''