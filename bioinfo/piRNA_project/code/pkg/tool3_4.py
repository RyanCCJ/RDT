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

# main
def tool3_4(data):
    all_mRNA = pd.read_csv(os.path.abspath(__file__+ '/../../../../data/mRNA_WS275_IDtoName.csv'))
    gene_path = data['gene_path']
    title_map_gene = data['title_map_gene']
    title_map_rna = data['title_map_rna']

    for ANALYZE_TYPE in data['ANALYZE_TYPE']:
        print('===={}===='.format(ANALYZE_TYPE))
        if ANALYZE_TYPE == 'TARGET_SCORE':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME1, FILE_NAME2 in zip(data['COMPARE_DATA1'], data['COMPARE_DATA2']):
                    for EXTEND in data['EXTEND_RANGE']:
                        for gene in data['ANALYZE_GROUP']:
                            df = pd.read_csv(os.path.abspath(__file__+ '/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME1+'_all_mRNA_tool3_START_'+str(EXTEND)+'.csv'))
                            df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])

                            df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME1+'_all_mRNA_tool3_STOP_'+str(EXTEND)+'.csv'))
                            df_withid2 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid2['Gene ID'] = df_withid2['Gene ID'].apply(lambda x:x.split('=')[1])

                            df = pd.read_csv(os.path.abspath(__file__+ '/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME2+'_all_mRNA_tool3_START_'+str(EXTEND)+'.csv'))
                            df_withid3 = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid3['Gene ID'] = df_withid3['Gene ID'].apply(lambda x:x.split('=')[1])

                            df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME2+'_all_mRNA_tool3_STOP_'+str(EXTEND)+'.csv'))
                            df_withid4 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid4['Gene ID'] = df_withid4['Gene ID'].apply(lambda x:x.split('=')[1])

                            print(rc, FILE_NAME1, FILE_NAME2, EXTEND, gene)
                            gene_name = check_new(data, gene)
                            non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste = calculate_line(df_withid, rc, FILE_NAME1, EXTEND, 'START', gene_path, gene, gene_name)
                            non_zero2, mRNA_num2, avg2, med2, avg_plus_ste2, avg_minus_ste2, std2, ste2 = calculate_line(df_withid2, rc, FILE_NAME1, EXTEND, 'STOP', gene_path, gene, gene_name)
                            non_zero3, mRNA_num3, avg3, med3, avg_plus_ste3, avg_minus_ste3, std3, ste3 = calculate_line(df_withid3, rc, FILE_NAME2, EXTEND, 'START', gene_path, gene, gene_name)
                            non_zero4, mRNA_num4, avg4, med4, avg_plus_ste4, avg_minus_ste4, std4, ste4 = calculate_line(df_withid4, rc, FILE_NAME2, EXTEND, 'STOP', gene_path, gene, gene_name)

                            codon_plot(rc, FILE_NAME1, title_map_rna[FILE_NAME1], FILE_NAME2, title_map_rna[FILE_NAME2], EXTEND, gene, title_map_gene[gene], non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste,
                            non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2,
                            non_zero3, mRNA_num3, avg3, avg_plus_ste3, avg_minus_ste3,
                            non_zero4, mRNA_num4, avg4, avg_plus_ste4, avg_minus_ste4)
                                    
            
        elif ANALYZE_TYPE == 'G22':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME1, FILE_NAME2 in zip(data['COMPARE_DATA1'], data['COMPARE_DATA2']):
                    for EXTEND in data['EXTEND_RANGE']:
                        for gene in data['ANALYZE_GROUP']:
                            df = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME1+'_all_mRNA_tool3_START_'+str(EXTEND)+'.csv'))
                            df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])

                            df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME1+'_all_mRNA_tool3_STOP_'+str(EXTEND)+'.csv'))
                            df_withid2 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid2['Gene ID'] = df_withid2['Gene ID'].apply(lambda x:x.split('=')[1])

                            df = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME2+'_all_mRNA_tool3_START_'+str(EXTEND)+'.csv'))
                            df_withid3 = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid3['Gene ID'] = df_withid3['Gene ID'].apply(lambda x:x.split('=')[1])

                            df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME2+'_all_mRNA_tool3_STOP_'+str(EXTEND)+'.csv'))
                            df_withid4 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid4['Gene ID'] = df_withid4['Gene ID'].apply(lambda x:x.split('=')[1])

                            print(rc, FILE_NAME1, FILE_NAME2, EXTEND, gene)
                            gene_name = check_new(data, gene)
                            non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste = calculate_line_g22(df_withid, rc, FILE_NAME1, EXTEND, 'START', gene_path, gene, gene_name)
                            non_zero2, mRNA_num2, avg2, med2, avg_plus_ste2, avg_minus_ste2, std2, ste2 = calculate_line_g22(df_withid2, rc, FILE_NAME1, EXTEND, 'STOP', gene_path, gene, gene_name)
                            non_zero3, mRNA_num3, avg3, med3, avg_plus_ste3, avg_minus_ste3, std3, ste3 = calculate_line_g22(df_withid3, rc, FILE_NAME2, EXTEND, 'START', gene_path, gene, gene_name)
                            non_zero4, mRNA_num4, avg4, med4, avg_plus_ste4, avg_minus_ste4, std4, ste4 = calculate_line_g22(df_withid4, rc, FILE_NAME2, EXTEND, 'STOP', gene_path, gene, gene_name)

                            codon_plot_g22(rc, FILE_NAME1, title_map_rna[FILE_NAME1], FILE_NAME2, title_map_rna[FILE_NAME2], EXTEND, gene, title_map_gene[gene], non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste,
                            non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2,
                            non_zero3, mRNA_num3, avg3, avg_plus_ste3, avg_minus_ste3,
                            non_zero4, mRNA_num4, avg4, avg_plus_ste4, avg_minus_ste4)
                                    
                            #### HEAD TAIL
                            #df = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME1+'_all_mRNA_tool3_HEAD_'+str(EXTEND)+'.csv'))
                            #df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])

                            #df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME1+'_all_mRNA_tool3_TAIL_'+str(EXTEND)+'.csv'))
                            #df_withid2 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid2['Gene ID'] = df_withid2['Gene ID'].apply(lambda x:x.split('=')[1])

                            #df = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME2+'_all_mRNA_tool3_HEAD_'+str(EXTEND)+'.csv'))
                            #df_withid3 = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid3['Gene ID'] = df_withid3['Gene ID'].apply(lambda x:x.split('=')[1])

                            #df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/'+FILE_NAME2+'_all_mRNA_tool3_TAIL_'+str(EXTEND)+'.csv'))
                            #df_withid4 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid4['Gene ID'] = df_withid4['Gene ID'].apply(lambda x:x.split('=')[1])
                                    
                            #non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste = calculate_line_g22_HT(df_withid, rc, FILE_NAME1, EXTEND, 'HEAD', gene)
                            #non_zero2, mRNA_num2, avg2, med2, avg_plus_ste2, avg_minus_ste2, std2, ste2 = calculate_line_g22_HT(df_withid2, rc, FILE_NAME1, EXTEND, 'TAIL', gene)
                            #non_zero3, mRNA_num3, avg3, med3, avg_plus_ste3, avg_minus_ste3, std3, ste3 = calculate_line_g22_HT(df_withid3, rc, FILE_NAME2, EXTEND, 'HEAD', gene)
                            #non_zero4, mRNA_num4, avg4, med4, avg_plus_ste4, avg_minus_ste4, std4, ste4 = calculate_line_g22_HT(df_withid4, rc, FILE_NAME2, EXTEND, 'TAIL', gene)

                            #codon_plot_g22_HT(rc, FILE_NAME1, FILE_NAME2, EXTEND, gene, non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste,
                            #non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2,
                            #non_zero3, mRNA_num3, avg3, avg_plus_ste3, avg_minus_ste3,
                            #non_zero4, mRNA_num4, avg4, avg_plus_ste4, avg_minus_ste4)
                
# tool3 requirement
def calculate_line(df_withid, rc, FILE_NAME, EXTEND, SIDE, gene_path, gene, title_map_gene):

    #df_filter = mRNA_list(df_withid, gene)
    df_filter = add_two_mRNA_list(df_withid, gene_path, gene, title_map_gene)
    # calculate mRNA number
    show_range = list(np.arange(-EXTEND,EXTEND+1))
    mRNA_num = []
    avg = []
    med = []
    ste = []
    std = []
    non_zero = []
    for ind in show_range:
        a = list(df_filter[df_filter[str(ind)]!='N'][str(ind)])
        a = [float(i) for i in a]
        mRNA_num.append(len(a))
        avg.append(np.mean(a))
        med.append(np.median(a))
        std.append(np.std(a))
        ste.append(stats.sem(a))
        non_zero.append(len([nz for nz in a if nz != 0]))

    avg_plus_ste = [ i+j for i,j in zip(avg,ste)]
    avg_minus_ste = [ i-j for i,j in zip(avg,ste)]
    return non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste

def codon_plot(rc, FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, EXTEND, gene, title_map_gene, non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste,
               non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2,
               non_zero3, mRNA_num3, avg3, avg_plus_ste3, avg_minus_ste3,
               non_zero4, mRNA_num4, avg4, avg_plus_ste4, avg_minus_ste4):
    
    #res_path = os.path.abspath(__file__+ '/../../../usetRNA_target_score_'+rc+'_output/codon_fig/COMPARE_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(EXTEND)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/codon_fig/compare/codon_COMPARE_'+FILE_NAME1+'_'+FILE_NAME2+'_G'+str(gene)+'_L'+str(EXTEND)+rc+'.png')
    fig ,((ax01,ax03),(ax1,ax3),(ax2,ax4),(ax22,ax42)) = plt.subplots(4,2, sharex=True, figsize=(22,14),dpi=200)
    # start
    SIDE = 'START'
    #### old figure  
    ax01.set_title(title_map_gene+'_'+SIDE+str(EXTEND)+' codon')
    ax01.set_ylabel('number of mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax01.axvline(x=0,c='k',linestyle='dashed')
    ax01.plot(range_show,mRNA_num,color='#AA7700')
    ax01.plot(range_show,non_zero,color='#CC6600')
    ax01.plot(range_show,mRNA_num3,color='#668800')
    ax01.plot(range_show,non_zero3,color='green')
    ax01.tick_params(axis='y', labelsize=20)
    #ax01.legend([SIDE+' codon',title_map_rna1+' mRNA have position',title_map_rna1+' mRNA have site',
               #title_map_rna2+' mRNA have position',title_map_rna2+' mRNA have site'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    #### ratio figure
    ratio_data = [non_zero[i]/mRNA_num[i] for i in range(len(mRNA_num))]
    ratio_data3 = [non_zero3[i]/mRNA_num3[i] for i in range(len(mRNA_num3))]
    
    #plot mRNA info
    ax1.set_ylabel('ratio mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax1.axvline(x=0,c='k',linestyle='dashed')
    
    ax1.plot(range_show,ratio_data,color='#CC6600')
    ax1.plot(range_show,ratio_data3,color='green')
    ax1.tick_params(axis='y', labelsize=20)
    #ax1.legend([SIDE+' codon',title_map_rna1+' mRNA site ratio\n(site/position)',title_map_rna2+' mRNA site ratio\n(site/position)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    #ax1.legend([SIDE+' codon',title_map_rna1+' mRNA have position',title_map_rna1+' mRNA have site',
               #title_map_rna2+' mRNA have position',title_map_rna2+' mRNA have site'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

          
    #### codon figure
     
    MAX_Yc = max(max(avg_plus_ste),max(avg_plus_ste2), max(avg_plus_ste3), max(avg_plus_ste4))
    if rc == 'no_need_rc':
        ax2.set_ylabel('number of sites',fontsize=20)
    else:
        ax2.set_ylabel('read counts', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax2.axvline(x=0,c='k',linestyle='dashed')
    ax2.plot(range_show,avg,color='#CC6600')
    ax2.plot(range_show,avg3,color='green')
                
    ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
    ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
    ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)    
    
    ax2.plot(range_show,avg_plus_ste3,color='#227700',linestyle='dashed')
    ax2.plot(range_show,avg_minus_ste3,color='#227700',linestyle='dashed')
    ax2.fill_between(range_show, avg_plus_ste3, avg_minus_ste3, color='#227700',alpha=0.3)
    
    #ax2.legend([SIDE+' codon',title_map_rna1+' AVG(+/- 1STE)',title_map_rna2+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax2.tick_params(axis='y', labelsize=20)
    ax2.tick_params(axis='x', labelsize=20)
    #ax2.set_xlabel('nt', fontsize=20)
    ax2.set_ylim(0, MAX_Yc+0.1*MAX_Yc)
    
    #### codon distribution figure
    avg_sum = sum(avg)
    avg3_sum = sum(avg3)
    avg_dist = [i/avg_sum for i in avg]
    avg3_dist = [i/avg3_sum for i in avg3]

    avg_plus_ste_dist = [i/avg_sum for i in avg_plus_ste]
    avg_minus_ste_dist = [i/avg_sum for i in avg_minus_ste]

    avg_plus_ste3_dist = [i/avg3_sum for i in avg_plus_ste3]
    avg_minus_ste3_dist = [i/avg3_sum for i in avg_minus_ste3] 
    
    avg2_sum = sum(avg2)
    avg4_sum = sum(avg4)
    avg2_dist = [i/avg2_sum for i in avg2]
    avg4_dist = [i/avg4_sum for i in avg4]

    avg_plus_ste2_dist = [i/avg2_sum for i in avg_plus_ste2]
    avg_minus_ste2_dist = [i/avg2_sum for i in avg_minus_ste2]

    avg_plus_ste4_dist = [i/avg4_sum for i in avg_plus_ste4]
    avg_minus_ste4_dist = [i/avg4_sum for i in avg_minus_ste4] 
    
    MAX_Y = max(max(avg_plus_ste_dist),max(avg_plus_ste2_dist), max(avg_plus_ste3_dist), max(avg_plus_ste4_dist))
    
    if rc == 'no_need_rc':
        ylabel = 'number of sites'
    else:
        ylabel = 'read counts'
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax22.set_ylabel(ylabel+' distribution', fontsize=20)
    ax22.axvline(x=0,c='k',linestyle='dashed')
    ax22.plot(range_show,avg_dist,color='#CC6600')
    ax22.plot(range_show,avg3_dist,color='green')
                
    ax22.plot(range_show,avg_plus_ste_dist,color='orange',linestyle='dashed')
    ax22.plot(range_show,avg_minus_ste_dist,color='orange',linestyle='dashed')
    ax22.fill_between(range_show, avg_plus_ste_dist, avg_minus_ste_dist, color='orange',alpha=0.3)    
    
    ax22.plot(range_show,avg_plus_ste3_dist,color='#227700',linestyle='dashed')
    ax22.plot(range_show,avg_minus_ste3_dist,color='#227700',linestyle='dashed')
    ax22.fill_between(range_show, avg_plus_ste3_dist, avg_minus_ste3_dist, color='#227700',alpha=0.3)
    
    #ax22.legend([SIDE+' codon',title_map_rna1+' AVG(+/- 1STE)',title_map_rna2+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax22.tick_params(axis='y', labelsize=20)
    ax22.tick_params(axis='x', labelsize=20)
    ax22.set_xlabel('nt', fontsize=20)
    ax22.set_ylim(0, MAX_Y+0.1*MAX_Y)
    
    # STOP
    SIDE = 'STOP'
    #### old figure  
    ax03.set_title(title_map_gene+'_'+SIDE+str(EXTEND)+' codon')
    #ax03.set_ylabel('number of mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax03.axvline(x=0,c='k',linestyle='dashed')
    ax03.plot(range_show,mRNA_num2,color='#AA7700')
    ax03.plot(range_show,non_zero2,color='#CC6600')
    ax03.plot(range_show,mRNA_num4,color='#668800')
    ax03.plot(range_show,non_zero4,color='green')
    ax03.tick_params(axis='y', labelsize=20)
    ax03.legend([SIDE+' codon',title_map_rna1+' mRNA have position',title_map_rna1+' mRNA have site',
               title_map_rna2+' mRNA have position',title_map_rna2+' mRNA have site'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    #### ratio figure       
    ratio_data2 = [non_zero2[i]/mRNA_num2[i] for i in range(len(mRNA_num2))]
    ratio_data4 = [non_zero4[i]/mRNA_num4[i] for i in range(len(mRNA_num4))]
    #plot mRNA info
    #plt.figure(figsize=(14,8),dpi=150)
    #ax3.set_title(title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    
    #ax3.set_ylabel('ratio mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax3.axvline(x=0,c='k',linestyle='dashed')
    ax3.plot(range_show,ratio_data2,color='#CC6600')
    ax3.plot(range_show,ratio_data4,color='green')
    ax3.tick_params(axis='y', labelsize=20)
    ax3.legend([SIDE+' codon',title_map_rna1+' mRNA site ratio\n(site/position)',title_map_rna2+' mRNA site ratio\n(site/position)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    #### codon figure
    #plot site
    '''if rc == 'no_need_rc':
        ax4.set_ylabel('number of sites',fontsize=20)
    else:
        ax4.set_ylabel('read counts', fontsize=20)
    '''
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax4.axvline(x=0,c='k',linestyle='dashed')
    ax4.plot(range_show,avg2,color='#CC6600')
    ax4.plot(range_show,avg4,color='green')
    ax4.plot(range_show,avg_plus_ste2,color='orange',linestyle='dashed')
    ax4.plot(range_show,avg_minus_ste2,color='orange',linestyle='dashed')
    ax4.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='orange',alpha=0.3)
    
    ax4.plot(range_show,avg_plus_ste4,color='#227700',linestyle='dashed')
    ax4.plot(range_show,avg_minus_ste4,color='#227700',linestyle='dashed')
    ax4.fill_between(range_show, avg_plus_ste4, avg_minus_ste4, color='#227700',alpha=0.3)
    ax4.legend([SIDE+' codon',title_map_rna1+' AVG(+/- 1STE)',title_map_rna2+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax4.tick_params(axis='y', labelsize=20)
    ax4.set_ylim(0, MAX_Yc+0.1*MAX_Yc)
    
    #### codon distribution figure
   
    if rc == 'no_need_rc':
        ylabel = 'number of sites'
    else:
        ylabel = 'read counts'
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    #ax42.set_ylabel('number of sites',fontsize=20)
    ax42.axvline(x=0,c='k',linestyle='dashed')
    ax42.plot(range_show,avg2_dist,color='#CC6600')
    ax42.plot(range_show,avg4_dist,color='green')
    ax42.plot(range_show,avg_plus_ste2_dist,color='orange',linestyle='dashed')
    ax42.plot(range_show,avg_minus_ste2_dist,color='orange',linestyle='dashed')
    ax42.fill_between(range_show, avg_plus_ste2_dist, avg_minus_ste2_dist, color='orange',alpha=0.3)
    
    ax42.plot(range_show,avg_plus_ste4_dist,color='#227700',linestyle='dashed')
    ax42.plot(range_show,avg_minus_ste4_dist,color='#227700',linestyle='dashed')
    ax42.fill_between(range_show, avg_plus_ste4_dist, avg_minus_ste4_dist, color='#227700',alpha=0.3)
    ax42.legend([SIDE+' codon',title_map_rna1+' AVG(+/- 1STE)',title_map_rna2+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax42.tick_params(axis='y', labelsize=20)
    ax42.set_ylim(0, MAX_Y+0.1*MAX_Y)
    
    plt.xlabel('nt',fontsize=20)
    plt.xticks(fontsize=20)    
    plt.tight_layout()
    plt.savefig(res_path)
    tmp2 = pd.DataFrame({'number of mRNA position(start)':mRNA_num,
                         'number of mRNA2 position(start)':mRNA_num3,
                         'number of mRNA site(start)':non_zero,
                         'number of mRNA2 site(start)':non_zero3,
                         'avg(start)':avg,'avg2(start)':avg3,
                         'avg_plus_ste(start)':avg_plus_ste,'avg_plus_ste2(start)':avg_plus_ste3,
                         'avg_minus_ste(start)':avg_minus_ste,'avg_minus_ste2(start)':avg_minus_ste3,
                         'avg_dist_dist(start)':avg_dist,'avg2_dist(start)':avg3_dist,
                         'avg_plus_ste_dist(start)':avg_plus_ste_dist,'avg_plus_ste2_dist(start)':avg_plus_ste3_dist,
                         'avg_minus_ste_dist(start)':avg_minus_ste_dist,'avg_minus_ste2_dist(start)':avg_minus_ste3_dist,
                         
                         'number of mRNA position(stop)':mRNA_num2,
                         'number of mRNA2 position(stop)':mRNA_num4,
                         'number of mRNA site(stop)':non_zero2,
                         'number of mRNA2 site(stastoprt)':non_zero4,
                         'avg(stop)':avg2,'avg2(stop)':avg4,
                         'avg_plus_ste(stop)':avg_plus_ste2,'avg_plus_ste2(stop)':avg_plus_ste4,
                         'avg_minus_ste(stop)':avg_minus_ste2,'avg_minus_ste2(stop)':avg_minus_ste4,
                         'avg_dist(stop)':avg2_dist,'avg2_dist(stop)':avg4_dist,
                         'avg_plus_ste_dist(stop)':avg_plus_ste2_dist,'avg_plus_ste2_dist(stop)':avg_plus_ste4_dist,
                         'avg_minus_ste_dist(stop)':avg_minus_ste2_dist,'avg_minus_ste2_dist(stop)':avg_minus_ste4_dist,
    })
    #tmp2.to_csv(res_path.replace('png', 'csv'),index=False)
    plt.clf()
    plt.close()
    del non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste
    gc.collect()

def calculate_line_g22(df_withid, rc, FILE_NAME, EXTEND, SIDE, gene_path, gene, title_map_gene):

    df_filter = add_two_mRNA_list(df_withid, gene_path, gene, title_map_gene)

    # calculate mRNA number
    show_range = list(np.arange(-EXTEND,EXTEND+1))
    mRNA_num = []
    avg = []
    med = []
    ste = []
    std = []
    non_zero = []
    for ind in show_range:
        a = list(df_filter[df_filter[str(ind)]!='N'][str(ind)])
        a = [float(i) for i in a]
        mRNA_num.append(len(a))
        avg.append(np.mean(a))
        med.append(np.median(a))
        std.append(np.std(a))
        ste.append(stats.sem(a))
        non_zero.append(len([nz for nz in a if nz != 0]))

    avg_plus_ste = [ i+j for i,j in zip(avg,ste)]
    avg_minus_ste = [ i-j for i,j in zip(avg,ste)]
    return non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste

def calculate_line_g22_HT(df_withid, rc, FILE_NAME, EXTEND, SIDE, gene_path, gene, title_map_gene):

    df_filter = add_two_mRNA_list(df_withid, gene_path, gene, title_map_gene)

    # calculate mRNA number
    show_range = list(np.arange(EXTEND))
    mRNA_num = []
    avg = []
    med = []
    ste = []
    std = []
    non_zero = []
    for ind in show_range:
        a = list(df_filter[df_filter[str(ind)]!='N'][str(ind)])
        a = [float(i) for i in a]
        mRNA_num.append(len(a))
        avg.append(np.mean(a))
        med.append(np.median(a))
        std.append(np.std(a))
        ste.append(stats.sem(a))
        non_zero.append(len([nz for nz in a if nz != 0]))

    avg_plus_ste = [ i+j for i,j in zip(avg,ste)]
    avg_minus_ste = [ i-j for i,j in zip(avg,ste)]
    return non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste

def codon_plot_g22(rc, FILE_NAME1, title_map_rna1, FILE_NAME2, title_map_rna2, EXTEND, gene, title_map_gene, non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste,
               non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2,
               non_zero3, mRNA_num3, avg3, avg_plus_ste3, avg_minus_ste3,
               non_zero4, mRNA_num4, avg4, avg_plus_ste4, avg_minus_ste4):
    
    #res_path = os.path.abspath(__file__+ '/../../../usemiRNA_22G_output/codon_fig/COMPARE_'+FILE_NAME1+'_'+FILE_NAME2+'G'+str(gene)+'_L'+str(EXTEND)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/codon_fig/compare/codon_COMPARE_'+FILE_NAME1+'_'+FILE_NAME2+'_G'+str(gene)+'_L'+str(EXTEND)+'.png')
    fig ,((ax01,ax03),(ax1,ax3),(ax2,ax4),(ax22,ax42)) = plt.subplots(4,2, sharex=True, figsize=(22,14),dpi=200)
    # start
    SIDE = 'START'
    #### old figure  
    ax01.set_title(title_map_gene+'_'+SIDE+str(EXTEND)+' codon')
    ax01.set_ylabel('number of mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax01.axvline(x=0,c='k',linestyle='dashed')
    ax01.plot(range_show,mRNA_num,color='#AA7700')
    ax01.plot(range_show,non_zero,color='#CC6600')
    ax01.plot(range_show,mRNA_num3,color='#668800')
    ax01.plot(range_show,non_zero3,color='green')
    ax01.tick_params(axis='y', labelsize=20)
    #ax01.legend([SIDE+' codon',title_map_rna1+' mRNA have position',title_map_rna1+' mRNA have site',
               #title_map_rna2+' mRNA have position',title_map_rna2+' mRNA have site'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    #### ratio figure
    ratio_data = [non_zero[i]/mRNA_num[i] for i in range(len(mRNA_num))]
    ratio_data3 = [non_zero3[i]/mRNA_num3[i] for i in range(len(mRNA_num3))]
    
    #plot mRNA info
    ax1.set_ylabel('ratio mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax1.axvline(x=0,c='k',linestyle='dashed')
    
    ax1.plot(range_show,ratio_data,color='#CC6600')
    ax1.plot(range_show,ratio_data3,color='green')
    ax1.tick_params(axis='y', labelsize=20)
    #ax1.legend([SIDE+' codon',title_map_rna1+' mRNA site ratio\n(site/position)',title_map_rna2+' mRNA site ratio\n(site/position)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    #ax1.legend([SIDE+' codon',title_map_rna1+' mRNA have position',title_map_rna1+' mRNA have site',
               #title_map_rna2+' mRNA have position',title_map_rna2+' mRNA have site'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

          
    #### codon figure
     
    MAX_Yc = max(max(avg_plus_ste),max(avg_plus_ste2), max(avg_plus_ste3), max(avg_plus_ste4))
    if rc == 'no_need_rc':
        ax2.set_ylabel('number of sites',fontsize=20)
    else:
        ax2.set_ylabel('read counts', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax2.axvline(x=0,c='k',linestyle='dashed')
    ax2.plot(range_show,avg,color='#CC6600')
    ax2.plot(range_show,avg3,color='green')
                
    ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
    ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
    ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)    
    
    ax2.plot(range_show,avg_plus_ste3,color='#227700',linestyle='dashed')
    ax2.plot(range_show,avg_minus_ste3,color='#227700',linestyle='dashed')
    ax2.fill_between(range_show, avg_plus_ste3, avg_minus_ste3, color='#227700',alpha=0.3)
    
    #ax2.legend([SIDE+' codon',title_map_rna1+' AVG(+/- 1STE)',title_map_rna2+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax2.tick_params(axis='y', labelsize=20)
    ax2.tick_params(axis='x', labelsize=20)
    #ax2.set_xlabel('nt', fontsize=20)
    ax2.set_ylim(0, MAX_Yc+0.1*MAX_Yc)
    
    #### codon distribution figure
    avg_sum = sum(avg)
    avg3_sum = sum(avg3)
    avg_dist = [i/avg_sum for i in avg]
    avg3_dist = [i/avg3_sum for i in avg3]

    avg_plus_ste_dist = [i/avg_sum for i in avg_plus_ste]
    avg_minus_ste_dist = [i/avg_sum for i in avg_minus_ste]

    avg_plus_ste3_dist = [i/avg3_sum for i in avg_plus_ste3]
    avg_minus_ste3_dist = [i/avg3_sum for i in avg_minus_ste3] 
    
    avg2_sum = sum(avg2)
    avg4_sum = sum(avg4)
    avg2_dist = [i/avg2_sum for i in avg2]
    avg4_dist = [i/avg4_sum for i in avg4]

    avg_plus_ste2_dist = [i/avg2_sum for i in avg_plus_ste2]
    avg_minus_ste2_dist = [i/avg2_sum for i in avg_minus_ste2]

    avg_plus_ste4_dist = [i/avg4_sum for i in avg_plus_ste4]
    avg_minus_ste4_dist = [i/avg4_sum for i in avg_minus_ste4] 
    
    MAX_Y = max(max(avg_plus_ste_dist),max(avg_plus_ste2_dist), max(avg_plus_ste3_dist), max(avg_plus_ste4_dist))
    
    if rc == 'no_need_rc':
        ylabel = 'number of sites'
    else:
        ylabel = 'read counts'
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax22.set_ylabel(ylabel+' distribution', fontsize=20)
    ax22.axvline(x=0,c='k',linestyle='dashed')
    ax22.plot(range_show,avg_dist,color='#CC6600')
    ax22.plot(range_show,avg3_dist,color='green')
                
    ax22.plot(range_show,avg_plus_ste_dist,color='orange',linestyle='dashed')
    ax22.plot(range_show,avg_minus_ste_dist,color='orange',linestyle='dashed')
    ax22.fill_between(range_show, avg_plus_ste_dist, avg_minus_ste_dist, color='orange',alpha=0.3)    
    
    ax22.plot(range_show,avg_plus_ste3_dist,color='#227700',linestyle='dashed')
    ax22.plot(range_show,avg_minus_ste3_dist,color='#227700',linestyle='dashed')
    ax22.fill_between(range_show, avg_plus_ste3_dist, avg_minus_ste3_dist, color='#227700',alpha=0.3)
    
    #ax22.legend([SIDE+' codon',title_map_rna1+' AVG(+/- 1STE)',title_map_rna2+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax22.tick_params(axis='y', labelsize=20)
    ax22.tick_params(axis='x', labelsize=20)
    ax22.set_xlabel('nt', fontsize=20)
    ax22.set_ylim(0, MAX_Y+0.1*MAX_Y)
    
    # STOP
    SIDE = 'STOP'
    #### old figure  
    ax03.set_title(title_map_gene+'_'+SIDE+str(EXTEND)+' codon')
    #ax03.set_ylabel('number of mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax03.axvline(x=0,c='k',linestyle='dashed')
    ax03.plot(range_show,mRNA_num2,color='#AA7700')
    ax03.plot(range_show,non_zero2,color='#CC6600')
    ax03.plot(range_show,mRNA_num4,color='#668800')
    ax03.plot(range_show,non_zero4,color='green')
    ax03.tick_params(axis='y', labelsize=20)
    ax03.legend([SIDE+' codon',title_map_rna1+' mRNA have position',title_map_rna1+' mRNA have site',
               title_map_rna2+' mRNA have position',title_map_rna2+' mRNA have site'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    #### ratio figure       
    ratio_data2 = [non_zero2[i]/mRNA_num2[i] for i in range(len(mRNA_num2))]
    ratio_data4 = [non_zero4[i]/mRNA_num4[i] for i in range(len(mRNA_num4))]
    #plot mRNA info
    #plt.figure(figsize=(14,8),dpi=150)
    #ax3.set_title(title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    
    #ax3.set_ylabel('ratio mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax3.axvline(x=0,c='k',linestyle='dashed')
    ax3.plot(range_show,ratio_data2,color='#CC6600')
    ax3.plot(range_show,ratio_data4,color='green')
    ax3.tick_params(axis='y', labelsize=20)
    ax3.legend([SIDE+' codon',title_map_rna1+' mRNA site ratio\n(site/position)',title_map_rna2+' mRNA site ratio\n(site/position)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    #### codon figure
    #plot site
    '''if rc == 'no_need_rc':
        ax4.set_ylabel('number of sites',fontsize=20)
    else:
        ax4.set_ylabel('read counts', fontsize=20)
    '''
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax4.axvline(x=0,c='k',linestyle='dashed')
    ax4.plot(range_show,avg2,color='#CC6600')
    ax4.plot(range_show,avg4,color='green')
    ax4.plot(range_show,avg_plus_ste2,color='orange',linestyle='dashed')
    ax4.plot(range_show,avg_minus_ste2,color='orange',linestyle='dashed')
    ax4.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='orange',alpha=0.3)
    
    ax4.plot(range_show,avg_plus_ste4,color='#227700',linestyle='dashed')
    ax4.plot(range_show,avg_minus_ste4,color='#227700',linestyle='dashed')
    ax4.fill_between(range_show, avg_plus_ste4, avg_minus_ste4, color='#227700',alpha=0.3)
    ax4.legend([SIDE+' codon',title_map_rna1+' AVG(+/- 1STE)',title_map_rna2+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax4.tick_params(axis='y', labelsize=20)
    ax4.set_ylim(0, MAX_Yc+0.1*MAX_Yc)
    
    #### codon distribution figure
   
    if rc == 'no_need_rc':
        ylabel = 'number of sites'
    else:
        ylabel = 'read counts'
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    #ax42.set_ylabel('number of sites',fontsize=20)
    ax42.axvline(x=0,c='k',linestyle='dashed')
    ax42.plot(range_show,avg2_dist,color='#CC6600')
    ax42.plot(range_show,avg4_dist,color='green')
    ax42.plot(range_show,avg_plus_ste2_dist,color='orange',linestyle='dashed')
    ax42.plot(range_show,avg_minus_ste2_dist,color='orange',linestyle='dashed')
    ax42.fill_between(range_show, avg_plus_ste2_dist, avg_minus_ste2_dist, color='orange',alpha=0.3)
    
    ax42.plot(range_show,avg_plus_ste4_dist,color='#227700',linestyle='dashed')
    ax42.plot(range_show,avg_minus_ste4_dist,color='#227700',linestyle='dashed')
    ax42.fill_between(range_show, avg_plus_ste4_dist, avg_minus_ste4_dist, color='#227700',alpha=0.3)
    ax42.legend([SIDE+' codon',title_map_rna1+' AVG(+/- 1STE)',title_map_rna2+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax42.tick_params(axis='y', labelsize=20)
    ax42.set_ylim(0, MAX_Y+0.1*MAX_Y)
    
    plt.xlabel('nt',fontsize=20)
    plt.xticks(fontsize=20)    
    plt.tight_layout()
    plt.savefig(res_path)
    tmp2 = pd.DataFrame({'number of mRNA position(start)':mRNA_num,
                         'number of mRNA2 position(start)':mRNA_num3,
                         'number of mRNA site(start)':non_zero,
                         'number of mRNA2 site(start)':non_zero3,
                         'avg(start)':avg,'avg2(start)':avg3,
                         'avg_plus_ste(start)':avg_plus_ste,'avg_plus_ste2(start)':avg_plus_ste3,
                         'avg_minus_ste(start)':avg_minus_ste,'avg_minus_ste2(start)':avg_minus_ste3,
                         'avg_dist_dist(start)':avg_dist,'avg2_dist(start)':avg3_dist,
                         'avg_plus_ste_dist(start)':avg_plus_ste_dist,'avg_plus_ste2_dist(start)':avg_plus_ste3_dist,
                         'avg_minus_ste_dist(start)':avg_minus_ste_dist,'avg_minus_ste2_dist(start)':avg_minus_ste3_dist,
                         
                         'number of mRNA position(stop)':mRNA_num2,
                         'number of mRNA2 position(stop)':mRNA_num4,
                         'number of mRNA site(stop)':non_zero2,
                         'number of mRNA2 site(stastoprt)':non_zero4,
                         'avg(stop)':avg2,'avg2(stop)':avg4,
                         'avg_plus_ste(stop)':avg_plus_ste2,'avg_plus_ste2(stop)':avg_plus_ste4,
                         'avg_minus_ste(stop)':avg_minus_ste2,'avg_minus_ste2(stop)':avg_minus_ste4,
                         'avg_dist(stop)':avg2_dist,'avg2_dist(stop)':avg4_dist,
                         'avg_plus_ste_dist(stop)':avg_plus_ste2_dist,'avg_plus_ste2_dist(stop)':avg_plus_ste4_dist,
                         'avg_minus_ste_dist(stop)':avg_minus_ste2_dist,'avg_minus_ste2_dist(stop)':avg_minus_ste4_dist,
    })
    #tmp2.to_csv(res_path.replace('png', 'csv'),index=False)
    plt.clf()
    plt.close()
    del non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste
    gc.collect()
