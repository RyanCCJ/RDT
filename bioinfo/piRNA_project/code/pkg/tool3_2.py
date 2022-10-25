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
def tool3_2(data):
    all_mRNA = pd.read_csv(os.path.abspath(__file__+ '/../../../../data/mRNA_WS275_IDtoName.csv'))
    gene_path = data['gene_path']
    title_map_gene = data['title_map_gene']
    title_map_rna = data['title_map_rna']

    for ANALYZE_TYPE in data['ANALYZE_TYPE']:
        print('===={}===='.format(ANALYZE_TYPE))
        if ANALYZE_TYPE == 'TARGET_SCORE':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['TARGET_SCORE_FILE']:
                    for EXTEND in data['EXTEND_RANGE']:

                        df = pd.read_csv(os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_START_'+str(EXTEND)+'.csv'))
                        df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                        #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])
                        df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_STOP_'+str(EXTEND)+'.csv'))
                        df_withid2 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                        #df_withid2['Gene ID'] = df_withid2['Gene ID'].apply(lambda x:x.split('=')[1])
                        for gene in data['ANALYZE_GROUP']:#[0,1,2,8]:
                            print(rc, FILE_NAME,EXTEND, gene)
                            gene_name = check_new(data, gene)
                            non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste = calculate_line(df_withid, rc, FILE_NAME, EXTEND, 'START', gene_path, gene, gene_name)
                            non_zero2, mRNA_num2, avg2, med2, avg_plus_ste2, avg_minus_ste2, std2, ste2 = calculate_line(df_withid2, rc, FILE_NAME, EXTEND, 'STOP', gene_path, gene, gene_name)
                            codon_plot_target_score(rc, FILE_NAME, title_map_rna[FILE_NAME], EXTEND, gene, title_map_gene[gene], non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste,
                                    non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2, data['img_title'])

        elif ANALYZE_TYPE == 'PIRSCAN':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['PIR_FILE']:
                    for EXTEND in data['EXTEND_RANGE']:

                        df = pd.read_csv(os.path.abspath(__file__+ '/../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_START_'+str(EXTEND)+'.csv'))
                        df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                        #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])
                        df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_STOP_'+str(EXTEND)+'.csv'))
                        df_withid2 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                        #df_withid2['Gene ID'] = df_withid2['Gene ID'].apply(lambda x:x.split('=')[1])
                        for gene in data['ANALYZE_GROUP']:#[0,1,2,8]:
                            print(rc, FILE_NAME,EXTEND, gene)
                            gene_name = check_new(data, gene)
                            non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste = calculate_line(df_withid, rc, FILE_NAME, EXTEND, 'START', gene_path, gene, gene_name)
                            non_zero2, mRNA_num2, avg2, med2, avg_plus_ste2, avg_minus_ste2, std2, ste2 = calculate_line(df_withid2, rc, FILE_NAME, EXTEND, 'STOP', gene_path, gene, gene_name)
                            codon_plot_pirscan(rc, FILE_NAME, title_map_rna[FILE_NAME], EXTEND, gene, title_map_gene[gene], non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste,
                                    non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2)

        elif ANALYZE_TYPE == 'RNAUP':
            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['RNAUP_FILE']:
                    for EXTEND in data['EXTEND_RANGE']:

                        df = pd.read_csv(os.path.abspath(__file__+ '/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_START_'+str(EXTEND)+'.csv'))
                        df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                        #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])
                        df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_STOP_'+str(EXTEND)+'.csv'))
                        df_withid2 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                        #df_withid2['Gene ID'] = df_withid2['Gene ID'].apply(lambda x:x.split('=')[1])
                        for gene in data['ANALYZE_GROUP']:
                            print(rc, FILE_NAME,EXTEND, gene)
                            gene_name = check_new(data, gene)
                            non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste = calculate_line(df_withid, rc, FILE_NAME, EXTEND, 'START', gene_path, gene, gene_name)
                            non_zero2, mRNA_num2, avg2, med2, avg_plus_ste2, avg_minus_ste2, std2, ste2 = calculate_line(df_withid2, rc, FILE_NAME, EXTEND, 'STOP', gene_path, gene, gene_name)
                            codon_plot_rnaup(rc, FILE_NAME, title_map_rna[FILE_NAME], EXTEND, gene, title_map_gene[gene], non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste,
                                    non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2, data['img_title'])

        elif ANALYZE_TYPE == 'G22':
            for rc in data['READ_COUNT_TYPE']:

                for FILE_NAME in data['G22_FILE']:
                    for EXTEND in data['EXTEND_RANGE']:

                        df = pd.read_csv(os.path.abspath(__file__+ '/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_START_'+str(EXTEND)+'.csv'))
                        df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                        #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])
                        df2 = pd.read_csv(os.path.abspath(__file__+ '/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_STOP_'+str(EXTEND)+'.csv'))
                        df_withid2 = pd.merge(df2, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                        #df_withid2['Gene ID'] = df_withid2['Gene ID'].apply(lambda x:x.split('=')[1])
                        for gene in data['ANALYZE_GROUP']:#[0,1,2,8]:
                            print(rc, FILE_NAME,EXTEND, gene)
                            gene_name = check_new(data, gene)
                            non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste = calculate_line(df_withid, rc, FILE_NAME, EXTEND, 'START', gene_path, gene, gene_name)
                            non_zero2, mRNA_num2, avg2, med2, avg_plus_ste2, avg_minus_ste2, std2, ste2 = calculate_line(df_withid2, rc, FILE_NAME, EXTEND, 'STOP', gene_path, gene, gene_name)
                            codon_plot_g22(rc, FILE_NAME, title_map_rna[FILE_NAME], EXTEND, gene, title_map_gene[gene], non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste,
                                    non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2)
                            
    #del df_withid, df_withid2, df, df2


def calculate_line(df_withid, rc, FILE_NAME, EXTEND, SIDE, gene_path, gene, title_map_gene):

    #df_withid = mRNA_list(df_withid, gene)
    df_withid = add_two_mRNA_list(df_withid, gene_path, gene, title_map_gene)
    # calculate mRNA number
    show_range = list(np.arange(-EXTEND,EXTEND+1))
    mRNA_num = []
    avg = []
    med = []
    ste = []
    std = []
    non_zero = []
    for ind in show_range:
        a = list(df_withid[df_withid[str(ind)]!='N'][str(ind)])
        a = [float(i) for i in a]
        mRNA_num.append(len(a))
        avg.append(np.mean(a))
        med.append(np.median(a))
        std.append(np.std(a))
        ste.append(stats.sem(a))
        non_zero.append(len([nz for nz in a if nz != 0]))
#         if len(a) != 0:
#             avg.append(np.mean(a))
#             med.append(np.median(a))
#             if str(np.std(a)) != 'nan':
#                 std.append(np.std(a))
#             else:
#                 std.append(0)
#             if str(stats.sem(a)) != 'nan':
#                 ste.append(stats.sem(a))
#             else:
#                 ste.append(0)
#             non_zero.append(len([nz for nz in a if nz != 0]))
#         else:
#             avg.append(0)
#             med.append(0)
#             std.append(0)
#             ste.append(0)
#             non_zero.append(0)

    avg_plus_ste = [ i+j for i,j in zip(avg,ste)]
    avg_minus_ste = [ i-j for i,j in zip(avg,ste)]
    return non_zero, mRNA_num, avg, med, avg_plus_ste, avg_minus_ste, std, ste

def codon_plot_target_score(rc, FILE_NAME, title_map_rna, EXTEND, gene, title_map_gene, non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste, non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2, img_title):
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/codon_fig/codon_'+FILE_NAME+'_G'+str(gene)+'_L'+str(EXTEND)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/codon_fig/single/codon_'+FILE_NAME+'_G'+str(gene)+'_L'+str(EXTEND)+rc+'.png')
    fig ,(ax2,ax4) = plt.subplots(1,2, sharex=True, figsize=(15,5),dpi=200)
    # start
    SIDE = 'START'
    #plot mRNA info
    #plt.figure(figsize=(14,8),dpi=150)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    
    #plot site
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/codon_fig/site_'+FILE_NAME+'G'+str(gene)+SIDE+'_L'+str(EXTEND)+'.png')
    
    #plt.figure(figsize=(14,8),dpi=150)
    #ax2.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax2.xlabel('nt',fontsize=20)
    MAX_Y = max(max(avg_plus_ste), max(avg_plus_ste2))
    if rc == 'no_need_rc':
        ax2.set_ylabel('sites / nucleotide bin', fontsize=20)
    else:
        ax2.set_ylabel('read counts / nucleotide bin', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax2.axvline(x=0,c='k',linestyle='dashed')
    ax2.plot(range_show,avg,color='#CC6600')
    ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
    ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
    ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)
    #ax2.legend([SIDE+' codon','AVG.','+1STE.','-1STE.'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax2.tick_params(axis='y', labelsize=20)
    ax2.tick_params(axis='x', labelsize=20)
    ax2.set_xlabel((' '*30)+'Start'+(' '*25)+'(nt)', fontsize=20)
    ax2.set_ylim(0, MAX_Y+0.1*MAX_Y)
    ax2.set_title(img_title, fontsize=25)
    # STOP
    SIDE = 'STOP'
    #plot mRNA info
    #plt.figure(figsize=(14,8),dpi=150)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    
    #plot site
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/codon_fig/site_'+FILE_NAME+'G'+str(gene)+SIDE+'_L'+str(EXTEND)+'.png')
    
    #plt.figure(figsize=(14,8),dpi=150)
    #ax2.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax2.xlabel('nt',fontsize=20)
    if rc == 'no_need_rc':
        ax4.set_ylabel('sites / nucleotide bin', fontsize=20)
    else:
        ax4.set_ylabel('read counts / nucleotide bin', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax4.axvline(x=0,c='k',linestyle='dashed')
    ax4.plot(range_show,avg2,color='#CC6600')
    ax4.plot(range_show,avg_plus_ste2,color='orange',linestyle='dashed')
    ax4.plot(range_show,avg_minus_ste2,color='orange',linestyle='dashed')
    ax4.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='orange',alpha=0.3)
    ax4.tick_params(axis='y', labelsize=20)
    ax4.set_xlabel((' '*30)+'Stop'+(' '*25)+'(nt)', fontsize=20)
    ax4.set_ylim(0, MAX_Y+0.1*MAX_Y)
    ax4.set_title(img_title, fontsize=25)

    plt.xticks(fontsize=20)    
    plt.tight_layout()
    plt.savefig(res_path)
    
    plt.clf()
    plt.close()
    del non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste
    gc.collect()

def codon_plot_pirscan(rc, FILE_NAME, title_map_rna, EXTEND, gene, title_map_gene, non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste, non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2):
    
    #res_path = os.path.abspath(__file__+ '/../../pirscan_'+rc+'_output/codon_fig/codon_'+FILE_NAME+'_G'+str(gene)+'_L'+str(EXTEND)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/codon_fig/single/codon_'+FILE_NAME+'_G'+str(gene)+'_L'+str(EXTEND)+rc+'.png')
    fig ,((ax1,ax3),(ax2,ax4)) = plt.subplots(2,2, sharex=True, figsize=(24,10),dpi=200)
    # start
    SIDE = 'START'
    #plot mRNA info
    #plt.figure(figsize=(14,8),dpi=150)
    ax1.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax1.xlabel('nt',fontsize=20)
    ax1.set_ylabel('mRNA number',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax1.axvline(x=0,c='k',linestyle='dashed')
    ax1.plot(range_show,mRNA_num,color='brown')
    ax1.plot(range_show,non_zero,color='#227700')
    ax1.tick_params(axis='y', labelsize=20)
    #ax1.legend([SIDE+' codon','mRNA have position','mRNA have site'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

    #ax1.savefig(res_path)
    
    #plot site
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/codon_fig/site_'+FILE_NAME+'G'+str(gene)+SIDE+'_L'+str(EXTEND)+'.png')
    
    #plt.figure(figsize=(14,8),dpi=150)
    #ax2.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax2.xlabel('nt',fontsize=20)
    MAX_Y = max(max(avg_plus_ste), max(avg_plus_ste2))
    if rc == 'no_need_rc':
        ax2.set_ylabel('sites number',fontsize=20)
    else:
        ax2.set_ylabel('read counts number', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax2.axvline(x=0,c='k',linestyle='dashed')
    ax2.plot(range_show,avg,color='#CC6600')
    ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
    ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
    ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)
    #ax2.legend([SIDE+' codon','AVG.','+1STE.','-1STE.'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax2.tick_params(axis='y', labelsize=20)
    ax2.tick_params(axis='x', labelsize=20)
    ax2.set_xlabel('nt', fontsize=20)
    ax2.set_ylim(0, MAX_Y+0.1*MAX_Y)
    # STOP
    SIDE = 'STOP'
    #plot mRNA info
    #plt.figure(figsize=(14,8),dpi=150)
    ax3.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax1.xlabel('nt',fontsize=20)
    ax3.set_ylabel('mRNA number',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax3.axvline(x=0,c='k',linestyle='dashed')
    ax3.plot(range_show,mRNA_num2,color='brown')
    ax3.plot(range_show,non_zero2,color='#227700')
    ax3.tick_params(axis='y', labelsize=20)
    ax3.legend([SIDE+' codon','mRNA have position','mRNA have site'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

    #ax1.savefig(res_path)
    
    #plot site
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/codon_fig/site_'+FILE_NAME+'G'+str(gene)+SIDE+'_L'+str(EXTEND)+'.png')
    
    #plt.figure(figsize=(14,8),dpi=150)
    #ax2.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax2.xlabel('nt',fontsize=20)
    if rc == 'no_need_rc':
        ax4.set_ylabel('sites number',fontsize=20)
    else:
        ax4.set_ylabel('read counts number', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax4.axvline(x=0,c='k',linestyle='dashed')
    ax4.plot(range_show,avg2,color='#CC6600')
    ax4.plot(range_show,avg_plus_ste2,color='orange',linestyle='dashed')
    ax4.plot(range_show,avg_minus_ste2,color='orange',linestyle='dashed')
    ax4.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='orange',alpha=0.3)
    ax4.legend([SIDE+' codon','AVG.','+1STE.','-1STE.'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax4.tick_params(axis='y', labelsize=20)
    ax4.set_ylim(0, MAX_Y+0.1*MAX_Y)
    plt.xlabel('nt',fontsize=20)
    plt.xticks(fontsize=20)    
    plt.tight_layout()
    plt.savefig(res_path)
    
    plt.clf()
    plt.close()
    del non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste
    gc.collect()

def codon_plot_rnaup(rc, FILE_NAME, title_map_rna, EXTEND, gene, title_map_gene, non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste, non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2, img_title):
    
    #res_path = os.path.abspath(__file__+ '/../../../RNAup_'+rc+'_output/codon_fig/codon_'+FILE_NAME+'_G'+str(gene)+'_L'+str(EXTEND)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/codon_fig/single/codon_'+FILE_NAME+'_G'+str(gene)+'_L'+str(EXTEND)+rc+'.png')
    fig ,(ax2,ax4) = plt.subplots(1,2, sharex=True, figsize=(15,5),dpi=200)
    # start
    SIDE = 'START'
    #plot mRNA info
    #plt.figure(figsize=(14,8),dpi=150)
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    #plot site
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/codon_fig/site_'+FILE_NAME+'G'+str(gene)+SIDE+'_L'+str(EXTEND)+'.png')
    
    #plt.figure(figsize=(14,8),dpi=150)
    #ax2.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax2.xlabel('nt',fontsize=20)
    MAX_Y = max(max(avg_plus_ste), max(avg_plus_ste2))
    
    if rc == 'no_need_rc':
        ax2.set_ylabel('sites / nucleotide bin',fontsize=20)
    else:
        ax2.set_ylabel('read counts / nucleotide bin', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax2.axvline(x=0,c='k',linestyle='dashed')
    ax2.plot(range_show,avg,color='#CC6600')
    ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
    ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
    ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)
    #ax2.legend([SIDE+' codon','AVG.','+1STE.','-1STE.'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax2.tick_params(axis='y', labelsize=20)
    ax2.tick_params(axis='x', labelsize=20)
    ax2.set_xlabel((' '*30)+'Start'+(' '*25)+'(nt)', fontsize=20)
    ax2.set_ylim(0, MAX_Y+0.1*MAX_Y)
    ax2.set_title(img_title, fontsize=30)
    # STOP
    SIDE = 'STOP'
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    #plot site
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/codon_fig/site_'+FILE_NAME+'G'+str(gene)+SIDE+'_L'+str(EXTEND)+'.png')
    
    #plt.figure(figsize=(14,8),dpi=150)
    #ax2.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax2.xlabel('nt',fontsize=20)
    if rc == 'no_need_rc':
        ax4.set_ylabel('sites / nucleotide bin',fontsize=20)
    else:
        ax4.set_ylabel('read counts / nucleotide bin', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax4.axvline(x=0,c='k',linestyle='dashed')
    ax4.plot(range_show,avg2,color='#CC6600')
    ax4.plot(range_show,avg_plus_ste2,color='orange',linestyle='dashed')
    ax4.plot(range_show,avg_minus_ste2,color='orange',linestyle='dashed')
    ax4.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='orange',alpha=0.3)
    ax4.tick_params(axis='y', labelsize=20)
    ax4.set_ylim(0, MAX_Y+0.1*MAX_Y)
    ax4.set_title(img_title, fontsize=30)

    plt.xlabel((' '*30)+'Stop'+(' '*25)+'(nt)',fontsize=20)
    plt.xticks(fontsize=20)    
    plt.tight_layout()
    plt.savefig(res_path)
    
    plt.clf()
    plt.close()
    del non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste
    gc.collect()   

def codon_plot_g22(rc, FILE_NAME, title_map_rna, EXTEND, gene, title_map_gene, non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste, non_zero2, mRNA_num2, avg2, avg_plus_ste2, avg_minus_ste2):
    
    #res_path = os.path.abspath(__file__+ '/../../../22G_output/codon_fig/codon_'+FILE_NAME+'_G'+str(gene)+'_L'+str(EXTEND)+'.png')
    res_path = os.path.abspath(__file__+ '/../../../../../static/paper/codon_fig/single/codon_'+FILE_NAME+'_G'+str(gene)+'_L'+str(EXTEND)+'.png')
    fig ,((ax1,ax3),(ax2,ax4)) = plt.subplots(2,2, sharex=True, figsize=(24,10),dpi=200)
    # start
    SIDE = 'START'
    #ratio_data1 = [non_zero[i]/mRNA_num[i] for i in range(len(mRNA_num))]
    #ratio_data2 = [non_zero2[i]/mRNA_num2[i] for i in range(len(mRNA_num2))]
    
    
    #plot mRNA info
    #plt.figure(figsize=(14,8),dpi=150)
    ax1.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax1.xlabel('nt',fontsize=20)
    ax1.set_ylabel('number of mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax1.axvline(x=0,c='k',linestyle='dashed')
    ax1.plot(range_show,mRNA_num,color='brown')
    ax1.plot(range_show,non_zero,color='#227700')
    #ax1.plot(range_show,ratio_data1,color='#227700') # ratio
    ax1.tick_params(axis='y', labelsize=20)
    #ax1.legend([SIDE+' codon','mRNA site ratio(position/site)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

    #ax1.savefig(res_path)
    
    #plot site
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/codon_fig/site_'+FILE_NAME+'G'+str(gene)+SIDE+'_L'+str(EXTEND)+'.png')
    
    #plt.figure(figsize=(14,8),dpi=150)
    #ax2.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax2.xlabel('nt',fontsize=20)
    MAX_Y = max(max(avg_plus_ste), max(avg_plus_ste2))
    if rc == 'no_need_rc':
        ax2.set_ylabel('number of sites',fontsize=20)
    else:
        ax2.set_ylabel('read counts', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax2.axvline(x=0,c='k',linestyle='dashed')
    ax2.plot(range_show,avg,color='#CC6600')
    ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
    ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
    ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)
    #ax2.legend([SIDE+' codon','AVG.','+1STE.','-1STE.'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax2.tick_params(axis='y', labelsize=20)
    ax2.tick_params(axis='x', labelsize=20)
    ax2.set_xlabel('nt', fontsize=20)
    ax2.set_ylim(0, MAX_Y+0.1*MAX_Y)
    # STOP
    SIDE = 'STOP'
    #plot mRNA info
    #plt.figure(figsize=(14,8),dpi=150)
    ax3.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax1.xlabel('nt',fontsize=20)
    ax3.set_ylabel('number of mRNA',fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax3.axvline(x=0,c='k',linestyle='dashed')
    ax3.plot(range_show,mRNA_num2,color='brown')
    ax3.plot(range_show,non_zero2,color='#227700')
    #ax3.plot(range_show,ratio_data2,color='#227700') # ratio
    ax3.tick_params(axis='y', labelsize=20)
    ax3.legend([SIDE+' codon','mRNA site ratio(position/site)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

    #ax1.savefig(res_path)
    
    #plot site
    #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/codon_fig/site_'+FILE_NAME+'G'+str(gene)+SIDE+'_L'+str(EXTEND)+'.png')
    
    #plt.figure(figsize=(14,8),dpi=150)
    #ax2.set_title(title_map_rna+'_'+title_map_gene+'_'+SIDE+str(EXTEND)+' codon', fontsize=20)
    #ax2.xlabel('nt',fontsize=20)
    if rc == 'no_need_rc':
        ax4.set_ylabel('number of sites',fontsize=20)
    else:
        ax4.set_ylabel('read counts', fontsize=20)
    
    range_show = np.arange(-1*EXTEND,EXTEND+1)
    ax4.axvline(x=0,c='k',linestyle='dashed')
    ax4.plot(range_show,avg2,color='#CC6600')
    ax4.plot(range_show,avg_plus_ste2,color='orange',linestyle='dashed')
    ax4.plot(range_show,avg_minus_ste2,color='orange',linestyle='dashed')
    ax4.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='orange',alpha=0.3)
    ax4.legend([SIDE+' codon','AVG.','+1STE.','-1STE.'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax4.tick_params(axis='y', labelsize=20)
    ax4.set_ylim(0, MAX_Y+0.1*MAX_Y)
    plt.xlabel('nt',fontsize=20)
    plt.xticks(fontsize=20)    
    plt.tight_layout()
    #plt.show()
    plt.savefig(res_path)
    
    plt.clf()
    plt.close()
    del non_zero, mRNA_num, avg, avg_plus_ste, avg_minus_ste
    gc.collect()

