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

# main
def tool2_3(data):
    all_mRNA = pd.read_csv(os.path.abspath(__file__+ '/../../../../data/mRNA_WS275_IDtoName.csv'))
    gene_path = data['gene_path']
    title_map_gene = data['title_map_gene']
    title_map_rna = data['title_map_rna']

    for ANALYZE_TYPE in data['ANALYZE_TYPE']:
        print('===={}===='.format(ANALYZE_TYPE))

        if ANALYZE_TYPE == 'TARGET_SCORE':
            for FILE_NAME in data['TARGET_SCORE_FILE']:
                for percent in data['PERCENT']:
                    for rc in data['READ_COUNT_TYPE']:
                        for gene1, gene2 in zip(data['GENE_LIST1'], data['GENE_LIST2']):
                            df = pd.read_csv(os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                            df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])


                            print(FILE_NAME,rc,gene1, gene2,percent)
                            #df_group = mRNA_list(df_withid, gene1)
                            gene_name1, gene_name2 = check_new(data, gene1, gene2)
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene1, gene_name1)
                            # calculate mRNA number gene1
                            show_range = list(np.arange(percent))
                            avg = []
                            med = []
                            ste = []
                            non_zero = []
                            all_data = []
                            for ind in show_range:
                                a = list(df_group['b'+str(ind+1)])
                                avg.append(np.mean(a))
                                med.append(np.median(a))
                                ste.append(stats.sem(a))
                                non_zero.append(len([nz for nz in a if nz != 0]))
                                all_data.append(len(a))

                            avg_plus_ste = [ i+j for i,j in zip(avg,ste)]
                            avg_minus_ste = [ i-j for i,j in zip(avg,ste)]


                            # calculate mRNA number gene2
                            #df_group = mRNA_list(df_withid, gene1)
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene2, gene_name2)

                            avg2 = []
                            med2 = []
                            ste2 = []
                            non_zero2 = []
                            all_data2 = []
                            for ind in show_range:
                                a = list(df_group['b'+str(ind+1)])
                                avg2.append(np.mean(a))
                                med2.append(np.median(a))
                                ste2.append(stats.sem(a))
                                non_zero2.append(len([nz for nz in a if nz != 0]))
                                all_data2.append(len(a))

                            avg_plus_ste2 = [ i+j for i,j in zip(avg2,ste2)]
                            avg_minus_ste2 = [ i-j for i,j in zip(avg2,ste2)]
                            
                            LENGHT1 = max(all_data)    # 方法可議
                            LENGHT2 = max(all_data2)
                            #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/meta_fig/COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(percent)+'.png')
                            res_path = os.path.abspath(__file__+ '/../../../../../static/paper/meta_fig/compare_list/metagene_COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(percent)+rc+'.png')
                            fig ,(ax1,ax2,ax3) = plt.subplots(3,1, sharex=True, figsize=(12,12),dpi=200)
                            # start

                            #plot mRNA info

                            ax1.set_title(title_map_rna[FILE_NAME] + '\n'+title_map_gene[gene1] +' N='+str(LENGHT1) +
                                        '\n' + title_map_gene[gene2] +' N='+str(LENGHT2), fontsize=20)
                            #ax1.set_ylabel('ratio mRNA',fontsize=20)
                            ax1.set_ylabel('mRNA number',fontsize=20)
                            range_show = np.arange(percent)
                            ax1.plot(range_show,non_zero,color='#CC6600')
                            ax1.plot(range_show,non_zero2,color='#227700')
                            #ax1.plot(range_show,ratio_data,color='#CC6600')
                            #ax1.plot(range_show,ratio_data2,color='#227700')

                            ax1.tick_params(axis='y', labelsize=20)
                            ax1.legend([title_map_gene[gene1]+' mRNA site ratio\n(site/position)',title_map_gene[gene2]+' mRNA site ratio\n(site/position)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

                            if rc == 'no_need_rc':
                                y_lab = 'number of site'
                            else:
                                y_lab = 'read counts'
                            ax2.set_ylabel(y_lab,fontsize=20)
                            ax2.plot(range_show,avg,color='#CC6600')
                            ax2.plot(range_show,avg2,color='green')

                            ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
                            ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)

                            ax2.plot(range_show,avg_plus_ste2,color='#227700',linestyle='dashed')
                            ax2.plot(range_show,avg_minus_ste2,color='#227700',linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='#227700',alpha=0.3)

                            ax2.legend([title_map_gene[gene1]+' AVG(+/- 1STE)',title_map_gene[gene2]+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
                            ax2.tick_params(axis='y', labelsize=20)
                            ax2.tick_params(axis='x', labelsize=20)
                            ax2.set_xlabel('region(%)', fontsize=20)
                            
                            # read count distribution
                            avg_sum = sum(avg)
                            avg2_sum = sum(avg2)
                            avg_dist = [i/avg_sum for i in avg]
                            avg2_dist = [i/avg2_sum for i in avg2]

                            avg_plus_ste_dist = [i/avg_sum for i in avg_plus_ste]
                            avg_minus_ste_dist = [i/avg_sum for i in avg_minus_ste]

                            avg_plus_ste2_dist = [i/avg2_sum for i in avg_plus_ste2]
                            avg_minus_ste2_dist = [i/avg2_sum for i in avg_minus_ste2]

                            ax3.set_ylabel(y_lab+' distribution',fontsize=20)
                            ax3.plot(range_show,avg_dist,color='#CC6600')
                            ax3.plot(range_show,avg2_dist,color='green')

                            ax3.plot(range_show,avg_plus_ste_dist,color='orange',linestyle='dashed')
                            ax3.plot(range_show,avg_minus_ste_dist,color='orange',linestyle='dashed')
                            ax3.fill_between(range_show, avg_plus_ste_dist, avg_minus_ste_dist, color='orange',alpha=0.3)

                            ax3.plot(range_show,avg_plus_ste2_dist,color='#227700',linestyle='dashed')
                            ax3.plot(range_show,avg_minus_ste2_dist,color='#227700',linestyle='dashed')
                            ax3.fill_between(range_show, avg_plus_ste2_dist, avg_minus_ste2_dist, color='#227700',alpha=0.3)

                            ax3.legend([title_map_gene[gene1]+' AVG(+/- 1STE)',title_map_gene[gene2]+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
                            ax3.tick_params(axis='y', labelsize=20)
                            ax3.tick_params(axis='x', labelsize=20)
                            ax3.set_xlabel('region(%)', fontsize=20)
                            
                            
                            plt.tight_layout()
                            plt.savefig(res_path)

                            plt.clf()
                            plt.close()
                            del non_zero, avg, avg_plus_ste, avg_minus_ste
                            gc.collect()

        elif ANALYZE_TYPE == 'PIRSCAN':
            for FILE_NAME in data['PIR_FILE']:
                for percent in data['PERCENT']:
                    for rc in data['READ_COUNT_TYPE']:
                        for gene1, gene2 in zip(data['GENE_LIST1'], data['GENE_LIST2']):
                            df = pd.read_csv(os.path.abspath(__file__+ '/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                            df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])


                            print(FILE_NAME,rc,gene1, gene2,percent)
                            gene_name1, gene_name2 = check_new(data, gene1, gene2)
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene1, gene_name1)

                            # calculate mRNA number gene1
                            show_range = list(np.arange(percent))
                            avg = []
                            med = []
                            ste = []
                            non_zero = []
                            all_data = []
                            for ind in show_range:
                                a = list(df_group['b'+str(ind+1)])
                                avg.append(np.mean(a))
                                med.append(np.median(a))
                                ste.append(stats.sem(a))
                                non_zero.append(len([nz for nz in a if nz != 0]))
                                all_data.append(len(a))

                            avg_plus_ste = [ i+j for i,j in zip(avg,ste)]
                            avg_minus_ste = [ i-j for i,j in zip(avg,ste)]


                            # calculate mRNA number gene2
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene2, gene_name2)

                            avg2 = []
                            med2 = []
                            ste2 = []
                            non_zero2 = []
                            all_data2 = []
                            for ind in show_range:
                                a = list(df_group['b'+str(ind+1)])
                                avg2.append(np.mean(a))
                                med2.append(np.median(a))
                                ste2.append(stats.sem(a))
                                non_zero2.append(len([nz for nz in a if nz != 0]))
                                all_data2.append(len(a))

                            avg_plus_ste2 = [ i+j for i,j in zip(avg2,ste2)]
                            avg_minus_ste2 = [ i-j for i,j in zip(avg2,ste2)]
                            
                            LENGHT1 = max(all_data)    # 方法可議
                            LENGHT2 = max(all_data2)
                            #res_path = os.path.abspath(__file__+ '/../../../pirscan_'+rc+'_output/meta_fig/COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(percent)+'.png')
                            res_path = os.path.abspath(__file__+ '/../../../../../static/paper/meta_fig/compare_list/metagene_COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(percent)+rc+'.png')
                            fig ,(ax1,ax2,ax3) = plt.subplots(3,1, sharex=True, figsize=(12,12),dpi=200)
                            # start

                            #plot mRNA info

                            ax1.set_title(title_map_rna[FILE_NAME] + '\n'+title_map_gene[gene1] +' N='+str(LENGHT1) +
                                        '\n' + title_map_gene[gene2] +' N='+str(LENGHT2), fontsize=20)
                            #ax1.set_ylabel('ratio mRNA',fontsize=20)
                            ax1.set_ylabel('mRNA number',fontsize=20)
                            range_show = np.arange(percent)
                            ax1.plot(range_show,non_zero,color='#CC6600')
                            ax1.plot(range_show,non_zero2,color='#227700')
                            #ax1.plot(range_show,ratio_data,color='#CC6600')
                            #ax1.plot(range_show,ratio_data2,color='#227700')

                            ax1.tick_params(axis='y', labelsize=20)
                            ax1.legend([title_map_gene[gene1]+' mRNA site ratio\n(site/position)',title_map_gene[gene2]+' mRNA site ratio\n(site/position)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

                            if rc == 'no_need_rc':
                                y_lab = 'number of site'
                            else:
                                y_lab = 'read counts'
                            ax2.set_ylabel(y_lab,fontsize=20)
                            ax2.plot(range_show,avg,color='#CC6600')
                            ax2.plot(range_show,avg2,color='green')

                            ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
                            ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)

                            ax2.plot(range_show,avg_plus_ste2,color='#227700',linestyle='dashed')
                            ax2.plot(range_show,avg_minus_ste2,color='#227700',linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='#227700',alpha=0.3)

                            ax2.legend([title_map_gene[gene1]+' AVG(+/- 1STE)',title_map_gene[gene2]+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
                            ax2.tick_params(axis='y', labelsize=20)
                            ax2.tick_params(axis='x', labelsize=20)
                            ax2.set_xlabel('region(%)', fontsize=20)
                            
                            # read count distribution
                            avg_sum = sum(avg)
                            avg2_sum = sum(avg2)
                            avg_dist = [i/avg_sum for i in avg]
                            avg2_dist = [i/avg2_sum for i in avg2]

                            avg_plus_ste_dist = [i/avg_sum for i in avg_plus_ste]
                            avg_minus_ste_dist = [i/avg_sum for i in avg_minus_ste]

                            avg_plus_ste2_dist = [i/avg2_sum for i in avg_plus_ste2]
                            avg_minus_ste2_dist = [i/avg2_sum for i in avg_minus_ste2]

                            ax3.set_ylabel(y_lab+' distribution',fontsize=20)
                            ax3.plot(range_show,avg_dist,color='#CC6600')
                            ax3.plot(range_show,avg2_dist,color='green')

                            ax3.plot(range_show,avg_plus_ste_dist,color='orange',linestyle='dashed')
                            ax3.plot(range_show,avg_minus_ste_dist,color='orange',linestyle='dashed')
                            ax3.fill_between(range_show, avg_plus_ste_dist, avg_minus_ste_dist, color='orange',alpha=0.3)

                            ax3.plot(range_show,avg_plus_ste2_dist,color='#227700',linestyle='dashed')
                            ax3.plot(range_show,avg_minus_ste2_dist,color='#227700',linestyle='dashed')
                            ax3.fill_between(range_show, avg_plus_ste2_dist, avg_minus_ste2_dist, color='#227700',alpha=0.3)

                            ax3.legend([title_map_gene[gene1]+' AVG(+/- 1STE)',title_map_gene[gene2]+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
                            ax3.tick_params(axis='y', labelsize=20)
                            ax3.tick_params(axis='x', labelsize=20)
                            ax3.set_xlabel('region(%)', fontsize=20)
                            
                            
                            plt.tight_layout()
                            plt.savefig(res_path)

                            plt.clf()
                            plt.close()
                            del non_zero, avg, avg_plus_ste, avg_minus_ste
                            gc.collect()

        elif ANALYZE_TYPE == 'RNAUP':
            for FILE_NAME in data['RNAUP_FILE']:
                for percent in data['PERCENT']:
                    for rc in data['READ_COUNT_TYPE']:
                        for gene1, gene2 in zip(data['GENE_LIST1'], data['GENE_LIST2']):
                            df = pd.read_csv(os.path.abspath(__file__+ '/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                            df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])


                            print(FILE_NAME,rc,gene1, gene2,percent)
                            gene_name1, gene_name2 = check_new(data, gene1, gene2)
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene1, gene_name1)

                            # calculate mRNA number gene1
                            show_range = list(np.arange(percent))
                            avg = []
                            med = []
                            ste = []
                            non_zero = []
                            all_data = []
                            for ind in show_range:
                                a = list(df_group['b'+str(ind+1)])
                                avg.append(np.mean(a))
                                med.append(np.median(a))
                                ste.append(stats.sem(a))
                                non_zero.append(len([nz for nz in a if nz != 0]))
                                all_data.append(len(a))
                                
                            avg_plus_ste = [ i+j for i,j in zip(avg,ste)]
                            avg_minus_ste = [ i-j for i,j in zip(avg,ste)]


                            # calculate mRNA number gene2
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene2, gene_name2)

                            avg2 = []
                            med2 = []
                            ste2 = []
                            non_zero2 = []
                            all_data2 = []
                            for ind in show_range:
                                a = list(df_group['b'+str(ind+1)])
                                avg2.append(np.mean(a))
                                med2.append(np.median(a))
                                ste2.append(stats.sem(a))
                                non_zero2.append(len([nz for nz in a if nz != 0]))
                                all_data2.append(len(a))
                                
                            avg_plus_ste2 = [ i+j for i,j in zip(avg2,ste2)]
                            avg_minus_ste2 = [ i-j for i,j in zip(avg2,ste2)]
                            
                            LENGHT1 = max(all_data)    # 方法可議
                            LENGHT2 = max(all_data2)
                            #res_path = os.path.abspath(__file__+ '/../../../RNAup_'+rc+'_output/meta_fig/COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(percent)+'.png')
                            res_path = os.path.abspath(__file__+ '/../../../../static/paper/meta_fig/compare_list/metagene_COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(percent)+rc+'.png')
                            fig ,(ax1,ax2,ax3) = plt.subplots(3,1, sharex=True, figsize=(12,12),dpi=200)
                            # start

                            #plot mRNA info

                            ax1.set_title(title_map_rna[FILE_NAME] + '\n'+title_map_gene[gene1] +' N='+str(LENGHT1) +
                                        '\n' + title_map_gene[gene2] +' N='+str(LENGHT2), fontsize=20)
                            #ax1.set_ylabel('ratio mRNA',fontsize=20)
                            ax1.set_ylabel('mRNA number',fontsize=20)
                            range_show = np.arange(percent)
                            ax1.plot(range_show,non_zero,color='#CC6600')
                            ax1.plot(range_show,non_zero2,color='#227700')
                            #ax1.plot(range_show,ratio_data,color='#CC6600')
                            #ax1.plot(range_show,ratio_data2,color='#227700')

                            ax1.tick_params(axis='y', labelsize=20)
                            ax1.legend([title_map_gene[gene1]+' mRNA site ratio\n(site/position)',title_map_gene[gene2]+' mRNA site ratio\n(site/position)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

                            if rc == 'no_need_rc':
                                y_lab = 'number of site'
                            else:
                                y_lab = 'read counts'
                            ax2.set_ylabel(y_lab,fontsize=20)
                            ax2.plot(range_show,avg,color='#CC6600')
                            ax2.plot(range_show,avg2,color='green')

                            ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
                            ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)

                            ax2.plot(range_show,avg_plus_ste2,color='#227700',linestyle='dashed')
                            ax2.plot(range_show,avg_minus_ste2,color='#227700',linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='#227700',alpha=0.3)

                            ax2.legend([title_map_gene[gene1]+' AVG(+/- 1STE)',title_map_gene[gene2]+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
                            ax2.tick_params(axis='y', labelsize=20)
                            ax2.tick_params(axis='x', labelsize=20)
                            ax2.set_xlabel('region(%)', fontsize=20)

                            # read count distribution
                            avg_sum = sum(avg)
                            avg2_sum = sum(avg2)
                            avg_dist = [i/avg_sum for i in avg]
                            avg2_dist = [i/avg2_sum for i in avg2]

                            avg_plus_ste_dist = [i/avg_sum for i in avg_plus_ste]
                            avg_minus_ste_dist = [i/avg_sum for i in avg_minus_ste]

                            avg_plus_ste2_dist = [i/avg2_sum for i in avg_plus_ste2]
                            avg_minus_ste2_dist = [i/avg2_sum for i in avg_minus_ste2]

                            ax3.set_ylabel(y_lab+' distribution',fontsize=20)
                            ax3.plot(range_show,avg_dist,color='#CC6600')
                            ax3.plot(range_show,avg2_dist,color='green')

                            ax3.plot(range_show,avg_plus_ste_dist,color='orange',linestyle='dashed')
                            ax3.plot(range_show,avg_minus_ste_dist,color='orange',linestyle='dashed')
                            ax3.fill_between(range_show, avg_plus_ste_dist, avg_minus_ste_dist, color='orange',alpha=0.3)

                            ax3.plot(range_show,avg_plus_ste2_dist,color='#227700',linestyle='dashed')
                            ax3.plot(range_show,avg_minus_ste2_dist,color='#227700',linestyle='dashed')
                            ax3.fill_between(range_show, avg_plus_ste2_dist, avg_minus_ste2_dist, color='#227700',alpha=0.3)

                            ax3.legend([title_map_gene[gene1]+' AVG(+/- 1STE)',title_map_gene[gene2]+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
                            ax3.tick_params(axis='y', labelsize=20)
                            ax3.tick_params(axis='x', labelsize=20)
                            ax3.set_xlabel('region(%)', fontsize=20)
                            
                            plt.clf()
                            plt.close()
                            del non_zero, avg, avg_plus_ste, avg_minus_ste
                            gc.collect()

        elif ANALYZE_TYPE == 'G22':
            for FILE_NAME in data['G22_FILE']:
                for percent in data['PERCENT']:
                    for rc in data['READ_COUNT_TYPE']:
                        for gene1, gene2 in zip(data['GENE_LIST1'], data['GENE_LIST2']):
                            df = pd.read_csv(os.path.abspath(__file__+ '/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                            df_withid = pd.merge(df, all_mRNA,left_on='Unnamed: 0', right_on='Gene name',how='left')
                            #df_withid['Gene ID'] = df_withid['Gene ID'].apply(lambda x:x.split('=')[1])


                            print(FILE_NAME,rc,gene1, gene2,percent)
                            gene_name1, gene_name2 = check_new(data, gene1, gene2)
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene1, gene_name1)

                            # calculate mRNA number gene1
                            show_range = list(np.arange(percent))
                            avg = []
                            med = []
                            ste = []
                            non_zero = []
                            all_data = []
                            for ind in show_range:
                                a = list(df_group['b'+str(ind+1)])
                                avg.append(np.mean(a))
                                med.append(np.median(a))
                                ste.append(stats.sem(a))
                                non_zero.append(len([nz for nz in a if nz != 0]))
                                all_data.append(len(a))

                            avg_plus_ste = [ i+j for i,j in zip(avg,ste)]
                            avg_minus_ste = [ i-j for i,j in zip(avg,ste)]
                            ratio_data = [non_zero[i]/all_data[i] for i in range(len(all_data))]

                            # calculate mRNA number gene2
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene2, gene_name2)

                            avg2 = []
                            med2 = []
                            ste2 = []
                            non_zero2 = []
                            all_data2 = []
                            for ind in show_range:
                                a = list(df_group['b'+str(ind+1)])
                                avg2.append(np.mean(a))
                                med2.append(np.median(a))
                                ste2.append(stats.sem(a))
                                non_zero2.append(len([nz for nz in a if nz != 0]))
                                all_data2.append(len(a))

                            avg_plus_ste2 = [ i+j for i,j in zip(avg2,ste2)]
                            avg_minus_ste2 = [ i-j for i,j in zip(avg2,ste2)]
                            ratio_data2 = [non_zero2[i]/all_data2[i] for i in range(len(all_data2))]


                            LENGHT1 = max(all_data)    # 方法可議
                            LENGHT2 = max(all_data2)

                            #res_path = os.path.abspath(__file__+ '/../../../22G_output/meta_fig/COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(percent)+'.png')
                            res_path = os.path.abspath(__file__+ '/../../../../../static/paper/meta_fig/compare_list/metagene_COMPARE_'+FILE_NAME+'_G'+str(gene1)+'_'+str(gene2)+'_L'+str(percent)+'.png')
                            fig ,(ax1,ax2,ax3) = plt.subplots(3,1, sharex=True, figsize=(12,12),dpi=200)
                            # start

                            #plot mRNA info

                            ax1.set_title(title_map_rna[FILE_NAME] + '\n'+title_map_gene[gene1] +' N='+str(LENGHT1) +
                                        '\n' + title_map_gene[gene2] +' N='+str(LENGHT2), fontsize=20)
                            #ax1.set_ylabel('ratio mRNA',fontsize=20)
                            ax1.set_ylabel('mRNA number',fontsize=20)
                            range_show = np.arange(percent)
                            ax1.plot(range_show,non_zero,color='#CC6600')
                            ax1.plot(range_show,non_zero2,color='#227700')
                            #ax1.plot(range_show,ratio_data,color='#CC6600')
                            #ax1.plot(range_show,ratio_data2,color='#227700')

                            ax1.tick_params(axis='y', labelsize=20)
                            ax1.legend([title_map_gene[gene1]+' mRNA site ratio\n(site/position)',title_map_gene[gene2]+' mRNA site ratio\n(site/position)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')

                            if rc == 'no_need_rc':
                                y_lab = 'number of site'
                            else:
                                y_lab = 'read counts'
                            ax2.set_ylabel(y_lab,fontsize=20)
                            ax2.plot(range_show,avg,color='#CC6600')
                            ax2.plot(range_show,avg2,color='green')

                            ax2.plot(range_show,avg_plus_ste,color='orange',linestyle='dashed')
                            ax2.plot(range_show,avg_minus_ste,color='orange',linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste, avg_minus_ste, color='orange',alpha=0.3)

                            ax2.plot(range_show,avg_plus_ste2,color='#227700',linestyle='dashed')
                            ax2.plot(range_show,avg_minus_ste2,color='#227700',linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste2, avg_minus_ste2, color='#227700',alpha=0.3)

                            ax2.legend([title_map_gene[gene1]+' AVG(+/- 1STE)',title_map_gene[gene2]+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
                            ax2.tick_params(axis='y', labelsize=20)
                            ax2.tick_params(axis='x', labelsize=20)
                            ax2.set_xlabel('region(%)', fontsize=20)

                            # read count distribution
                            avg_sum = sum(avg)
                            avg2_sum = sum(avg2)
                            avg_dist = [i/avg_sum for i in avg]
                            avg2_dist = [i/avg2_sum for i in avg2]

                            avg_plus_ste_dist = [i/avg_sum for i in avg_plus_ste]
                            avg_minus_ste_dist = [i/avg_sum for i in avg_minus_ste]

                            avg_plus_ste2_dist = [i/avg2_sum for i in avg_plus_ste2]
                            avg_minus_ste2_dist = [i/avg2_sum for i in avg_minus_ste2]

                            ax3.set_ylabel(y_lab+' distribution',fontsize=20)
                            ax3.plot(range_show,avg_dist,color='#CC6600')
                            ax3.plot(range_show,avg2_dist,color='green')

                            ax3.plot(range_show,avg_plus_ste_dist,color='orange',linestyle='dashed')
                            ax3.plot(range_show,avg_minus_ste_dist,color='orange',linestyle='dashed')
                            ax3.fill_between(range_show, avg_plus_ste_dist, avg_minus_ste_dist, color='orange',alpha=0.3)

                            ax3.plot(range_show,avg_plus_ste2_dist,color='#227700',linestyle='dashed')
                            ax3.plot(range_show,avg_minus_ste2_dist,color='#227700',linestyle='dashed')
                            ax3.fill_between(range_show, avg_plus_ste2_dist, avg_minus_ste2_dist, color='#227700',alpha=0.3)

                            ax3.legend([title_map_gene[gene1]+' AVG(+/- 1STE)',title_map_gene[gene2]+' AVG(+/- 1STE)'], fontsize=18, bbox_to_anchor=(1.05, 1.0), loc='upper left')
                            ax3.tick_params(axis='y', labelsize=20)
                            ax3.tick_params(axis='x', labelsize=20)
                            ax3.set_xlabel('region(%)', fontsize=20)
                            
                            plt.tight_layout()
                            plt.savefig(res_path)

                            plt.clf()
                            plt.close()
                            del non_zero, avg, avg_plus_ste, avg_minus_ste
                            gc.collect()
    #del df_withid, df, df_group