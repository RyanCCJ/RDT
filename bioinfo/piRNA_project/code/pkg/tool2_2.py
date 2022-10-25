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

def check_new(data, gene):
    tran_name = 'new_transcript'
    if gene == 0:
        gene_name = data['NEW_ANALYZE_GROUP']
        data['title_map_gene'][gene] = tran_name
    else:
        gene_name = data['title_map_gene'][gene]
    return gene_name

# main
def tool2_2(data):
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
                        df = pd.read_csv(os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/' +
                                        FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                        df_withid = pd.merge(
                            df, all_mRNA, left_on='Unnamed: 0', right_on='Gene name', how='left')
                        #df_withid['Gene ID'] = df_withid['Gene ID'].apply(
                        #    lambda x: x.split('=')[1])

                        for gene in data['ANALYZE_GROUP']:  # [1,2,8]:
                            print(FILE_NAME, rc, gene, percent)
                            #df_group = mRNA_list(df_withid, gene)
                            gene_name = check_new(data, gene)
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene, gene_name)

                            # calculate mRNA number
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

                            avg_plus_ste = [i+j for i, j in zip(avg, ste)]
                            avg_minus_ste = [i-j for i, j in zip(avg, ste)]
                            ratio_data = [non_zero[i]/all_data[i] for i in range(len(all_data))]

                            #res_path = os.path.abspath(__file__+ '/../../../target_score_'+rc+'_output/meta_fig/metagene_' + \
                            #    FILE_NAME+'_G'+str(gene)+'_L'+str(percent)+'.png')
                            res_path = os.path.abspath(__file__+ '/../../../../../static/paper/meta_fig/single/metagene_' + \
                                FILE_NAME+'_G'+str(gene)+'_L'+str(percent)+rc+'.png')
                            fig, (ax2) = plt.subplots(1, 1, sharex=True, figsize=(8, 5), dpi=200)
                            # start
                            range_show = np.arange(percent)
                            
                            if rc == 'no_need_rc':
                                ax2.set_ylabel('sites / bin', fontsize=20)
                            else:
                                ax2.set_ylabel('read counts / bin', fontsize=20)

                            ax2.plot(range_show, avg, color='#CC6600')
                            ax2.plot(range_show, avg_plus_ste,
                                    color='orange', linestyle='dashed')
                            ax2.plot(range_show, avg_minus_ste,
                                    color='orange', linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste,
                                            avg_minus_ste, color='orange', alpha=0.3)
                            ax2.tick_params(axis='y', labelsize=20)
                            ax2.tick_params(axis='x', labelsize=20)
                            ax2.set_xlabel('region(%)', fontsize=20)
                            ax2.set_title(data['img_title'], fontsize=25)

                            plt.tight_layout()
                            plt.savefig(res_path)
                            tmp2 = pd.DataFrame({'number of mRNA':pd.Series(ratio_data),
                                                'avg':pd.Series(avg),
                                                'avg_plus_ste':pd.Series(avg_plus_ste),
                                                'avg_minus_ste':pd.Series(avg_minus_ste)})
                            #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
                            plt.clf()
                            plt.close()
                            del non_zero, avg, avg_plus_ste, avg_minus_ste
                            gc.collect()

        elif ANALYZE_TYPE == 'PIRSCAN':
            for FILE_NAME in data['PIR_FILE']:
                for percent in data['PERCENT']:
                    for rc in data['READ_COUNT_TYPE']:
                        df = pd.read_csv(
                            os.path.abspath(__file__+ '/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                        df_withid = pd.merge(
                            df, all_mRNA, left_on='Unnamed: 0', right_on='Gene name', how='left')
                        #df_withid['Gene ID'] = df_withid['Gene ID'].apply(
                        #    lambda x: x.split('=')[1])

                        for gene in data['ANALYZE_GROUP']:
                            print(FILE_NAME, rc, gene, percent)
                            gene_name = check_new(data, gene)
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene, gene_name)

                            # calculate mRNA number
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

                            avg_plus_ste = [i+j for i, j in zip(avg, ste)]
                            avg_minus_ste = [i-j for i, j in zip(avg, ste)]
                            ratio_data = [non_zero[i]/all_data[i] for i in range(len(all_data))]
                            
                            #res_path = os.path.abspath(__file__+ '/../../../pirscan_'+rc+'_output/meta_fig/metagene_' + \
                            #    FILE_NAME+'_G'+str(gene)+'_L'+str(percent)+'.png')
                            res_path = os.path.abspath(__file__+ '/../../../../../static/paper/meta_fig/single/metagene_' + \
                                FILE_NAME+'_G'+str(gene)+'_L'+str(percent)+rc+'.png')
                            fig, (ax2) = plt.subplots(1, 1, sharex=True, figsize=(8, 5), dpi=200)
                            # start

                            # plot mRNA info

                            range_show = np.arange(percent)

                            if rc == 'no_need_rc':
                                ax2.set_ylabel('sites / bin', fontsize=20)
                            else:
                                ax2.set_ylabel('read counts / bin', fontsize=20)

                            ax2.plot(range_show, avg, color='#CC6600')
                            ax2.plot(range_show, avg_plus_ste,
                                    color='orange', linestyle='dashed')
                            ax2.plot(range_show, avg_minus_ste,
                                    color='orange', linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste,
                                            avg_minus_ste, color='orange', alpha=0.3)
                            ax2.tick_params(axis='y', labelsize=20)
                            ax2.tick_params(axis='x', labelsize=20)
                            ax2.set_xlabel('region(%)', fontsize=20)
                            ax2.set_title(data['img_title'], fontsize=25)

                            plt.tight_layout()
                            plt.savefig(res_path)
                            tmp2 = pd.DataFrame({'number of mRNA':pd.Series(ratio_data),
                                                'avg':pd.Series(avg),
                                                'avg_plus_ste':pd.Series(avg_plus_ste),
                                                'avg_minus_ste':pd.Series(avg_minus_ste)})
                            #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
                            plt.clf()
                            plt.close()
                            del non_zero, avg, avg_plus_ste, avg_minus_ste
                            gc.collect()

        elif ANALYZE_TYPE == 'RNAUP':
            for FILE_NAME in data['RNAUP_FILE']:
                for percent in data['PERCENT']:
                    for rc in data['READ_COUNT_TYPE']:
                        df = pd.read_csv(
                            os.path.abspath(__file__+ '/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                        df_withid = pd.merge(
                            df, all_mRNA, left_on='Unnamed: 0', right_on='Gene name', how='left')
                        #df_withid['Gene ID'] = df_withid['Gene ID'].apply(
                        #    lambda x: x.split('=')[1])

                        for gene in data['ANALYZE_GROUP']:
                            print(FILE_NAME, rc, gene, percent)
                            gene_name = check_new(data, gene)
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene, gene_name)

                            # calculate mRNA number
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
                                
                            avg_plus_ste = [i+j for i, j in zip(avg, ste)]
                            avg_minus_ste = [i-j for i, j in zip(avg, ste)]
                            ratio_data = [non_zero[i]/all_data[i] for i in range(len(all_data))]
                            
                            #res_path = os.path.abspath(__file__+ '/../../../RNAup_'+rc+'_output/meta_fig/metagene_' + \
                            #     FILE_NAME+'_G'+str(gene)+'_L'+str(percent)+'.png')
                            res_path = os.path.abspath(__file__+ '/../../../../../static/paper/meta_fig/single/metagene_' + \
                                FILE_NAME+'_G'+str(gene)+'_L'+str(percent)+rc+'.png')
                            fig, ax2 = plt.subplots(1, 1, sharex=True, figsize=(8, 5), dpi=200)

                            range_show = np.arange(percent)

                            if rc == 'no_need_rc':
                                ax2.set_ylabel('sites / bin', fontsize=20)
                            else:
                                ax2.set_ylabel('read counts number', fontsize=20)

                            ax2.plot(range_show, avg, color='#CC6600')
                            ax2.plot(range_show, avg_plus_ste,
                                    color='orange', linestyle='dashed')
                            ax2.plot(range_show, avg_minus_ste,
                                    color='orange', linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste,
                                            avg_minus_ste, color='orange', alpha=0.3)
                            ax2.tick_params(axis='y', labelsize=20)
                            ax2.tick_params(axis='x', labelsize=20)
                            ax2.set_xlabel('region(%)', fontsize=20)
                            ax2.set_title(data['img_title'], fontsize=25)

                            plt.tight_layout()
                            plt.savefig(res_path)
                            tmp2 = pd.DataFrame({'number of mRNA':pd.Series(ratio_data),
                                                'avg':pd.Series(avg),
                                                'avg_plus_ste':pd.Series(avg_plus_ste),
                                                'avg_minus_ste':pd.Series(avg_minus_ste)})
                            #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
                            plt.clf()
                            plt.close()
                            del non_zero, avg, avg_plus_ste, avg_minus_ste
                            gc.collect()

        elif ANALYZE_TYPE == 'G22':
            for FILE_NAME in data['G22_FILE']:
                for percent in data['PERCENT']:
                    for rc in data['READ_COUNT_TYPE']:
                        df = pd.read_csv(os.path.abspath(__file__+ '/../../../22G_output/'+FILE_NAME +
                                        '_all_mRNA_tool2_'+str(percent)+'.csv'))
                        df_withid = pd.merge(
                            df, all_mRNA, left_on='Unnamed: 0', right_on='Gene name', how='left')
                        #df_withid['Gene ID'] = df_withid['Gene ID'].apply(
                        #    lambda x: x.split('=')[1])

                        for gene in data['ANALYZE_GROUP']:  # [1,2,8]:
                            print(FILE_NAME, rc, gene, percent)
                            gene_name = check_new(data, gene)
                            df_group = add_two_mRNA_list(df_withid, gene_path, gene, gene_name)

                            # calculate mRNA number
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

                            avg_plus_ste = [i+j for i, j in zip(avg, ste)]
                            avg_minus_ste = [i-j for i, j in zip(avg, ste)]
                            ratio_data = [non_zero[i]/all_data[i] for i in range(len(all_data))]
                            
                            #res_path = os.path.abspath(__file__+ '/../../../22G_output/meta_fig/metagene_'+FILE_NAME+'_G'+str(gene)+'_L'+str(percent)+'.png')
                            res_path = os.path.abspath(__file__+ '/../../../../../static/paper/meta_fig/single/metagene_'+FILE_NAME+'_G'+str(gene)+'_L'+str(percent)+'.png')
                            fig, (ax2) = plt.subplots(1, 1, sharex=True, figsize=(8, 5), dpi=200)
                            # start

                            range_show = np.arange(percent)
                            
                            if rc == 'no_need_rc':
                                ax2.set_ylabel('sites / bin', fontsize=20)
                            else:
                                ax2.set_ylabel('read counts / bin', fontsize=20)

                            ax2.plot(range_show, avg, color='#CC6600')
                            ax2.plot(range_show, avg_plus_ste,
                                    color='orange', linestyle='dashed')
                            ax2.plot(range_show, avg_minus_ste,
                                    color='orange', linestyle='dashed')
                            ax2.fill_between(range_show, avg_plus_ste,
                                            avg_minus_ste, color='orange', alpha=0.3)

                            ax2.tick_params(axis='y', labelsize=20)
                            ax2.tick_params(axis='x', labelsize=20)
                            ax2.set_xlabel('region(%)', fontsize=20)
                            ax2.set_title(data['img_title'], fontsize=25)
                            
                            tmp2 = pd.DataFrame({'number of mRNA':pd.Series(ratio_data),
                                                'avg':pd.Series(avg),
                                                'avg_plus_ste':pd.Series(avg_plus_ste),
                                                'avg_minus_ste':pd.Series(avg_minus_ste)})
                            #tmp2.to_csv(res_path.replace('png', 'csv'), index=False)
                            plt.tight_layout()
                            plt.savefig(res_path)
                            #plt.show()
                            plt.clf()
                            plt.close()
                            del non_zero, avg, avg_plus_ste, avg_minus_ste
                            gc.collect()

    #del df_withid, df, df_group