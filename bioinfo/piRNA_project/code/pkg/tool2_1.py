import math
import os
import numpy as np
import pandas as pd

# main
def tool2_1(data):
    percent = data['PERCENT'][0]
    all_mRNA = pd.read_csv(os.path.abspath(__file__+'/../../../../data/mRNA_WS275_WITHregion_v3.csv'))
    all_mRNA = all_mRNA[['Gene ID','sequence_length','Gene name','CDS start','CDS end']]
    all_mRNA['len5'] = all_mRNA['CDS start'] - 1
    all_mRNA['lencds'] = all_mRNA['CDS end'] - all_mRNA['CDS start'] + 1
    all_mRNA['len3'] = all_mRNA['sequence_length'] - all_mRNA['CDS end']
    output = []

    for ANALYZE_TYPE in data['ANALYZE_TYPE']:
        #print('===={}===='.format(ANALYZE_TYPE))
        if ANALYZE_TYPE =='TARGET_SCORE':
            for FILE_NAME in data['TARGET_SCORE_FILE']:
                for rc in data['READ_COUNT_TYPE']:
                    # target score
                    if data['DATA'] == 1:
                        df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/target_score_'+FILE_NAME+'_'+rc+'norm_2step.csv'))
                        df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc']]
                        #df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc','xgu_inseed', 'gu_inseed', 'xgu_innon-seed',
                        #'gu_innon-seed', 'totalmismatch', 'xGU_mispos', 'GU_mispos']]
                        df['middle_pos'] = df['target_score_pos'].apply(lambda x: int((int(x.split('-')[0])+int(x.split('-')[1]))/2))
                        df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                        df_withid = pd.merge(df, all_mRNA,left_on='transcript_name', right_on='Gene name',how='left')

                        df_withid['big_rel'] = df_withid['middle_pos'] / df_withid['sequence_length']
                        df_withid = df_withid[df_withid['big_rel'] <= 1].reset_index(drop=True)

                        all_mRNA_meta = []
                        single_mRNA_meta = []
                        mRNA_name = []
                        M = 0


                        for mRNA, mlen in zip(list(all_mRNA['Gene name']),list(all_mRNA['sequence_length'])):

                            print(M, end='\r')
                            mRNA_name.append(mRNA)
                            tmp = df_withid[df_withid['transcript_name'].isin([mRNA])]

                            single_mRNA_meta += [0] * percent

                            for i,row in tmp.iterrows():
                                pos_bin = math.ceil(row['big_rel']*100 * (percent/100) )  
                                single_mRNA_meta[pos_bin-1] += round(row['norm_rc'],3) #因為第一項是[0]

                            single_mRNA_meta.append(mlen/percent)    
                            all_mRNA_meta.append(single_mRNA_meta)   
                            single_mRNA_meta = []
                            M += 1

                        df = pd.DataFrame(all_mRNA_meta, columns =['b'+str(k+1) for k in range(100)]+['bin_size'], index=mRNA_name)
                        df.to_csv(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                        output.append(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))

                    # CROSS CONDITION 
                    if data['DATA'] == 2:
                        df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/target_score_'+FILE_NAME+'_'+rc+'norm_usetRNA_2step.csv'))
                        df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc']]
                        #df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc','xgu_inseed', 'gu_inseed', 'xgu_innon-seed',
                        #'gu_innon-seed', 'totalmismatch', 'xGU_mispos', 'GU_mispos']]
                        df['middle_pos'] = df['target_score_pos'].apply(lambda x: int((int(x.split('-')[0])+int(x.split('-')[1]))/2))
                        df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                        df_withid = pd.merge(df, all_mRNA,left_on='transcript_name', right_on='Gene name',how='left')

                        df_withid['big_rel'] = df_withid['middle_pos'] / df_withid['sequence_length']
                        df_withid = df_withid[df_withid['big_rel'] <= 1].reset_index(drop=True)

                        all_mRNA_meta = []
                        single_mRNA_meta = []
                        mRNA_name = []
                        M = 0


                        for mRNA, mlen in zip(list(all_mRNA['Gene name']),list(all_mRNA['sequence_length'])):

                            print(M, end='\r')
                            mRNA_name.append(mRNA)
                            tmp = df_withid[df_withid['transcript_name'].isin([mRNA])]

                            single_mRNA_meta += [0] * percent

                            for i,row in tmp.iterrows():
                                pos_bin = math.ceil(row['big_rel']*100 * (percent/100) )  
                                single_mRNA_meta[pos_bin-1] += round(row['norm_rc'],3) #因為第一項是[0]

                            single_mRNA_meta.append(mlen/percent)    
                            all_mRNA_meta.append(single_mRNA_meta)   
                            single_mRNA_meta = []
                            M += 1

                        df = pd.DataFrame(all_mRNA_meta, columns =['b'+str(k+1) for k in range(100)]+['bin_size'], index=mRNA_name)
                        df.to_csv(os.path.abspath(__file__+'/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                        output.append(os.path.abspath(__file__+'/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))

        # pirscan
        elif ANALYZE_TYPE == 'PIRSCAN':
            for FILE_NAME in data['PIR_FILE']:
                for rc in data['READ_COUNT_TYPE']:
                    #print(FILE_NAME, rc)

                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/'+FILE_NAME+'_pirscan_withup.csv'))
                    df = df.rename(columns={'rem_tran_target_pos':'target_score_pos','score':'targeting_score'})
                    df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc']]
                    df['middle_pos'] = df['target_score_pos'].apply(lambda x: int((int(x.split('-')[0])+int(x.split('-')[1]))/2))
                    df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                    df_withid = pd.merge(df, all_mRNA,left_on='transcript_name', right_on='Gene name',how='left')

                    df_withid['big_rel'] = df_withid['middle_pos'] / df_withid['sequence_length']
                    df_withid = df_withid[df_withid['big_rel'] <= 1].reset_index(drop=True)

                    all_mRNA_meta = []
                    single_mRNA_meta = []
                    mRNA_name = []
                    M = 0


                    for mRNA, mlen in zip(list(all_mRNA['Gene name']),list(all_mRNA['sequence_length'])):

                        print(M, end='\r')
                        mRNA_name.append(mRNA)
                        tmp = df_withid[df_withid['transcript_name'].isin([mRNA])]

                        single_mRNA_meta += [0] * percent

                        for i,row in tmp.iterrows():
                            pos_bin = math.ceil(row['big_rel']*100 * (percent/100) )  
                            single_mRNA_meta[pos_bin-1] += round(row['norm_rc'],3) #因為第一項是[0]

                        single_mRNA_meta.append(mlen/percent)    
                        all_mRNA_meta.append(single_mRNA_meta)   
                        single_mRNA_meta = []
                        M += 1

                    df = pd.DataFrame(all_mRNA_meta, columns =['b'+str(k+1) for k in range(100)]+['bin_size'], index=mRNA_name)
                    df.to_csv(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                    output.append(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))

        # rnaup
        elif ANALYZE_TYPE == 'RNAUP':
            for FILE_NAME in data['RNAUP_FILE']:
                for rc in data['READ_COUNT_TYPE']:
                    #print(FILE_NAME, rc)

                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/RNAup_'+FILE_NAME+'_'+rc+'norm_2step.csv'))
                    df = df[['RNAup_pos','RNAup_score','regulator_name','transcript_name','norm_rc']]
                    df['middle_pos'] = df['RNAup_pos'].apply(lambda x: int((int(x.split('-')[0])+int(x.split('-')[1]))/2))
                    df['site_len'] = df['RNAup_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                    df_withid = pd.merge(df, all_mRNA,left_on='transcript_name', right_on='Gene name',how='left')

                    df_withid['big_rel'] = df_withid['middle_pos'] / df_withid['sequence_length']
                    df_withid = df_withid[df_withid['big_rel'] <= 1].reset_index(drop=True)

                    all_mRNA_meta = []
                    single_mRNA_meta = []
                    mRNA_name = []
                    M = 0


                    for mRNA, mlen in zip(list(all_mRNA['Gene name']),list(all_mRNA['sequence_length'])):

                        print(M, end='\r')
                        mRNA_name.append(mRNA)
                        tmp = df_withid[df_withid['transcript_name'].isin([mRNA])]

                        single_mRNA_meta += [0] * percent

                        for i,row in tmp.iterrows():
                            pos_bin = math.ceil(row['big_rel']*100 * (percent/100) )  
                            single_mRNA_meta[pos_bin-1] += round(row['norm_rc'],3) #因為第一項是[0]

                        single_mRNA_meta.append(mlen/percent)    
                        all_mRNA_meta.append(single_mRNA_meta)   
                        single_mRNA_meta = []
                        M += 1

                    df = pd.DataFrame(all_mRNA_meta, columns =['b'+str(k+1) for k in range(100)]+['bin_size'], index=mRNA_name)
                    df.to_csv(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                    output.append(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))

        # G22
        elif ANALYZE_TYPE == 'G22':
            for FILE_NAME in data['G22_FILE']:
                if data['DATA'] == 1:
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/22G/'+FILE_NAME+'_'+'norm_1step.csv'))
                    df = df.rename(columns={'evenly_rc':'norm_rc','rem_tran_target_pos':'target_score_pos','ref_id':'transcript_name','input_id':'regulator_name'})
                    df = df[['norm_rc','target_score_pos','regulator_name','transcript_name']]
                    df['middle_pos'] = df['target_score_pos'].apply(lambda x: int((int(x.split('-')[0])+int(x.split('-')[1]))/2))
                    df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                    df_withid = pd.merge(df, all_mRNA,left_on='transcript_name', right_on='Gene name',how='left')

                    df_withid['big_rel'] = df_withid['middle_pos'] / df_withid['sequence_length']
                    df_withid = df_withid[df_withid['big_rel'] <= 1].reset_index(drop=True)

                    all_mRNA_meta = []
                    single_mRNA_meta = []
                    mRNA_name = []
                    M = 0


                    for mRNA, mlen in zip(list(all_mRNA['Gene name']),list(all_mRNA['sequence_length'])):

                        print(M, end='\r')
                        mRNA_name.append(mRNA)
                        tmp = df_withid[df_withid['transcript_name'].isin([mRNA])]

                        single_mRNA_meta += [0] * percent

                        for i,row in tmp.iterrows():
                            pos_bin = math.ceil(row['big_rel']*100 * (percent/100) )  
                            single_mRNA_meta[pos_bin-1] += round(row['norm_rc'],3) #因為第一項是[0]

                        single_mRNA_meta.append(mlen/percent)    
                        all_mRNA_meta.append(single_mRNA_meta)   
                        single_mRNA_meta = []
                        M += 1

                    df = pd.DataFrame(all_mRNA_meta, columns =['b'+str(k+1) for k in range(100)]+['bin_size'], index=mRNA_name)
                    df.to_csv(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                    output.append(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))

                # NOR_G22
                if data['DATA'] == 2:
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/22G/'+FILE_NAME+'_'+'EGL17M0_usemiRNA_norm_1step.csv'))  #### piRNA -> 21
                    df = df.rename(columns={'evenly_rc':'norm_rc','rem_tran_target_pos':'target_score_pos','ref_id':'transcript_name','input_id':'regulator_name'})
                    df = df[['norm_rc','target_score_pos','regulator_name','transcript_name']]
                    df['middle_pos'] = df['target_score_pos'].apply(lambda x: int((int(x.split('-')[0])+int(x.split('-')[1]))/2))
                    df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                    df_withid = pd.merge(df, all_mRNA,left_on='transcript_name', right_on='Gene name',how='left')

                    df_withid['big_rel'] = df_withid['middle_pos'] / df_withid['sequence_length']
                    df_withid = df_withid[df_withid['big_rel'] <= 1].reset_index(drop=True)

                    all_mRNA_meta = []
                    single_mRNA_meta = []
                    mRNA_name = []
                    M = 0


                    for mRNA, mlen in zip(list(all_mRNA['Gene name']),list(all_mRNA['sequence_length'])):

                        print(M, end='\r')
                        mRNA_name.append(mRNA)
                        tmp = df_withid[df_withid['transcript_name'].isin([mRNA])]

                        single_mRNA_meta += [0] * percent

                        for i,row in tmp.iterrows():
                            pos_bin = math.ceil(row['big_rel']*100 * (percent/100) )  
                            single_mRNA_meta[pos_bin-1] += round(row['norm_rc'],3) #因為第一項是[0]

                        single_mRNA_meta.append(mlen/percent)    
                        all_mRNA_meta.append(single_mRNA_meta)   
                        single_mRNA_meta = []
                        M += 1

                    df = pd.DataFrame(all_mRNA_meta, columns =['b'+str(k+1) for k in range(100)]+['bin_size'], index=mRNA_name)
                    df.to_csv(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
                    output.append(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool2_'+str(percent)+'.csv'))
    return output