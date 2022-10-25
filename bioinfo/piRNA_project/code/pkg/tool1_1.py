import os
import pandas as pd

# main
def tool1_1(data):
    all_mRNA = pd.read_csv(os.path.abspath(__file__+'/../../../../data/mRNA_WS275_WITHregion_v3.csv'))
    all_mRNA = all_mRNA[['Gene ID','sequence_length','Gene name','CDS start','CDS end']]
    all_mRNA['len5'] = all_mRNA['CDS start'] - 1
    all_mRNA['lencds'] = all_mRNA['CDS end'] - all_mRNA['CDS start'] + 1
    all_mRNA['len3'] = all_mRNA['sequence_length'] - all_mRNA['CDS end']
    output = []

    for ANALYZE_TYPE in data['ANALYZE_TYPE']:
        #print('===={}===='.format(ANALYZE_TYPE)) 
        if ANALYZE_TYPE == 'TARGET_SCORE':
            for FILE_NAME in data['TARGET_SCORE_FILE']:
                for rc in data['READ_COUNT_TYPE']:
                    # TARGET SCORE SITE
                    if data['DATA'] == 1:
                        df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/target_score_'+FILE_NAME+'_'+rc+'norm_2step.csv'))
                        df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc','xgu_inseed', 'gu_inseed', 'xgu_innon-seed',
                                    'gu_innon-seed', 'totalmismatch', 'xGU_mispos', 'GU_mispos']]
                        df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                        all_mRNA_result = tool1_main(df, all_mRNA)
                        all_mRNA_result.to_csv(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'),index=None)
                        output.append(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))

                    # TARGET SCORE SITE CROSS CONDITION 
                    if data['DATA'] == 2:
                        df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/target_score_'+FILE_NAME+'_'+rc+'norm_usetRNA_2step.csv'))
                        df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc','xgu_inseed', 'gu_inseed', 'xgu_innon-seed',
                                    'gu_innon-seed', 'totalmismatch', 'xGU_mispos', 'GU_mispos']]
                        df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                        all_mRNA_result = tool1_main(df, all_mRNA)
                        all_mRNA_result.to_csv(os.path.abspath(__file__+'/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'),index=None)
                        output.append(os.path.abspath(__file__+'/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))

        # pirscan
        elif ANALYZE_TYPE == 'PIRSCAN':
            for FILE_NAME in data['PIR_FILE']:
                for rc in data['READ_COUNT_TYPE']:
                    #print(FILE_NAME, rc)
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/'+FILE_NAME+'_pirscan_withup.csv'))
                    df = df.rename(columns={'rem_tran_target_pos':'target_score_pos','score':'targeting_score'})
                    df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc']]
                    df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                    all_mRNA_result = tool1_main(df, all_mRNA)
                    all_mRNA_result.to_csv(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'),index=None)
                    output.append(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))

        #RNAUP
        elif ANALYZE_TYPE == 'RNAUP':
            for FILE_NAME in data['RNAUP_FILE']:
                for rc in data['READ_COUNT_TYPE']: 
                    #print(FILE_NAME, rc)
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/RNAup_'+FILE_NAME+'_'+rc+'norm_2step.csv'))
                    df = df[['RNAup_pos','RNAup_score','regulator_name','transcript_name','norm_rc']]
                    #print(len(df))
                    df['site_len'] = df['RNAup_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                    all_mRNA_result = tool1_FOR_RNAUP(df, all_mRNA)
                    all_mRNA_result.to_csv(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'),index=None)
                    output.append(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))

        # G22
        elif ANALYZE_TYPE == 'G22':
            for FILE_NAME in data['G22_FILE']:
                if data['DATA'] == 1:
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/22G/'+FILE_NAME+'_'+'norm_1step.csv'))
                    df = df.rename(columns={'evenly_rc':'norm_rc','rem_tran_target_pos':'target_score_pos','ref_id':'transcript_name','input_id':'regulator_name'})
                    df = df[['norm_rc','target_score_pos','regulator_name','transcript_name']]
                    df['site_len'] = df['target_score_pos'].apply(lambda x: int(1+abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                    all_mRNA_result = tool1_main_G22(df, all_mRNA)
                    all_mRNA_result.to_csv(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'),index=None)
                    output.append(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))

                # G22 CROSS CONDITION
                if data['DATA'] == 2:
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/22G/'+FILE_NAME+'_'+'EGL17M0_usemiRNA_norm_1step.csv')) #### piRNA -> 21
                    df = df.rename(columns={'evenly_rc':'norm_rc','rem_tran_target_pos':'target_score_pos','ref_id':'transcript_name','input_id':'regulator_name'})
                    df = df[['norm_rc','target_score_pos','regulator_name','transcript_name']]
                    df['site_len'] = df['target_score_pos'].apply(lambda x: int(1+abs((int(x.split('-')[0])-int(x.split('-')[1])))))
                    all_mRNA_result = tool1_main_G22(df, all_mRNA)
                    all_mRNA_result.to_csv(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'),index=None)
                    output.append(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool1.csv'))
    return output


def tool1_main(data, all_mRNA):
    all_utr5_count = []
    all_cds_count = []
    all_utr3_count = []
    m_count = 0
    for name, cds_1, cds_2 in zip(all_mRNA['Gene name'], all_mRNA['CDS start'], all_mRNA['CDS end']):
        print(m_count, end='\r')
        m_count += 1
        tmp = data[data['transcript_name'].isin([name])].reset_index(drop=True)
        UTR5 = [] #FOR THIS TRANSCRIPTS
        CDS = []
        UTR3 = []
        for i,row in tmp.iterrows():
            init_pos = int( row['target_score_pos'].split('-')[0].split('.')[0] )
            end_pos = int( row['target_score_pos'].split('-')[1].split('.')[0] )
            read_count = row['norm_rc']

            if end_pos < cds_1 :
                #score_5 += [score] * math.ceil(read_count)
                UTR5.append(1*read_count)
                CDS.append(0)
                UTR3.append(0)
            elif init_pos < cds_1 <= end_pos:
                #score_5cds += [score] * math.ceil(read_count)
                UTR5.append( (cds_1 - init_pos)/row['site_len'] *read_count )
                CDS.append( (end_pos - cds_1 + 1)/row['site_len'] *read_count)
                UTR3.append(0)
            elif init_pos >= cds_1 and end_pos <= cds_2:
                #score_cds += [score] * math.ceil(read_count)
                UTR5.append(0)
                CDS.append(1*read_count)
                UTR3.append(0)
            elif init_pos <= cds_2 < end_pos:
                #score_3cds += [score] * math.ceil(read_count)
                UTR5.append(0)
                CDS.append((cds_2 - init_pos + 1)/row['site_len'] *read_count)
                UTR3.append((end_pos - cds_2)/row['site_len'] *read_count)
            elif cds_2 < init_pos:
                #score_3 += [score] * math.ceil(read_count)
                UTR5.append(0)
                CDS.append(0)
                UTR3.append(1*read_count)
            else:
                UTR5.append('N')
                CDS.append('N')
                UTR3.append('N')
        all_utr5_count.append(sum(UTR5))
        all_cds_count.append(sum(CDS))
        all_utr3_count.append(sum(UTR3))
    all_mRNA['count5'] = all_utr5_count
    all_mRNA['countcds'] = all_cds_count
    all_mRNA['count3'] = all_utr3_count
    return all_mRNA

def tool1_FOR_RNAUP(data, all_mRNA):
    all_utr5_count = []
    all_cds_count = []
    all_utr3_count = []
    m_count = 0
    for name, cds_1, cds_2 in zip(all_mRNA['Gene name'], all_mRNA['CDS start'], all_mRNA['CDS end']):
        print(m_count, end='\r')
        m_count += 1
        tmp = data[data['transcript_name'].isin([name])].reset_index(drop=True)
        UTR5 = [] #FOR THIS TRANSCRIPTS
        CDS = []
        UTR3 = []
        for i,row in tmp.iterrows():
            init_pos = int( row['RNAup_pos'].split('-')[0].split('.')[0] )
            end_pos = int( row['RNAup_pos'].split('-')[1].split('.')[0] )
            read_count = row['norm_rc']

            if end_pos < cds_1 :
                #score_5 += [score] * math.ceil(read_count)
                UTR5.append(1*read_count)
                CDS.append(0)
                UTR3.append(0)
            elif init_pos < cds_1 <= end_pos:
                #score_5cds += [score] * math.ceil(read_count)
                UTR5.append( (cds_1 - init_pos)/row['site_len'] *read_count )
                CDS.append( (end_pos - cds_1 + 1)/row['site_len'] *read_count)
                UTR3.append(0)
            elif init_pos >= cds_1 and end_pos <= cds_2:
                #score_cds += [score] * math.ceil(read_count)
                UTR5.append(0)
                CDS.append(1*read_count)
                UTR3.append(0)
            elif init_pos <= cds_2 < end_pos:
                #score_3cds += [score] * math.ceil(read_count)
                UTR5.append(0)
                CDS.append((cds_2 - init_pos + 1)/row['site_len'] *read_count)
                UTR3.append((end_pos - cds_2)/row['site_len'] *read_count)
            elif cds_2 < init_pos:
                #score_3 += [score] * math.ceil(read_count)
                UTR5.append(0)
                CDS.append(0)
                UTR3.append(1*read_count)
            else:
                UTR5.append('N')
                CDS.append('N')
                UTR3.append('N')
        all_utr5_count.append(sum(UTR5))
        all_cds_count.append(sum(CDS))
        all_utr3_count.append(sum(UTR3))
    all_mRNA['count5'] = all_utr5_count
    all_mRNA['countcds'] = all_cds_count
    all_mRNA['count3'] = all_utr3_count
    return all_mRNA

def tool1_main_G22(data, all_mRNA):
    all_utr5_count = []
    all_cds_count = []
    all_utr3_count = []
    m_count = 0
    for name, cds_1, cds_2 in zip(all_mRNA['Gene name'], all_mRNA['CDS start'], all_mRNA['CDS end']):
        print(m_count, end='\r')
        m_count += 1
        tmp = data[data['transcript_name'].isin([name])].reset_index(drop=True)
        UTR5 = [] #FOR THIS TRANSCRIPTS
        CDS = []
        UTR3 = []
        for i,row in tmp.iterrows():
            init_pos = int( row['target_score_pos'].split('-')[0].split('.')[0] )
            end_pos = int( row['target_score_pos'].split('-')[1].split('.')[0] )
            read_count = row['norm_rc']

            if end_pos < cds_1 :
                #score_5 += [score] * math.ceil(read_count)
                UTR5.append(1*read_count)
                CDS.append(0)
                UTR3.append(0)
            elif init_pos < cds_1 <= end_pos:
                #score_5cds += [score] * math.ceil(read_count)
                UTR5.append( (cds_1 - init_pos)/row['site_len'] *read_count )
                CDS.append( (end_pos - cds_1 + 1)/row['site_len'] *read_count)
                UTR3.append(0)
            elif init_pos >= cds_1 and end_pos <= cds_2:
                #score_cds += [score] * math.ceil(read_count)
                UTR5.append(0)
                CDS.append(1*read_count)
                UTR3.append(0)
            elif init_pos <= cds_2 < end_pos:
                #score_3cds += [score] * math.ceil(read_count)
                UTR5.append(0)
                CDS.append((cds_2 - init_pos + 1)/row['site_len'] *read_count)
                UTR3.append((end_pos - cds_2)/row['site_len'] *read_count)
            elif cds_2 < init_pos:
                #score_3 += [score] * math.ceil(read_count)
                UTR5.append(0)
                CDS.append(0)
                UTR3.append(1*read_count)
            else:
                UTR5.append('N')
                CDS.append('N')
                UTR3.append('N')
        all_utr5_count.append(sum(UTR5))
        all_cds_count.append(sum(CDS))
        all_utr3_count.append(sum(UTR3))
    all_mRNA['count5'] = all_utr5_count
    all_mRNA['countcds'] = all_cds_count
    all_mRNA['count3'] = all_utr3_count
    return all_mRNA

def tool1_main_G22_per_site(data, all_mRNA):
    all_utr5_count = []
    all_cds_count = []
    all_utr3_count = []
    cross_utr5_cds = 0
    cross_cds_utr3 = 0
    m_count = 0
    for name, cds_1, cds_2 in zip(all_mRNA['Gene name'], all_mRNA['CDS start'], all_mRNA['CDS end']):
        print(m_count, end='\r')
        m_count += 1
        tmp = data[data['transcript_name'].isin([name])].reset_index(drop=True)
        UTR5 = [] #FOR THIS TRANSCRIPTS
        CDS = []
        UTR3 = []
        for i,row in tmp.iterrows():
            init_pos = int( row['site_pos'].split('-')[0].split('.')[0] )
            end_pos = int( row['site_pos'].split('-')[1].split('.')[0] )
            read_count = row['norm_rc']

            if end_pos < cds_1 :
                UTR5.append(1*read_count)
                CDS.append(0)
                UTR3.append(0)
            elif init_pos < cds_1 <= end_pos:
                UTR5.append(0)
                CDS.append(0)
                UTR3.append(0)
                cross_utr5_cds += 1
            elif init_pos >= cds_1 and end_pos <= cds_2:
                UTR5.append(0)
                CDS.append(1*read_count)
                UTR3.append(0)
            elif init_pos <= cds_2 < end_pos:
                UTR5.append(0)
                CDS.append(0)
                UTR3.append(0)
                cross_cds_utr3 += 1
            elif cds_2 < init_pos:
                UTR5.append(0)
                CDS.append(0)
                UTR3.append(1*read_count)
            else:
                UTR5.append('N')
                CDS.append('N')
                UTR3.append('N')
        all_utr5_count.append(sum(UTR5))
        all_cds_count.append(sum(CDS))
        all_utr3_count.append(sum(UTR3))
    all_mRNA['count5'] = all_utr5_count
    all_mRNA['countcds'] = all_cds_count
    all_mRNA['count3'] = all_utr3_count
    #print('cross_utr5_cds: {}\ncross_cds_utr3: {}'.format(cross_utr5_cds, cross_cds_utr3))
    return all_mRNA


