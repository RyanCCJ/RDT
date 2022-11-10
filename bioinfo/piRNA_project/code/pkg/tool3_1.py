import json
import os
import numpy as np
import pandas as pd

# main
def tool3_1(data):
    all_mRNA = pd.read_csv(os.path.abspath(__file__+'/../../../../data/mRNA_WS275_WITHregion_v3.csv'))
    all_mRNA = all_mRNA[['Gene ID','sequence_length','Gene name','CDS start','CDS end']]
    all_mRNA['len5'] = all_mRNA['CDS start'] - 1
    all_mRNA['lencds'] = all_mRNA['CDS end'] - all_mRNA['CDS start'] + 1
    all_mRNA['len3'] = all_mRNA['sequence_length'] - all_mRNA['CDS end']
    output = []

    for ANALYZE_TYPE in data['ANALYZE_TYPE']:
        #print('===={}===='.format(ANALYZE_TYPE))
        # target
        if ANALYZE_TYPE == 'TARGET_SCORE':
            if data['DATA'] == 1:
                for FILE_NAME in data['TARGET_SCORE_FILE']:
                    for rc in data['READ_COUNT_TYPE']:
                        df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/target_score_'+FILE_NAME+'_'+rc+'norm_2step.csv'))
                        df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc']]
                        #df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc','xgu_inseed', 'gu_inseed', 'xgu_innon-seed',
                        #    'gu_innon-seed', 'totalmismatch', 'xGU_mispos', 'GU_mispos']]
                        df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))

                        all_pos_START_record = {}
                        all_pos_END_record = {}
                        m_count = 0
                        for name,seq_len,cds_start,cds_end in zip(all_mRNA['Gene name'], all_mRNA['sequence_length'], all_mRNA['CDS start'], all_mRNA['CDS end']):
                            print(m_count, end='\r')
                            m_count += 1
                            tmp = df[df['transcript_name'].isin([name])].reset_index(drop=True)
                            this_trans_START_pos = {}
                            this_trans_END_pos = {}
                            for i, row in tmp.iterrows():
                                read_count = row['norm_rc']
                                start = int(row['target_score_pos'].split('-')[0]) #rnaup coordinate need plus one
                                end = row['target_score_pos'].split('-')[1]
                                for p in range(int(start), int(end)):
                                    relate_start = p - cds_start
                                    if relate_start not in this_trans_START_pos:
                                        this_trans_START_pos[relate_start] = read_count
                                    else:
                                        this_trans_START_pos[relate_start] += read_count
                                    this_trans_START_pos[relate_start] = round(this_trans_START_pos[relate_start],3)

                                    relate_end = p - cds_end
                                    if relate_end not in this_trans_END_pos:
                                        this_trans_END_pos[relate_end] = read_count
                                    else:
                                        this_trans_END_pos[relate_end] += read_count
                                    this_trans_END_pos[relate_end] = round(this_trans_END_pos[relate_end],3)
                            if 1-cds_start not in this_trans_START_pos:
                                this_trans_START_pos[1-cds_start] = 0
                            if int(seq_len) - cds_start not in this_trans_START_pos:
                                this_trans_START_pos[int(seq_len) - cds_start] = 0

                            if 1-cds_end not in this_trans_END_pos:
                                this_trans_END_pos[1-cds_end] = 0
                            if int(seq_len) -cds_end not in this_trans_END_pos:
                                this_trans_END_pos[int(seq_len) - cds_end] = 0    

                            this_trans_START_pos = dict(sorted(this_trans_START_pos.items(), key=lambda e:e[0]))
                            this_trans_END_pos = dict(sorted(this_trans_END_pos.items(), key=lambda e:e[0]))
                            all_pos_START_record[name] = this_trans_START_pos
                            all_pos_END_record[name] = this_trans_END_pos

                        with open(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_START.txt'),'w') as f:
                            f.write(json.dumps(all_pos_START_record))
                        with open(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_STOP.txt'),'w') as f:
                            f.write(json.dumps(all_pos_END_record))

                for rc in data['READ_COUNT_TYPE']:
                    for FILE_NAME in data['TARGET_SCORE_FILE']:
                        for EXTEND in data['EXTEND_RANGE']:
                            for SIDE in ['START','STOP']:
                                #print(rc, FILE_NAME, EXTEND, SIDE)
                                make_length_file_target(FILE_NAME, rc, EXTEND, SIDE, output)

            # cross condition
            if data['DATA'] == 2:
                for FILE_NAME in data['TARGET_SCORE_FILE']:
                    for rc in data['READ_COUNT_TYPE']:
                        df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/target_score_'+FILE_NAME+'_'+rc+'norm_usetRNA_2step.csv'))
                        df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc']]
                        #df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc','xgu_inseed', 'gu_inseed', 'xgu_innon-seed',
                        #    'gu_innon-seed', 'totalmismatch', 'xGU_mispos', 'GU_mispos']]
                        df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))

                        all_pos_START_record = {}
                        all_pos_END_record = {}
                        m_count = 0
                        for name,seq_len,cds_start,cds_end in zip(all_mRNA['Gene name'], all_mRNA['sequence_length'], all_mRNA['CDS start'], all_mRNA['CDS end']):
                            print(m_count, end='\r')
                            m_count += 1
                            tmp = df[df['transcript_name'].isin([name])].reset_index(drop=True)
                            this_trans_START_pos = {}
                            this_trans_END_pos = {}
                            for i, row in tmp.iterrows():
                                read_count = row['norm_rc']
                                start = int(row['target_score_pos'].split('-')[0]) #rnaup coordinate need plus one
                                end = row['target_score_pos'].split('-')[1]
                                for p in range(int(start), int(end)):
                                    relate_start = p - cds_start
                                    if relate_start not in this_trans_START_pos:
                                        this_trans_START_pos[relate_start] = read_count
                                    else:
                                        this_trans_START_pos[relate_start] += read_count
                                    this_trans_START_pos[relate_start] = round(this_trans_START_pos[relate_start],3)

                                    relate_end = p - cds_end
                                    if relate_end not in this_trans_END_pos:
                                        this_trans_END_pos[relate_end] = read_count
                                    else:
                                        this_trans_END_pos[relate_end] += read_count
                                    this_trans_END_pos[relate_end] = round(this_trans_END_pos[relate_end],3)
                            if 1-cds_start not in this_trans_START_pos:
                                this_trans_START_pos[1-cds_start] = 0
                            if int(seq_len) - cds_start not in this_trans_START_pos:
                                this_trans_START_pos[int(seq_len) - cds_start] = 0

                            if 1-cds_end not in this_trans_END_pos:
                                this_trans_END_pos[1-cds_end] = 0
                            if int(seq_len) -cds_end not in this_trans_END_pos:
                                this_trans_END_pos[int(seq_len) - cds_end] = 0    

                            this_trans_START_pos = dict(sorted(this_trans_START_pos.items(), key=lambda e:e[0]))
                            this_trans_END_pos = dict(sorted(this_trans_END_pos.items(), key=lambda e:e[0]))
                            all_pos_START_record[name] = this_trans_START_pos
                            all_pos_END_record[name] = this_trans_END_pos

                        with open(os.path.abspath(__file__+'/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_START.txt'),'w') as f:
                            f.write(json.dumps(all_pos_START_record))
                        with open(os.path.abspath(__file__+'/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_STOP.txt'),'w') as f:
                            f.write(json.dumps(all_pos_END_record))

                for rc in data['READ_COUNT_TYPE']:
                    for FILE_NAME in data['TARGET_SCORE_FILE']:
                        for EXTEND in data['EXTEND_RANGE']:
                            for SIDE in ['START','STOP']:
                                #print(rc, FILE_NAME, EXTEND, SIDE)
                                make_length_file_nor(FILE_NAME, rc, EXTEND, SIDE, output)

        # pirscan
        elif ANALYZE_TYPE == 'PIRSCAN':
            for FILE_NAME in data['PIR_FILE']:
                for rc in data['READ_COUNT_TYPE']:
                    #print(FILE_NAME, rc)
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/'+FILE_NAME+'_pirscan_withup.csv'))
                    df = df.rename(columns={'rem_tran_target_pos':'target_score_pos','RNAup_score':'targeting_score'})

                    df = df[['target_score_pos','targeting_score','regulator_name','transcript_name','norm_rc']]
                    df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))

                    all_pos_START_record = {}
                    all_pos_END_record = {}
                    m_count = 0
                    for name,seq_len,cds_start,cds_end in zip(all_mRNA['Gene name'], all_mRNA['sequence_length'], all_mRNA['CDS start'], all_mRNA['CDS end']):
                        print(m_count, end='\r')
                        m_count += 1
                        tmp = df[df['transcript_name'].isin([name])].reset_index(drop=True)
                        this_trans_START_pos = {}
                        this_trans_END_pos = {}
                        for i, row in tmp.iterrows():
                            read_count = row['norm_rc']
                            start = int(row['target_score_pos'].split('-')[0]) #rnaup coordinate need plus one
                            end = row['target_score_pos'].split('-')[1]
                            for p in range(int(start), int(end)):
                                relate_start = p - cds_start
                                if relate_start not in this_trans_START_pos:
                                    this_trans_START_pos[relate_start] = read_count
                                else:
                                    this_trans_START_pos[relate_start] += read_count
                                this_trans_START_pos[relate_start] = round(this_trans_START_pos[relate_start],3)

                                relate_end = p - cds_end
                                if relate_end not in this_trans_END_pos:
                                    this_trans_END_pos[relate_end] = read_count
                                else:
                                    this_trans_END_pos[relate_end] += read_count
                                this_trans_END_pos[relate_end] = round(this_trans_END_pos[relate_end],3)
                        if 1-cds_start not in this_trans_START_pos:
                            this_trans_START_pos[1-cds_start] = 0
                        if int(seq_len) - cds_start not in this_trans_START_pos:
                            this_trans_START_pos[int(seq_len) - cds_start] = 0

                        if 1-cds_end not in this_trans_END_pos:
                            this_trans_END_pos[1-cds_end] = 0
                        if int(seq_len) -cds_end not in this_trans_END_pos:
                            this_trans_END_pos[int(seq_len) - cds_end] = 0    

                        this_trans_START_pos = dict(sorted(this_trans_START_pos.items(), key=lambda e:e[0]))
                        this_trans_END_pos = dict(sorted(this_trans_END_pos.items(), key=lambda e:e[0]))
                        all_pos_START_record[name] = this_trans_START_pos
                        all_pos_END_record[name] = this_trans_END_pos

                    with open(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_START.txt'),'w') as f:
                        f.write(json.dumps(all_pos_START_record))
                    with open(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_STOP.txt'),'w') as f:
                        f.write(json.dumps(all_pos_END_record))

            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['PIR_FILE']:
                    for EXTEND in data['EXTEND_RANGE']:
                        for SIDE in ['START','STOP']:
                            #print(rc, FILE_NAME, EXTEND, SIDE)
                            make_length_file_pir(FILE_NAME, rc, EXTEND, SIDE, output)

        #rnaup
        elif ANALYZE_TYPE == 'RNAUP':
            for FILE_NAME in data['RNAUP_FILE']:
                for rc in data['READ_COUNT_TYPE']:
                    #print(FILE_NAME, rc)
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/CLASH/RNAup_'+FILE_NAME+'_'+rc+'norm_2step.csv'))
                    df = df[['RNAup_pos','RNAup_score','regulator_name','transcript_name','norm_rc']]
                    df['site_len'] = df['RNAup_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))

                    all_pos_START_record = {}
                    all_pos_END_record = {}
                    m_count = 0
                    for name,seq_len,cds_start,cds_end in zip(all_mRNA['Gene name'], all_mRNA['sequence_length'], all_mRNA['CDS start'], all_mRNA['CDS end']):
                        print(m_count, end='\r')
                        m_count += 1
                        tmp = df[df['transcript_name'].isin([name])].reset_index(drop=True)
                        this_trans_START_pos = {}
                        this_trans_END_pos = {}
                        for i, row in tmp.iterrows():
                            read_count = row['norm_rc']
                            start = int(row['RNAup_pos'].split('-')[0]) #rnaup coordinate need plus one
                            end = row['RNAup_pos'].split('-')[1]
                            for p in range(int(start), int(end)):
                                relate_start = p - cds_start
                                if relate_start not in this_trans_START_pos:
                                    this_trans_START_pos[relate_start] = read_count
                                else:
                                    this_trans_START_pos[relate_start] += read_count
                                this_trans_START_pos[relate_start] = round(this_trans_START_pos[relate_start],3)

                                relate_end = p - cds_end
                                if relate_end not in this_trans_END_pos:
                                    this_trans_END_pos[relate_end] = read_count
                                else:
                                    this_trans_END_pos[relate_end] += read_count
                                this_trans_END_pos[relate_end] = round(this_trans_END_pos[relate_end],3)
                        if 1-cds_start not in this_trans_START_pos:
                            this_trans_START_pos[1-cds_start] = 0
                        if int(seq_len) - cds_start not in this_trans_START_pos:
                            this_trans_START_pos[int(seq_len) - cds_start] = 0

                        if 1-cds_end not in this_trans_END_pos:
                            this_trans_END_pos[1-cds_end] = 0
                        if int(seq_len) -cds_end not in this_trans_END_pos:
                            this_trans_END_pos[int(seq_len) - cds_end] = 0    

                        this_trans_START_pos = dict(sorted(this_trans_START_pos.items(), key=lambda e:e[0]))
                        this_trans_END_pos = dict(sorted(this_trans_END_pos.items(), key=lambda e:e[0]))
                        all_pos_START_record[name] = this_trans_START_pos
                        all_pos_END_record[name] = this_trans_END_pos

                    with open(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_START.txt'),'w') as f:
                        f.write(json.dumps(all_pos_START_record))
                    with open(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_STOP.txt'),'w') as f:
                        f.write(json.dumps(all_pos_END_record))

            for rc in data['READ_COUNT_TYPE']:
                for FILE_NAME in data['RNAUP_FILE']:
                    for EXTEND in data['EXTEND_RANGE']:
                        for SIDE in ['START','STOP']:
                            #print(rc, FILE_NAME, EXTEND, SIDE)
                            make_length_file_rnaup(FILE_NAME, rc, EXTEND, SIDE, output)

        # G22
        elif ANALYZE_TYPE == 'G22':
            if data['DATA'] == 1:
                for FILE_NAME in data['G22_FILE']:    
                    #print(FILE_NAME)
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/22G/'+FILE_NAME+'_'+'norm_1step.csv'))
                    df = df.rename(columns={'evenly_rc':'norm_rc','rem_tran_target_pos':'target_score_pos','ref_id':'transcript_name','input_id':'regulator_name'})
                    df = df[['norm_rc','target_score_pos','regulator_name','transcript_name']]
                    df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))

                    all_pos_START_record = {}
                    all_pos_END_record = {}
                    
                    all_pos_HEAD_record = {}
                    all_pos_TAIL_record = {}
                    m_count = 0
                    for name,seq_len,cds_start,cds_end in zip(all_mRNA['Gene name'], all_mRNA['sequence_length'], all_mRNA['CDS start'], all_mRNA['CDS end']):
                        print(m_count, end='\r')
                        m_count += 1
                        tmp = df[df['transcript_name'].isin([name])].reset_index(drop=True)
                        this_trans_START_pos = {}
                        this_trans_END_pos = {}
                        
                        this_trans_HEAD_pos = {}
                        this_trans_TAIL_pos = {}
                        for i, row in tmp.iterrows():
                            read_count = row['norm_rc']
                            start = int(row['target_score_pos'].split('-')[0]) #rnaup coordinate need plus one
                            end = row['target_score_pos'].split('-')[1]
                            for p in range(int(start), int(end)):
                                relate_start = p - cds_start
                                relate_head = p
                                if relate_start not in this_trans_START_pos:
                                    this_trans_START_pos[relate_start] = read_count
                                else:
                                    this_trans_START_pos[relate_start] += read_count
                                this_trans_START_pos[relate_start] = round(this_trans_START_pos[relate_start],3)
                                if relate_head not in this_trans_HEAD_pos:
                                    this_trans_HEAD_pos[relate_head] = read_count
                                else:
                                    this_trans_HEAD_pos[relate_head] += read_count
                                this_trans_HEAD_pos[relate_head] = round(this_trans_HEAD_pos[relate_head],3)
                                
                                relate_end = p - cds_end
                                relate_tail = seq_len - p
                                if relate_end not in this_trans_END_pos:
                                    this_trans_END_pos[relate_end] = read_count
                                else:
                                    this_trans_END_pos[relate_end] += read_count
                                this_trans_END_pos[relate_end] = round(this_trans_END_pos[relate_end],3)
                                if relate_tail not in this_trans_TAIL_pos:
                                    this_trans_TAIL_pos[relate_tail] = read_count
                                else:
                                    this_trans_TAIL_pos[relate_tail] += read_count
                                this_trans_TAIL_pos[relate_tail] = round(this_trans_TAIL_pos[relate_tail],3)
                                
                        if 1-cds_start not in this_trans_START_pos:
                            this_trans_START_pos[1-cds_start] = 0
                        if int(seq_len) - cds_start not in this_trans_START_pos:
                            this_trans_START_pos[int(seq_len) - cds_start] = 0
                        if 1-cds_end not in this_trans_END_pos:
                            this_trans_END_pos[1-cds_end] = 0
                        if int(seq_len) -cds_end not in this_trans_END_pos:
                            this_trans_END_pos[int(seq_len) - cds_end] = 0
                            
                        if 1 not in this_trans_HEAD_pos:
                            this_trans_HEAD_pos[1] = 0
                        if int(seq_len) not in this_trans_HEAD_pos:
                            this_trans_HEAD_pos[int(seq_len)] = 0
                        if 1 not in this_trans_TAIL_pos:
                            this_trans_TAIL_pos[1] = 0
                        if int(seq_len) not in this_trans_TAIL_pos:
                            this_trans_TAIL_pos[int(seq_len)] = 0
                            
                        this_trans_START_pos = dict(sorted(this_trans_START_pos.items(), key=lambda e:e[0]))
                        this_trans_END_pos = dict(sorted(this_trans_END_pos.items(), key=lambda e:e[0]))
                        all_pos_START_record[name] = this_trans_START_pos
                        all_pos_END_record[name] = this_trans_END_pos
                        
                        this_trans_HEAD_pos = dict(sorted(this_trans_HEAD_pos.items(), key=lambda e:e[0]))
                        this_trans_TAIL_pos = dict(sorted(this_trans_TAIL_pos.items(), key=lambda e:e[0]))
                        all_pos_HEAD_record[name] = this_trans_HEAD_pos
                        all_pos_TAIL_record[name] = this_trans_TAIL_pos
                        
                    with open(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_START.txt'),'w') as f:
                        f.write(json.dumps(all_pos_START_record))
                    with open(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_STOP.txt'),'w') as f:
                        f.write(json.dumps(all_pos_END_record))
                    with open(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_HEAD.txt'),'w') as f:
                        f.write(json.dumps(all_pos_HEAD_record))
                    with open(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_TAIL.txt'),'w') as f:
                        f.write(json.dumps(all_pos_TAIL_record))
                        
                for FILE_NAME in data['G22_FILE']:    
                    for EXTEND in data['EXTEND_RANGE']:
                        for SIDE in ['START','STOP']:
                            #print(FILE_NAME, EXTEND, SIDE)
                            make_length_file_g22(FILE_NAME, EXTEND, SIDE, output)
                        for SIDE in ['HEAD','TAIL']:
                            #print(FILE_NAME, EXTEND, SIDE)
                            make_length_file_HT_g22(FILE_NAME, EXTEND, SIDE, output)

            # NOR_G22
            if data['DATA'] == 2:
                for FILE_NAME in data['G22_FILE']:    
                    #print(FILE_NAME)
                    df = pd.read_csv(os.path.abspath(__file__+'/../../../../data/22G/'+FILE_NAME+'_'+'EGL17M0_usemiRNA_norm_1step.csv')) #### piRNA -> 21
                    df = df.rename(columns={'evenly_rc':'norm_rc','rem_tran_target_pos':'target_score_pos','ref_id':'transcript_name','input_id':'regulator_name'})
                    df = df[['norm_rc','target_score_pos','regulator_name','transcript_name']]
                    df['site_len'] = df['target_score_pos'].apply(lambda x: int(abs((int(x.split('-')[0])-int(x.split('-')[1])))))

                    all_pos_START_record = {}
                    all_pos_END_record = {}
                    
                    all_pos_HEAD_record = {}
                    all_pos_TAIL_record = {}
                    m_count = 0
                    for name,seq_len,cds_start,cds_end in zip(all_mRNA['Gene name'], all_mRNA['sequence_length'], all_mRNA['CDS start'], all_mRNA['CDS end']):
                        print(m_count, end='\r')
                        m_count += 1
                        tmp = df[df['transcript_name'].isin([name])].reset_index(drop=True)
                        this_trans_START_pos = {}
                        this_trans_END_pos = {}
                        
                        this_trans_HEAD_pos = {}
                        this_trans_TAIL_pos = {}
                        
                        for i, row in tmp.iterrows():
                            read_count = row['norm_rc']
                            start = int(row['target_score_pos'].split('-')[0]) #rnaup coordinate need plus one
                            end = row['target_score_pos'].split('-')[1]
                            for p in range(int(start), int(end)):
                                
                                # START HEAD
                                relate_start = p - cds_start
                                relate_head = p
                                if relate_start not in this_trans_START_pos:
                                    this_trans_START_pos[relate_start] = read_count
                                else:
                                    this_trans_START_pos[relate_start] += read_count
                                this_trans_START_pos[relate_start] = round(this_trans_START_pos[relate_start],3)
                                
                                if relate_head not in this_trans_HEAD_pos:
                                    this_trans_HEAD_pos[relate_head] = read_count
                                else:
                                    this_trans_HEAD_pos[relate_head] += read_count
                                this_trans_HEAD_pos[relate_head] = round(this_trans_HEAD_pos[relate_head],3)
                                
                                # END TAIL
                                relate_end = p - cds_end
                                relate_tail = seq_len - p
                                if relate_end not in this_trans_END_pos:
                                    this_trans_END_pos[relate_end] = read_count
                                else:
                                    this_trans_END_pos[relate_end] += read_count
                                this_trans_END_pos[relate_end] = round(this_trans_END_pos[relate_end],3)
                                
                                if relate_tail not in this_trans_TAIL_pos:
                                    this_trans_TAIL_pos[relate_tail] = read_count
                                else:
                                    this_trans_TAIL_pos[relate_tail] += read_count
                                this_trans_TAIL_pos[relate_tail] = round(this_trans_TAIL_pos[relate_tail],3)
                                
                        if 1-cds_start not in this_trans_START_pos:
                            this_trans_START_pos[1-cds_start] = 0
                        if int(seq_len) - cds_start not in this_trans_START_pos:
                            this_trans_START_pos[int(seq_len) - cds_start] = 0

                        if 1-cds_end not in this_trans_END_pos:
                            this_trans_END_pos[1-cds_end] = 0
                        if int(seq_len) -cds_end not in this_trans_END_pos:
                            this_trans_END_pos[int(seq_len) - cds_end] = 0    
                        
                        if 1 not in this_trans_HEAD_pos:
                            this_trans_HEAD_pos[1] = 0
                        if int(seq_len) not in this_trans_HEAD_pos:
                            this_trans_HEAD_pos[int(seq_len)] = 0

                        if 1 not in this_trans_TAIL_pos:
                            this_trans_TAIL_pos[1] = 0
                        if int(seq_len) not in this_trans_TAIL_pos:
                            this_trans_TAIL_pos[int(seq_len)] = 0
                            
                        this_trans_START_pos = dict(sorted(this_trans_START_pos.items(), key=lambda e:e[0]))
                        this_trans_END_pos = dict(sorted(this_trans_END_pos.items(), key=lambda e:e[0]))
                        all_pos_START_record[name] = this_trans_START_pos
                        all_pos_END_record[name] = this_trans_END_pos
                        
                        this_trans_HEAD_pos = dict(sorted(this_trans_HEAD_pos.items(), key=lambda e:e[0]))
                        this_trans_TAIL_pos = dict(sorted(this_trans_TAIL_pos.items(), key=lambda e:e[0]))
                        all_pos_HEAD_record[name] = this_trans_HEAD_pos
                        all_pos_TAIL_record[name] = this_trans_TAIL_pos
                    
                    with open(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_START.txt'),'w') as f:
                        f.write(json.dumps(all_pos_START_record))
                    with open(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_STOP.txt'),'w') as f:
                        f.write(json.dumps(all_pos_END_record))
                        
                    with open(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_HEAD.txt'),'w') as f:
                        f.write(json.dumps(all_pos_HEAD_record))
                    with open(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_TAIL.txt'),'w') as f:
                        f.write(json.dumps(all_pos_TAIL_record))

                for FILE_NAME in data['G22_FILE']:    
                    for EXTEND in data['EXTEND_RANGE']:
                        for SIDE in ['START','STOP']:
                            #print(FILE_NAME, EXTEND, SIDE)
                            make_length_file_nor_g22(FILE_NAME, EXTEND, SIDE, output) 
                        for SIDE in ['HEAD','TAIL']:
                            #print(FILE_NAME, EXTEND, SIDE)
                            make_length_file_HT_nor_g22(FILE_NAME, EXTEND, SIDE, output)
    return output


def make_length_file_target(FILE_NAME, rc, EXTEND, SIDE, output):
    with open(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'.txt'),'r') as f:
        data = json.load(f)
    
    all_name = []
    all_start = []

    for name, value in data.items():

        # single transcript
        all_name.append(name)
        start_list = [0]*(EXTEND*2+1)
        for k,v in value.items():
            relate_pos = int(k)+EXTEND
            if 0 <= relate_pos <= EXTEND*2:
                start_list[relate_pos] += v
            elif relate_pos < 0:
                continue
            elif relate_pos > EXTEND*2:
                break

        tmp = list(value.items())
        #check head boundary
        head = int(tmp[0][0])
        if head > (-1)*EXTEND:
            N_list = ['N'] * abs(head - (-1)*EXTEND)
            start_list = N_list + start_list[head+EXTEND:]
        #check tail boundary
        tail = int(tmp[-1][0])
        if tail < EXTEND:
            N_list = ['N'] * abs(tail - EXTEND)
            start_list =  start_list[:tail+EXTEND+1] + N_list 
        all_start.append(start_list)
    start_df = pd.DataFrame(all_start, columns =list(np.arange(-EXTEND,EXTEND+1)), index=all_name)
    start_df.to_csv(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    output.append(os.path.abspath(__file__+'/../../../target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    
def make_length_file_nor(FILE_NAME, rc, EXTEND, SIDE, output):
    with open(os.path.abspath(__file__+'/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'.txt'),'r') as f:
        data = json.load(f)
    
    all_name = []
    all_start = []

    for name, value in data.items():

        # single transcript
        all_name.append(name)
        start_list = [0]*(EXTEND*2+1)
        for k,v in value.items():
            relate_pos = int(k)+EXTEND
            if 0 <= relate_pos <= EXTEND*2:
                start_list[relate_pos] += v
            elif relate_pos < 0:
                continue
            elif relate_pos > EXTEND*2:
                break

        tmp = list(value.items())
        #check head boundary
        head = int(tmp[0][0])
        if head > (-1)*EXTEND:
            N_list = ['N'] * abs(head - (-1)*EXTEND)
            start_list = N_list + start_list[head+EXTEND:]
        #check tail boundary
        tail = int(tmp[-1][0])
        if tail < EXTEND:
            N_list = ['N'] * abs(tail - EXTEND)
            start_list =  start_list[:tail+EXTEND+1] + N_list 
        all_start.append(start_list)
    start_df = pd.DataFrame(all_start, columns =list(np.arange(-EXTEND,EXTEND+1)), index=all_name)
    start_df.to_csv(os.path.abspath(__file__+'/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    output.append(os.path.abspath(__file__+'/../../../usetRNA_target_score_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    
def make_length_file_pir(FILE_NAME, rc, EXTEND, SIDE, output):
    with open(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'.txt'),'r') as f:
        data = json.load(f)
    
    all_name = []
    all_start = []

    for name, value in data.items():

        # single transcript
        all_name.append(name)
        start_list = [0]*(EXTEND*2+1)
        for k,v in value.items():
            relate_pos = int(k)+EXTEND
            if 0 <= relate_pos <= EXTEND*2:
                start_list[relate_pos] += v
            elif relate_pos < 0:
                continue
            elif relate_pos > EXTEND*2:
                break

        tmp = list(value.items())
        #check head boundary
        head = int(tmp[0][0])
        if head > (-1)*EXTEND:
            N_list = ['N'] * abs(head - (-1)*EXTEND)
            start_list = N_list + start_list[head+EXTEND:]
        #check tail boundary
        tail = int(tmp[-1][0])
        if tail < EXTEND:
            N_list = ['N'] * abs(tail - EXTEND)
            start_list =  start_list[:tail+EXTEND+1] + N_list 
        all_start.append(start_list)
    start_df = pd.DataFrame(all_start, columns =list(np.arange(-EXTEND,EXTEND+1)), index=all_name)
    start_df.to_csv(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    output.append(os.path.abspath(__file__+'/../../../pirscan_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    
def make_length_file_rnaup(FILE_NAME, rc, EXTEND, SIDE, output):
    with open(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'.txt'),'r') as f:
        data = json.load(f)
    
    all_name = []
    all_start = []

    for name, value in data.items():

        # single transcript
        all_name.append(name)
        start_list = [0]*(EXTEND*2+1)
        for k,v in value.items():
            relate_pos = int(k)+EXTEND
            if 0 <= relate_pos <= EXTEND*2:
                start_list[relate_pos] += v
            elif relate_pos < 0:
                continue
            elif relate_pos > EXTEND*2:
                break

        tmp = list(value.items())
        #check head boundary
        head = int(tmp[0][0])
        if head > (-1)*EXTEND:
            N_list = ['N'] * abs(head - (-1)*EXTEND)
            start_list = N_list + start_list[head+EXTEND:]
        #check tail boundary
        tail = int(tmp[-1][0])
        if tail < EXTEND:
            N_list = ['N'] * abs(tail - EXTEND)
            start_list =  start_list[:tail+EXTEND+1] + N_list 
        all_start.append(start_list)
    start_df = pd.DataFrame(all_start, columns =list(np.arange(-EXTEND,EXTEND+1)), index=all_name)
    start_df.to_csv(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    output.append(os.path.abspath(__file__+'/../../../RNAup_'+rc+'_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    
def make_length_file_g22(FILE_NAME, EXTEND, SIDE, output):
    with open(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'.txt'),'r') as f:
        data = json.load(f)
    
    all_name = []
    all_start = []

    for name, value in data.items():

        # single transcript
        all_name.append(name)
        start_list = [0]*(EXTEND*2+1)
        for k,v in value.items():
            relate_pos = int(k)+EXTEND
            if 0 <= relate_pos <= EXTEND*2:
                start_list[relate_pos] += v
            elif relate_pos < 0:
                continue
            elif relate_pos > EXTEND*2:
                break

        tmp = list(value.items())
        #check head boundary
        head = int(tmp[0][0])
        if head > (-1)*EXTEND:
            N_list = ['N'] * abs(head - (-1)*EXTEND)
            start_list = N_list + start_list[head+EXTEND:]
        #check tail boundary
        tail = int(tmp[-1][0])
        if tail < EXTEND:
            N_list = ['N'] * abs(tail - EXTEND)
            start_list =  start_list[:tail+EXTEND+1] + N_list 
        all_start.append(start_list)
    start_df = pd.DataFrame(all_start, columns =list(np.arange(-EXTEND,EXTEND+1)), index=all_name)
    start_df.to_csv(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    output.append(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    
def make_length_file_HT_g22(FILE_NAME, EXTEND, SIDE, output):
    with open(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'.txt'),'r') as f:
        data = json.load(f)
    
    all_name = []
    all_start = []

    for name, value in data.items():

        # single transcript
        all_name.append(name)
        start_list = [0]*EXTEND
        for k,v in value.items():
            relate_pos = int(k)
            if 0 <= relate_pos < EXTEND: # 0~99
                start_list[relate_pos] += v
            elif relate_pos < 0:
                continue
            elif relate_pos >= EXTEND:
                break
            
        all_start.append(start_list)
        
    start_df = pd.DataFrame(all_start, columns =list(np.arange(EXTEND)), index=all_name)
    start_df.to_csv(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    output.append(os.path.abspath(__file__+'/../../../22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    
def make_length_file_nor_g22(FILE_NAME, EXTEND, SIDE, output):
    with open(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'.txt'),'r') as f:
        data = json.load(f)
    
    all_name = []
    all_start = []

    for name, value in data.items():

        # single transcript
        all_name.append(name)
        start_list = [0]*(EXTEND*2+1)
        for k,v in value.items():
            relate_pos = int(k)+EXTEND
            if 0 <= relate_pos <= EXTEND*2:
                start_list[relate_pos] += v
            elif relate_pos < 0:
                continue
            elif relate_pos > EXTEND*2:
                break

        tmp = list(value.items())
        #check head boundary
        head = int(tmp[0][0])
        if head > (-1)*EXTEND:
            N_list = ['N'] * abs(head - (-1)*EXTEND)
            start_list = N_list + start_list[head+EXTEND:]
        #check tail boundary
        tail = int(tmp[-1][0])
        if tail < EXTEND:
            N_list = ['N'] * abs(tail - EXTEND)
            start_list =  start_list[:tail+EXTEND+1] + N_list 
        all_start.append(start_list)
    start_df = pd.DataFrame(all_start, columns =list(np.arange(-EXTEND,EXTEND+1)), index=all_name)
    start_df.to_csv(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    #start_df.to_csv(os.path.abspath(__file__+'/../../../../monday/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'),index=None)
    output.append(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))

def make_length_file_HT_nor_g22(FILE_NAME, EXTEND, SIDE, output):
    with open(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'.txt'),'r') as f:
        data = json.load(f)
    
    all_name = []
    all_start = []

    for name, value in data.items():

        # single transcript
        all_name.append(name)
        start_list = [0]*EXTEND
        for k,v in value.items():
            relate_pos = int(k)
            if 0 <= relate_pos < EXTEND: # 0~99
                start_list[relate_pos] += v
            elif relate_pos < 0:
                continue
            elif relate_pos >= EXTEND:
                break
            
        all_start.append(start_list)
        
    start_df = pd.DataFrame(all_start, columns =list(np.arange(EXTEND)), index=all_name)
    start_df.to_csv(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))
    #start_df.to_csv(os.path.abspath(__file__+'/../../../../monday/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'),index=None)
    output.append(os.path.abspath(__file__+'/../../../usemiRNA_22G_output/'+FILE_NAME+'_all_mRNA_tool3_'+SIDE+'_'+str(EXTEND)+'.csv'))