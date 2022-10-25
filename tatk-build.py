import argparse
import os
import sys
import yaml
from bioinfo.piRNA_project.code import tatk

__version__ = "version 1.0"


def base_name(file_path):
    # exmaple: /path/to/some/file.txt -> file
    return os.path.splitext(os.path.basename(file_path))[0]


if __name__ == '__main__':

    # arguments
    parser = argparse.ArgumentParser(
        description="This program is to build a reference for Transcriptome Analysis Toolkit (TATK).",
    )
    parser.add_argument("-v", "--version",
                        action="version",
                        version="%(prog)s " + __version__)
    parser.add_argument("-o", "--output",
                        help="output configuration in YAML format")
    parser.add_argument("-m", "--mode",
                        choices=["density", "metagene", "codon", "fold_change"],
                        help="choose a tool to analyze dataset")
    parser.add_argument("-t", "--type",
                        choices=["target_score", "RNAup", "pirScan", "22G"],
                        help="choose a type of the tool")
    parser.add_argument("-D",
                        action="store_true",
                        help="compare two datasets")
    parser.add_argument("--data",
                        help="An comma-separated file in specific format.")
    parser.add_argument("--data2",
                        help="Another comma-separated file in specific format.")
    parser.add_argument("-l", "--level",
                        choices=["read_count", "site"],
                        help='choose read-count level or site level')
    args = vars(parser.parse_args())

    data = {}

    # check if some necessary arguments are exist
    nec_args = args.copy()
    for arg in ['D', 'data2']:
        nec_args.pop(arg)
    program_exit = False
    for arg in nec_args:    
        if nec_args[arg] == None:
            print('[Error] Argument "--' + arg + '" should not left blank.')
            program_exit = True
    if program_exit:
        sys.exit()


    # mode
    if args['mode'] == 'density':
        data['TOOL'] = [1]
    elif args['mode'] == 'metagene':
        data['TOOL'] = [2]
    elif args['mode'] == 'codon':
        data['TOOL'] = [3]
    elif args['mode'] == 'fold_change':
        data['TOOL'] = [4]
    else:
        print('[Error] Unknown mode "' + args['mode'] + '". Please choose one from "density", "metagene", "codon", "fold_change".')
        sys.exit()


    # type
    if args['type'] == 'target_score':
        data['ANALYZE_TYPE'] = ['TARGET_SCORE']
    elif args['type'] == 'RNAup':
        data['ANALYZE_TYPE'] = ['RNAUP']
    elif args['type'] == 'pirScan':
        data['ANALYZE_TYPE'] = ['PIRSCAN']
    elif args['type'] == '22G':
        data['ANALYZE_TYPE'] = ['G22']
    else:
        print('[Error] Unknown mode "' + args['type'] + '". Please choose one from "target-score", "RNAup", "pirScan", "22G".')
        sys.exit()
    

    # data
    if args['D']:
        data['DATA'] = 2
        if args['data2'] == None:
            print('[Error] Argument "--data2" should not left blank while "-D" turns on.')
            sys.exit()
        if (data['ANALYZE_TYPE'] == ['RNAUP']) or (data['ANALYZE_TYPE'] == ['PIRSCAN']):
            print('[Error] Compare two datasets (-D) only when type is "target_score" or "22G".')
            sys.exit()
        elif data['ANALYZE_TYPE'] == ['TARGET_SCORE']:
            data['TARGET_SCORE_FILE'] = [base_name(args['data']), base_name(args['data2'])]
        elif data['ANALYZE_TYPE'] == ['G22']:
            data['G22_FILE'] = [base_name(args['data']), base_name(args['data2'])]
        data['COMPARE_DATA1'] = [base_name(args['data'])]
        data['COMPARE_DATA2'] = [base_name(args['data2'])]
        data['title_map_rna'] = {base_name(args['data']): base_name(args['data']),
                                 base_name(args['data2']): base_name(args['data2'])}
    else:
        data['DATA'] = 1
        if data['ANALYZE_TYPE'] == ['TARGET_SCORE']:
            data['TARGET_SCORE_FILE'] = [base_name(args['data'])]
        elif data['ANALYZE_TYPE'] == ['RNAUP']:
            data['RNAUP_FILE'] = [base_name(args['data'])]
        elif data['ANALYZE_TYPE'] == ['PIRSCAN']:
            data['PIR_FILE'] = [base_name(args['data'])]
        elif data['ANALYZE_TYPE'] == ['G22']:
            data['G22_FILE'] = [base_name(args['data'])]
        data['title_map_rna'] = {base_name(args['data']):base_name(args['data'])}

    
    # level
    if args['level'] == 'read_count':
        data['READ_COUNT_TYPE'] = ['need_rc']
    elif args['level'] == 'site':
        data['READ_COUNT_TYPE'] = ['no_need_rc']
    

    print("\n===========  tatk-build.py  ==============")
    for arg in args:
        print("{}: {}".format(arg,args[arg]))
    print("==========================================")


    # merge configurations and input arguments
    config_data = tatk.read_config()
    config_data.update(data)
    data = config_data 
    

    # copy dataset to directory
    if data['ANALYZE_TYPE'] == ['TARGET_SCORE']:
        if data['DATA'] == 1:
            os.system("cp " + args['data'] + " bioinfo/data/CLASH/target_score_" + data['TARGET_SCORE_FILE'][0] + "_" + data['READ_COUNT_TYPE'][0] + "norm_2step.csv")
        else:
            os.system("cp " + args['data'] + " bioinfo/data/CLASH/target_score_" + data['TARGET_SCORE_FILE'][0] + "_" + data['READ_COUNT_TYPE'][0] + "norm_usetRNA_2step.csv")
            os.system("cp " + args['data2'] + " bioinfo/data/CLASH/target_score_" + data['TARGET_SCORE_FILE'][1] + "_" + data['READ_COUNT_TYPE'][0] + "norm_usetRNA_2step.csv")
        
    elif data['ANALYZE_TYPE'] == ['RNAUP']:
        os.system("cp " + args['data'] + " bioinfo/data/CLASH/RNAup_" + data['RNAUP_FILE'][0] + "_" + data['READ_COUNT_TYPE'][0] + "norm_2step.csv")
        
    elif data['ANALYZE_TYPE'] == ['PIRSCAN']:
        os.system("cp " + args['data'] + " bioinfo/data/CLASH/" + data['PIR_FILE'][0] + "_pirscan_withup.csv")
        
    if data['ANALYZE_TYPE'] == ['G22']:
        if data['DATA'] == 1:
            os.system("cp " + args['data'] + " bioinfo/data/22G/" + data['G22_FILE'][0] + "_norm_1step.csv")
        else:
            os.system("cp " + args['data'] + " bioinfo/data/22G/" + data['G22_FILE'][0] + "_EGL17M0_usemiRNA_norm_1step.csv")
            os.system("cp " + args['data2'] + " bioinfo/data/22G/" + data['G22_FILE'][1] + "_EGL17M0_usemiRNA_norm_1step.csv")


    # analyze datasets
    print("\nStart building reference...")
    outputs = tatk.analyze(data)
    print("\nExport file:")
    for output in outputs:
        print(output)
    
    
    # export configuration
    if args['output'] != None:
        with open(args['output'], 'w') as f:
            yaml.dump(data, f, default_flow_style=False)
        with open(args['output'], 'r') as f:
            text = f.read()
        with open(args['output'], 'w') as f:
            f.write("#==============  arguments  ===============\n")
            for arg in args:
                f.write("# {}: {}\n".format(arg,args[arg]))
            f.write("#============  reference data  ============\n")
            for output in outputs:
                f.write("# " + output +"\n") 
            f.write("#============  configuration  =============\n")
            f.write(text)
            f.write("#==========================================\n")
        print(os.path.abspath(args['output']))

    
    print("\nProgram end with success.\n")