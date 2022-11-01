import argparse
import os
import sys
from bioinfo.piRNA_project.code import main

__version__ = "version 1.0"


def base_name(file_path):
    # exmaple: /path/to/some/file.txt -> file
    return os.path.splitext(os.path.basename(file_path))[0]


def collect_files(path):
    file_list = []
    for dirs, subdirs, files in os.walk(path):
        for f in files:
            file_list.append(os.path.join(dirs, f))
    return file_list


if __name__ == '__main__':

    # arguments
    parser = argparse.ArgumentParser(
        description="This program is to plot the result for RNA-seq Read Distribution Toolkit (RDT).",
    )
    parser.add_argument("-v", "--version",
                        action="version",
                        version="%(prog)s " + __version__)
    parser.add_argument("-i", "--input",
                        help="read metadata in YAML format")
    parser.add_argument("-o", "--output",
                        help="output figure in png format")
    parser.add_argument("-L",
                        action="store_true",
                        help="compare two transcripts")
    parser.add_argument("--list", 
                        help="file of transcript list")
    parser.add_argument("--list2",
                        help="file of another transcript list")
    parser.add_argument("-T", "--title",
                        help="set the title of image")
    #parser.add_argument("-s", "--svg"
    #                    action="store_true",
    #                    help="save image into Scalable Vector Graphics (svg) ")
    args = vars(parser.parse_args())

    data = {}

    # check if some necessary arguments are exist
    nec_args = args.copy()
    for arg in ['L', 'list', 'list2', 'title']:
        nec_args.pop(arg)
    program_exit = False
    for arg in nec_args:    
        if nec_args[arg] == None:
            print('[Error] Argument "--' + arg + '" should not left blank.')
            program_exit = True
    if program_exit:
        sys.exit()


    # list
    data['title_map_gene'] = { 1: 'all mRNAs' }
    if args['L'] == True:
        data['LIST'] = 2
        if args['list'] == None:
            data['GENE_LIST1'] = [1]
        else:
            data['GENE_LIST1'] = [2]
            data['title_map_gene'][2] = base_name(args['list'])
        if args['list2'] == None:
            data['GENE_LIST2'] = [1]
        else:
            data['GENE_LIST2'] = [3]
            data['title_map_gene'][3] = base_name(args['list2'])
    else:
        data['LIST'] = 1
        if args['list'] == None:
            data['ANALYZE_GROUP'] = [1]
        else:
            data['ANALYZE_GROUP'] = [2]
            data['title_map_gene'][2] = base_name(args['list'])


    # title
    if args['title'] == None:
        data['img_title'] = base_name(args['output'])
    else:
        data['img_title'] = args['title']


    # svg
    #if args['svg']:
    #    data['img_format'] = "svg"
    #else:
    #    data['img_format'] = "png"


    print("\n===========  rdt-plot.py  ==============")
    for arg in args:
        print("{}: {}".format(arg,args[arg]))
    print("==========================================")


    # merge configurations and input arguments
    config_data = main.read_config(args['input'])
    config_data.update(data)
    data = config_data 


    # copy list to directory
    data['gene_path'] = "../../data/gene_list"
    list_files = []
    for trans in [ args['list'], args['list2'] ]:
        if trans != None:
            os.system("cp " + trans + " bioinfo/data/gene_list/" + base_name(trans))
            list_files.append(os.path.abspath("bioinfo/data/gene_list/" + base_name(trans)))


    # make directory to store image
    if os.path.isdir("static"):
        os.system("rm -r static")
    os.mkdir("static")
    os.mkdir("static/paper")
    for dirs in ['density_fig','meta_fig','codon_fig','fold_change_fig']:
        os.mkdir(os.path.join("static/paper/",dirs))
        for subdirs in ['single','compare','compare_list']:
            os.mkdir(os.path.join("static/paper/",dirs,subdirs))


    # plot figure
    print("\nStart ploting figure...")
    main.plot(data)
    img_lst = collect_files("static")
    for img in img_lst:
        os.rename(img, args['output'])
    os.system("rm -r static")

    print("\nExport file:")
    for list_file in list_files:
        print(list_file)
    print(os.path.abspath(args['output']))

    print("\nProgram end with success.\n")