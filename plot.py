import argparse
import yaml
from bioinfo.piRNA_project.code import run

if __name__ == '__main__':
    
    version = "Beta 1.0"

    ### argument
    parser = argparse.ArgumentParser(
        prog="piRNA analysis",
        description="This is a tool to analyze piRNA mapping result.",
        #formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-v', '--version', action='version', version=version)
    parser.add_argument("fig", help="input directory")
    args = parser.parse_args()
    fig_name = args.fig


    ### configure
    config_path = "bioinfo/piRNA_project/code/config/"
    with open(config_path+'config.yml', 'r') as f:
        try:
            data = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(exc)

    with open(config_path+'config_gene.yml', 'r') as f:
        try:
            gene_data = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(exc)

    with open('examples/'+fig_name+'.yml', 'r') as f:
        try:
            example_data = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(exc)

    data.update(gene_data)
    data.update(example_data)
    

    ### plot fig
    run.plot(data)
