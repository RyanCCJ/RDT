import yaml

# density
from bioinfo.piRNA_project.code.pkg.tool1_1 import tool1_1
from bioinfo.piRNA_project.code.pkg.tool1_2 import tool1_2
from bioinfo.piRNA_project.code.pkg.tool1_3 import tool1_3
from bioinfo.piRNA_project.code.pkg.tool1_4 import tool1_4
# metagene
from bioinfo.piRNA_project.code.pkg.tool2_1 import tool2_1
from bioinfo.piRNA_project.code.pkg.tool2_2 import tool2_2
from bioinfo.piRNA_project.code.pkg.tool2_3 import tool2_3
from bioinfo.piRNA_project.code.pkg.tool2_4 import tool2_4
# codon
from bioinfo.piRNA_project.code.pkg.tool3_1 import tool3_1
from bioinfo.piRNA_project.code.pkg.tool3_2 import tool3_2
from bioinfo.piRNA_project.code.pkg.tool3_3 import tool3_3
from bioinfo.piRNA_project.code.pkg.tool3_4 import tool3_4
# fold change
from bioinfo.piRNA_project.code.pkg.tool4_4 import tool4_4


def analyze(data):
    for tool in data['TOOL']:
        if tool == 1:
            print("==== Density Analysis ====")
            output = tool1_1(data)
        elif tool == 2:
            print("==== Metagene Analysis ====")
            output = tool2_1(data)
        elif tool == 3:
            print("==== Codon Analysis ====")
            output = tool3_1(data)
        elif tool == 4:
            print("==== Fold Change Analysis ====")
            output = tool1_1(data)
    print("=== analyze done ===")
    return output


def plot(data):
    for tool in data['TOOL']:

        if tool == 1:
            print("==== Density Plot ====")
            if data['DATA'] == 2:
                tool1_4(data)
            elif data['LIST'] == 1:
                tool1_2(data)
            elif data['LIST'] == 2:
                tool1_3(data)

        elif tool == 2:
            print("==== Metagene Plot ====")
            if data['DATA'] == 2:
                tool2_4(data)
            elif data['LIST'] == 1:
                tool2_2(data)
            elif data['LIST'] == 2:
                tool2_3(data)

        elif tool == 3:
            print("==== Codon Plot ====")
            if data['DATA'] == 2:
                tool3_4(data)
            elif data['LIST'] == 1:
                tool3_2(data)
            elif data['LIST'] == 2:
                tool3_3(data)
        
        elif tool == 4:
            print("==== Fold Change Plot ====")
            if data['DATA'] == 2:
                tool4_4(data)
    
    print("=== plot done ===")


def read_config(new_config=None):
    
    with open('bioinfo/piRNA_project/code/config/config.yml', 'r') as f:
        try:
            data = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(exc)
    
    if new_config != None:
        with open(new_config, 'r') as f:
            try:
                new_data = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                print(exc)
        data.update(new_data)

    return data
