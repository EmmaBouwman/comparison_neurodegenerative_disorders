from Bio.KEGG import REST
import pandas as pd
import pywikipathways as pwpw
import xml.etree.ElementTree as ET 
from unipressed import IdMappingClient
import time

tree = ET.parse("../data/spinocerebellar ataxia/WP4760.gpml")
root = tree.getroot()



protein_name_list = []
for node in root:
    if 'TextLabel' in node.keys():
        for key, value in node.attrib.items():
            if key == 'TextLabel':
                protein_name_list.append(value)



new_protein_list = [protein.replace('\n', '') for protein in protein_name_list]

request = IdMappingClient.submit(
    source="GeneCards", dest="UniProtKB", ids=new_protein_list
)
time.sleep(1.0)

df = pd.DataFrame(data=request.each_result())

df.to_excel('wikipathway_to_uniprot_SCA.xlsx')

