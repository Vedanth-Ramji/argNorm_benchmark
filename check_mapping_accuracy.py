import os
import pandas as pd
from argnorm.drug_categorization import *
import pronto

ARO = pronto.Ontology('./data/aro.obo')

def get_loose_hits(df, db):
    loose = df[df['Cut_Off'].isin(['Loose'])]
    loose['Database'] = db
    return loose[['ORF_ID', 'ARO', 'Database']]

def generate_loose_hits_tsv():
    mappings = os.listdir('./rgi_mapping/')
    output = pd.DataFrame()

    for i in mappings:
        df = pd.read_csv('./rgi_mapping/' + i, sep='\t')
        db = i.split('_')[0]
        loose = get_loose_hits(df, db)
        output = pd.concat([output, loose])

    drugs_list = []
    drug_classes_list = [] 
    for i in range(output.shape[0]):
        aro = "ARO:" + str(output.iloc[i]['ARO'])
        drugs = confers_resistance_to(aro)
        drug_classes = drugs_to_drug_classes(drugs)
        
        drugs_list.append(list(map(lambda x: ARO[x].name, drugs)))
        drug_classes_list.append(list(map(lambda x: ARO[x].name, drug_classes)))
        
    output['Drugs'] = drugs_list
    output['Drug Classes'] = drug_classes_list

    drug_classes = []
    for i in range(output.shape[0]):
        if output.iloc[i]['Database'] == 'argannot':
            drug_class = output.iloc[i]['ORF_ID'].split(')')[0][1:]
        if output.iloc[i]['Database'] == 'ncbi':
            if output.iloc[i]['ORF_ID'].split('|')[-2]:
                drug_class = output.iloc[i]['ORF_ID'].split('|')[-2]
            else:
                drug_class = output.iloc[i]['ORF_ID']
        if output.iloc[i]['Database'] == 'deeparg':
            drug_class = output.iloc[i]['ORF_ID'].split('|')[-2]
        if output.iloc[i]['Database'] == 'resfinder':
            gene_name = output.iloc[i]['ORF_ID']
            resfinder_antibiotic_classes = pd.read_csv('./data/resfinder_antibiotic_classes.tsv', sep='\t')
            drug_class = resfinder_antibiotic_classes[resfinder_antibiotic_classes['Gene_accession no.'] == gene_name]['Class'].to_string().split(' ')[-1]
        if output.iloc[i]['Database'] == 'resfinderfg':
            drug_class = output.iloc[i]['ORF_ID'].split('|')[0]
        if output.iloc[i]['Database'] == 'megares':
            drug_class = output.iloc[i]['ORF_ID'].split('|')[2]
        if output.iloc[i]['Database'] == 'sarg':
            gene_name = output.iloc[i]['ORF_ID'].split(' ')[0]
            sarg_antibiotic_classes = pd.read_csv('./data/SARG_structure.tsv', sep='\t')
            drug_class = sarg_antibiotic_classes[sarg_antibiotic_classes['SARG.Seq.ID'] == gene_name]['Type'].to_string().split(' ')[-1]
        
        drug_classes.append(drug_class)

    output['Original Drug Classes'] = drug_classes
    output.to_csv('loose_hits.tsv', sep='\t', index=False)

def analyze_loose_hits_tsv():
    alternative_ids = {
        "beta-lactam antibiotic": ['Bla', 'beta_lactam', 'beta-lactam', 'Beta-lactamase', 'betalactams', 'beta-lactamase', 'Metallo-beta-lactamase'],
        "aminoglycoside antibiotic": ['Aminoglycoside', 'AMINOGLYCOSIDE', 'aminoglycoside', 'AGly', 'Aminoglycosides'],
        "macrolide antibiotic": ['MLS', 'MACROLIDE', 'Macrolide'],
        "tetracycline antibiotic": ['tetracycline', 'Tet', 'TETRACYCLINE', 'Tetracyclines', 'Tetracycline', 'Tetracycline resistance protein', 'Tetracycline repressor protein'],
        "phenicol antibiotic": ['Phe', 'chloramphenicol', 'amphenicol', 'Chloramphenicol acetyltransferase', 'Chloramphenicol acetyltransferase 2', 'florfenicol']
    }
    
    metals = ['mercury_resistance', 'multi-metal_resistance', 'tellurium_resistance', 'tellurium', 'arsenic', 'cadmium', 'copper', 'mercury', 'nickel', 'copper/silver', 'silver', 'cadmium/cobalt/nickel', 'chromate']
    virulence_genes_or_toxins = ['stx2', 'intimin', 'stx1']

    df = pd.read_csv('loose_hits.tsv', sep='\t')

    total = df.shape[0]
    correct = []
    incorrect = []
    metal_resistance_genes = []
    virulence_genes = []
    drug_and_biocide_genes = [] 
    
    for i in range(df.shape[0]):
        if str(df.iloc[i]['Original Drug Classes']).lower() == df.iloc[i]['Drug Classes'].lower():
            correct.append(df.iloc[i])
        else:
            if str(df.iloc[i]['Original Drug Classes']).lower() in df.iloc[i]['Drug Classes'].lower():
                correct.append(df.iloc[i])
            elif 'multidrug' in str(df.iloc[i]['Original Drug Classes']).lower():
                correct.append(df.iloc[i])
            elif 'macrolide-lincosamide-streptogramin' == str(df.iloc[i]['Original Drug Classes']).lower():
                matched = False
                for ii in 'macrolide-lincosamide-streptogramin'.split('-'):
                    if ii in df.iloc[i]['Drug Classes'].lower():
                        correct.append(df.iloc[i])
                        matched = True
                        break
                
                if not matched:
                    incorrect.append(df.iloc[i])
            elif 'sdia' in str(df.iloc[i]['ORF_ID']).lower() or 'cpxr' in str(df.iloc[i]['ORF_ID']).lower():
                correct.append(df.iloc[i])
            elif str(df.iloc[i]['Original Drug Classes']).lower() in metals:
                metal_resistance_genes.append(df.iloc[i])
            elif str(df.iloc[i]['Original Drug Classes']).lower() in virulence_genes_or_toxins:
                virulence_genes.append(df.iloc[i])
            elif str(df.iloc[i]['Original Drug Classes']) == "Drug_and_biocide_resistance":
                drug_and_biocide_genes.append(df.iloc[i])
            else:
                matched = False
                for id in alternative_ids:
                    if id in df.iloc[i]['Drug Classes']:
                        for alternative in alternative_ids[id]:
                            if str(alternative).lower() in str(df.iloc[i]['Original Drug Classes']).lower():
                                correct.append(df.iloc[i])
                                matched = True
                                break
                        
                if not matched:
                        incorrect.append(df.iloc[i])

    with open('loose_hit_outputs/incorrect_loose_hits.txt', 'w') as ofile:
        for i in incorrect:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\n")
            
    with open('loose_hit_outputs/correct_loose_hits.txt', 'w') as ofile:
        for i in correct:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\n")
            
    with open('loose_hit_outputs/loose_hit_metal_genes.txt', 'w') as ofile:
        for i in metal_resistance_genes:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\n")
            
    with open('loose_hit_outputs/loose_hit_virulence_genes.txt', 'w') as ofile:
        for i in virulence_genes:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\n")
            
    with open('loose_hit_outputs/loose_hit_drug_and_biocide_genes.txt', 'w') as ofile:
        for i in drug_and_biocide_genes:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\n")
    
analyze_loose_hits_tsv()