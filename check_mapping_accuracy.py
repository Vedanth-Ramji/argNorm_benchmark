import os
import pandas as pd
from argnorm.drug_categorization import *
import pronto

ARO = pronto.Ontology('./data/aro.obo')

def generate_hits_tsv():
    mappings = os.listdir('./rgi_mapping/')
    outputs = []

    for i in mappings:
        df = pd.read_csv('./rgi_mapping/' + i, sep='\t')
        db = i.split('_')[0]
        df['Database'] = db
        outputs.append(df[['ORF_ID', 'ARO', 'Database', 'Cut_Off']])
    
    output = pd.concat(outputs)
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
    
    resfinder_antibiotic_classes = pd.read_csv('./data/resfinder_antibiotic_classes.tsv', sep='\t')
    sarg_antibiotic_classes = pd.read_csv('./data/SARG_structure.tsv', sep='\t')

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
            drug_class = str(resfinder_antibiotic_classes[resfinder_antibiotic_classes['Gene_accession no.'] == gene_name]['Class'].values).replace("['", '').replace("']", '')
        if output.iloc[i]['Database'] == 'resfinderfg':
            drug_class = output.iloc[i]['ORF_ID'].split('|')[0]
        if output.iloc[i]['Database'] == 'megares':
            drug_class = output.iloc[i]['ORF_ID'].split('|')[2]
        if output.iloc[i]['Database'] == 'sarg':
            gene_name = output.iloc[i]['ORF_ID'].split(' ')[0]
            drug_class = str(sarg_antibiotic_classes[sarg_antibiotic_classes['SARG.Seq.ID'] == gene_name]['Type'].values).replace("['", '').replace("']", '')
        
        drug_classes.append(drug_class)

    output['Original Drug Classes'] = drug_classes
    output.to_csv('hits.tsv', sep='\t', index=False)

def analyze_hits_tsv():
    alternative_ids = {
        "beta-lactam antibiotic": ['Methicillin resistance mecR1 protein', 'Bla', 'beta_lactam', 'beta-lactam', 'Beta-lactamase', 'betalactams', 'beta-lactamase', 'Metallo-beta-lactamase', 'putative peptidoglycan D%2CD-transpeptidase PenA'],
        "aminoglycoside antibiotic": ['aminoglycoside', 'AGly', 'Aminoglycosides', "Streptomycin 3''-adenylyltransferase", 'Gentamicin 3-N-acetyltransferase', 'Bifunctional AAC/APH'],
        "macrolide antibiotic": ['MLS', 'MACROLIDE', 'Macrolide'],
        "tetracycline antibiotic": ['tetracycline', 'Tet', 'TETRACYCLINE', 'Tetracyclines', 'Tetracycline', 'Tetracycline resistance protein', 'Tetracycline repressor protein'],
        "phenicol antibiotic": ['Phe', 'chloramphenicol', 'amphenicol', 'Chloramphenicol acetyltransferase', 'Chloramphenicol acetyltransferase 2', 'florfenicol'],
        "phosphonic acid antibiotic": ['fosfomycin', 'Fosfomycin', 'Fcyn'],
        "sulfonamide antibiotic": ['Sulfonamides', 'Dihydropteroate synthase', 'Folate pathway antagonist'],
        "glycopeptide antibiotic": ['Glycopeptides', 'vancomycin', 'bleomycin', 'D-alanine--D-alanine ligase', 'D-alanine--D-alanine ligase B', 'D-alanine--D-alanine ligase A'],
        "diaminopyrimidine antibiotic": ['Dihydrofolate reductase', 'Trimethoprim', 'Folate pathway antagonist', 'Tmt'],
        "peptide antibiotic": ['bacitracin', 'polymyxin', 'other_peptide_antibiotics', 'COLISTIN', 'lipopeptides', 'COL', 'TUBERACTINOMYCIN', 'edeine', 'defensin', 'Cationic_antimicrobial_peptides'],
        "aminocoumarin antibiotic": ['novobiocin'],
        "nucleoside antibiotic": ['puromycin', 'Nucleosides', 'tunicamycin', 'streptothricin', 'aminoglycoside', 'AGly'],
        "rifamycin antibiotic": ['rifampin'],
        "lincosamide antibiotic": ['lincosamide', 'MLS'],
        "streptogramin antibiotic": ['streptogramin', 'streptogramin A', 'streptogramin B', 'MLS'],
        'nitroimidazole antibiotic': ['Metronidazole', 'Ntmdz'],
        'oxazolidinone antibiotic': ['Oxzln'],
        'fusidane antibiotic': ['fusidic_acid', 'fusaric-acid', 'FUSIDIC_ACID', 'fusidic-acid', 'Fcd'],
        'fluoroquinolone antibiotic': ['Flq', 'Fluoroquinolones', 'PHENICOL/QUINOLONE'],
        'pleuromutilin antibiotic': ['pleuromutilin_tiamulin', 'LINCOSAMIDE/PLEUROMUTILIN'],
    }
    
    metals = ['mercury_resistance', 'multi-metal_resistance', 'tellurium_resistance', 'tellurium', 'arsenic', 'cadmium', 'copper', 'mercury', 'nickel', 'copper/silver', 'silver', 'cadmium/cobalt/nickel', 'chromate', 'COPPER/GOLD', 'GOLD']
    virulence_genes_or_toxins = ['stx2', 'intimin', 'stx1']

    df = pd.read_csv('hits.tsv', sep='\t')

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
            elif 'multidrug' in str(df.iloc[i]['Original Drug Classes']).lower() or 'multi-drug_resistance' in str(df.iloc[i]['Original Drug Classes']).lower():
                correct.append(df.iloc[i])
            elif 'macrolide-lincosamide-streptogramin' in str(df.iloc[i]['Original Drug Classes']).lower():
                matched = False
                for ii in 'macrolide-lincosamide-streptogramin'.split('-'):
                    if ii in df.iloc[i]['Drug Classes'].lower():
                        correct.append(df.iloc[i])
                        matched = True
                        break
                
                if not matched:
                    incorrect.append(df.iloc[i])
            elif 'sdia' in str(df.iloc[i]['ORF_ID']).lower() or 'cpxr' in str(df.iloc[i]['ORF_ID']).lower() or 'rosA' in str(df.iloc[i]['ORF_ID']) or 'rosB' in str(df.iloc[i]['ORF_ID']):
                correct.append(df.iloc[i])
            elif 'penicillin-binding_protein_' in str(df.iloc[i]['ORF_ID']).lower() and 'beta-lactam antibiotic' in str(df.iloc[i]['Drug Classes']):
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

    with open('loose_hit_outputs/flagged_loose_hits.txt', 'w') as ofile:
        for i in incorrect:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\t{i['Cut_Off']}\n")
            
    with open('loose_hit_outputs/matched_loose_hits.txt', 'w') as ofile:
        for i in correct:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\t{i['Cut_Off']}\n")
            
    with open('loose_hit_outputs/loose_hit_metal_genes.txt', 'w') as ofile:
        for i in metal_resistance_genes:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\t{i['Cut_Off']}\n")
            
    with open('loose_hit_outputs/loose_hit_virulence_genes.txt', 'w') as ofile:
        for i in virulence_genes:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\t{i['Cut_Off']}\n")
            
    with open('loose_hit_outputs/loose_hit_drug_and_biocide_genes.txt', 'w') as ofile:
        for i in drug_and_biocide_genes:
            ofile.write(f"{i['ORF_ID']}\t{i['Drug Classes']}\t{i['Original Drug Classes']}\t{i['Cut_Off']}\n")

generate_hits_tsv()
analyze_hits_tsv()