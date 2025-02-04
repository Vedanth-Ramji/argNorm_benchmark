import os
import pandas as pd
from argnorm.drug_categorization import confers_resistance_to, drugs_to_drug_classes
from argnorm.lib import get_aro_ontology

ARO = get_aro_ontology()

def generate_hits_tsv():
    mappings = sorted(os.listdir('./rgi_mapping/'))
    output = []
    for i in mappings:
        df = pd.read_csv('./rgi_mapping/' + i, sep='\t')
        db = i.split('_')[0]
        df['Database'] = db
        output.append(df[['ORF_ID', 'ARO', 'Database', 'Cut_Off']])
    output = pd.concat(output).reset_index(drop=True)

    drugs = ('ARO:' + output['ARO'].map(str)).map(confers_resistance_to)
    drug_classes = drugs.map(drugs_to_drug_classes)

    output['Drugs'] = drugs.map(lambda x: list(map(lambda y: ARO[y].name, x)))
    output['Drug Classes'] = drug_classes.map(lambda x: list(map(lambda y: ARO[y].name, x)))

    original_drug_classes = []

    argannot = output.query('Database == "argannot"')
    original_drug_classes.append(argannot['ORF_ID'].str.split(')').str[0].str[1:])

    ncbi = output.query('Database == "ncbi"')
    ncbi_drug_classes = ncbi['ORF_ID'].map(lambda x: x.split('|')[-2] if '|' in x else x)
    ncbi_drug_classes.where(ncbi_drug_classes != '', ncbi['ORF_ID'], inplace=True)
    original_drug_classes.append(ncbi_drug_classes)

    deeparg = output.query('Database == "deeparg"')
    original_drug_classes.append(deeparg['ORF_ID'].str.split('|').str[-2])

    resfinder = output.query('Database == "resfinder"')
    resfinder_antibiotic_classes = pd.read_csv('./data/resfinder_antibiotic_classes.tsv', sep='\t')
    resfinder_antibiotic_classes.set_index('Gene_accession no.', inplace=True)
    resfinder_antibiotic_classes = resfinder_antibiotic_classes['Class'].to_dict()
    resfinder_drug_classes = resfinder['ORF_ID'].map(resfinder_antibiotic_classes)
    original_drug_classes.append(resfinder_drug_classes)

    resfinderfg = output.query('Database == "resfinderfg"')
    original_drug_classes.append(resfinderfg['ORF_ID'].str.split('|').str[0])

    megares = output.query('Database == "megares"')
    original_drug_classes.append(megares['ORF_ID'].str.split('|').str[2])

    sarg = output.query('Database == "sarg"')
    gene_name = sarg["ORF_ID"].str.split(' ').str[0]
    sarg_antibiotic_classes = pd.read_csv('./data/SARG_structure.tsv', sep='\t')
    sarg_antibiotic_classes.set_index('SARG.Seq.ID', inplace=True)
    sarg_antibiotic_classes = sarg_antibiotic_classes['Type'].to_dict()
    sarg_drug_classes = gene_name.map(sarg_antibiotic_classes)
    original_drug_classes.append(sarg_drug_classes)

    original_drug_classes = pd.concat(original_drug_classes)

    output['Original Drug Classes'] = original_drug_classes.reindex(output.index)
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
        "nucleoside antibiotic": ['puromycin', 'Nucleosides', 'tunicamycin'],
        "rifamycin antibiotic": ['rifampin'],
        "lincosamide antibiotic": ['lincosamide', 'MLS'],
        "streptogramin antibiotic": ['streptogramin', 'streptogramin A', 'streptogramin B', 'MLS'],
        'nitroimidazole antibiotic': ['Metronidazole', 'Ntmdz'],
        'oxazolidinone antibiotic': ['Oxzln'],
        'fusidane antibiotic': ['fusidic_acid', 'fusaric-acid', 'FUSIDIC_ACID', 'fusidic-acid', 'Fcd'],
        'fluoroquinolone antibiotic': ['Flq', 'Fluoroquinolones', 'PHENICOL/QUINOLONE'],
        'pleuromutilin antibiotic': ['pleuromutilin_tiamulin', 'LINCOSAMIDE/PLEUROMUTILIN']
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