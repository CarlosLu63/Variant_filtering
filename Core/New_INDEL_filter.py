import pandas as pd
from alive_progress import alive_bar
import time

class INDEL_filtering(object):
    
    def __init__(self):
        pass
    
    def filter1(self, file):
        filtered_file = file[(file.QUAL >= 200) & (file.FILTER == 'PASS')]
        return filtered_file

    def filter2(self, file):
        filtered_file = file[file.Consequence.str.contains('stop|missense|start|frameshift|inframe|splice')]
        return filtered_file

    def filter3(self, file):
        filtered_file = file[(file.EAS_AF <= 0.01) | (pd.isna(file.EAS_AF))]
        return filtered_file

    def filter4(self, file):
        filtered_file = file[(file.gnomADe_EAS_AF <= 0.01) | (pd.isna(file.gnomADe_EAS_AF))]
        return filtered_file

    def filter5(self, file):
        filtered_file = file[(file.CADD_PHRED >= 20)| (pd.isna(file.CADD_PHRED))]
        return filtered_file

    def filter6(self, file):
        filtered_file = file[(file.TWFQ <= 0.05) | (pd.isna(file.TWFQ))]
        return filtered_file

    def filter7(self, file):
        filtered_file = file[(file.CLIN_SIG_20230626.str.contains('Benign|benign')==False) | (pd.isna(file.CLIN_SIG_20230626))]
        return filtered_file

    def filter8(self, file):
        filtered_file = file[pd.notna(file.Phenotype_MIM_number)]
        return filtered_file

    def filter9(self, file):
        filtered_file = file[(file.AF_Phalanx <= 0.05) | (pd.isna(file.AF_Phalanx))]
        return filtered_file

    def filter10(self, file):
        filtered_file = file[(file.SIFT.str.contains('tolerated')==False) | (pd.isna(file.SIFT))]
        return filtered_file

    def filter11(self, file):
        filtered_file = file[(file.PolyPhen.str.contains('benign')==False) | (pd.isna(file.PolyPhen))]
        return filtered_file

    def Variants_filter(self, new_indel):
        filters = [
            self.filter1, 
            self.filter2, 
            self.filter3, 
            self.filter4, 
            self.filter5, 
            self.filter6, 
            self.filter7, 
            self.filter8, 
            self.filter9
        ]

        result = new_indel

        with alive_bar(len(filters), title = 'INDEL Filtering', bar = 'checks', spinner = 'brackets') as bar:
            for idx, filter_func in enumerate(filters, start=1):
                filtered_variants = filter_func(result)
                #print(f"Filter {idx}")
                #print(f"Remaining {len(filtered_variants)} variants\n")
                time.sleep(0.1)
                result = filtered_variants
                if len(result) == 0:
                    print(f"Stop filtering at filter {idx}")
                    status = "PASS"
                    break
                bar()
            else:
                print("Filtering Done!\n")
                status = "PASS"
        
        return result, status
