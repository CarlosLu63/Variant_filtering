import sys
import os
import pandas as pd
import numpy as np
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from SNV_filter import SNV_filtering
from INDEL_filter import INDEL_filtering

#wdir = '/share/NGS_Project/1O923090801_woBED/'

class Variant_filtering(object):
    
    def main_flow(self, case_id):

        #sys.stdout = open(self.result_dir + case_id + '_filter_log.txt', 'w')
        
        os.chdir(self.wdir)
        print(f"START SNV Filtering...\n")
        print(f"Reading {case_id} raw excel file...\n")
        raw_snv = pd.read_excel(case_id + '.xlsx', sheet_name = 'SNV', skiprows = 1)
        new_snv = raw_snv.replace('.', np.nan)
        try:  
            SNV_result, status = SNV_filtering().Variants_filter(new_snv)
        except Exception as e:
             print(e)
        #final_data.to_excel(case_id + "_filtered.xlsx", index = False)

        print(f"START INDEL Filtering...\n")
        print(f"Reading {case_id} raw excel file...\n")
        raw_indel = pd.read_excel(case_id + '.xlsx', sheet_name = 'INDEL', skiprows = 1)
        new_indel = raw_indel.replace('.', np.nan)
        try:
            INDEL_result, status = INDEL_filtering().Variants_filter(new_indel)
        except Exception as e:
             print(e)
        #final_data.to_excel(case_id + "_filtered.xlsx", index = False)

        #print(f"Saving results to excel file...\n")
        wb = Workbook()
        wb.remove(wb['Sheet'])
        for item, df in zip(("SNV", "INDEL"), (SNV_result, INDEL_result)):
            print(f"Saving {item} results to excel file...\n")
            wb.create_sheet(item)
            for row in dataframe_to_rows(df, index = False):
                wb[item].append(row)
            wb.save(self.result_dir + case_id + '_filtered.xlsx')

        #sys.stdout.close()

    def __init__(self):
        self.wdir = 'D:/WES_LDTS/報告解讀/test_data/'
        self.result_dir = 'D:/WES_LDTS/報告解讀/Filter_results/'

def main():	 
	case_id = sys.argv[1]  
	main_obj = Variant_filtering()
	main_obj.main_flow(case_id)

if __name__ == "__main__":
### simple main
	main()

#################################################################################
#################################################################################
#################################################################################
"""
fil_1 = new_snv[(new_snv.QUAL >= 200) & (new_snv.FILTER == 'PASS')]
print(f"Remaining {len(fil_1)} variants...")

fil_2 = fil_1[fil_1.Consequence.str.contains('stop|missense|start|frameshift|inframe|splice')]
print(f"Remaining {len(fil_2)} variants...")

fil_3 = fil_2[(fil_2.EAS_AF <= 0.01) | (pd.isna(fil_2.EAS_AF))]
print(f"Remaining {len(fil_3)} variants...")

fil_4 = fil_3[(fil_3.gnomADe_EAS_AF <= 0.01) | (pd.isna(fil_3.gnomADe_EAS_AF))]
print(f"Remaining {len(fil_4)} variants...")

fil_5 = fil_4[(fil_4.CADD_PHRED >= 20)| (pd.isna(fil_4.CADD_PHRED))]
print(f"Remaining {len(fil_5)} variants...")

fil_6 = fil_5[(fil_5.TWFQ <= 0.05) | (pd.isna(fil_5.TWFQ))]
print(f"Remaining {len(fil_6)} variants...")

fil_7 = fil_6[(fil_6.CLIN_SIG_20230626.str.contains('Benign|benign')==False) | (pd.isna(fil_6.CLIN_SIG_20230626))]
print(f"Remaining {len(fil_7)} variants...")

fil_8 = fil_7[pd.notna(fil_7.Phenotype_MIM_number)]
print(f"Remaining {len(fil_8)} variants...")

fil_9 = fil_8[(fil_8.AF_Phalanx <= 0.05) | (pd.isna(fil_8.AF_Phalanx))]
print(f"Remaining {len(fil_9)} variants...")

fil_10 = fil_9[(fil_9.SIFT.str.contains('tolerated')==False) | (pd.isna(fil_9.SIFT))]
print(f"Remaining {len(fil_10)} variants...")

fil_11 = fil_10[(fil_10.PolyPhen.str.contains('benign')==False) | (pd.isna(fil_10.PolyPhen))]
print(f"Remaining {len(fil_11)} variants...")
"""
