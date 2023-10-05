import pandas as pd

class INDEL_filtering(object):
    
    def __init__(self):
        pass
    
    def filter_1(self, file):
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

    def Variants_filter(self, new_snv):
        if len(self.filter_1(new_snv)) > 0:
            result = self.filter_1(new_snv)
            print(f"Filter 1")
            print(f"Remaining {len(result)} variants\nGoing to filter 2\n")    
            if len(self.filter2(result)) > 0:
                result = self.filter2(result)
                print(f"Filter 2")
                print(f"Remaining {len(result)} variants\nGoing to filter 3\n")      
                if len(self.filter3(result)) > 0:
                    result = self.filter3(result)
                    print(f"Filter 3")
                    print(f"Remaining {len(result)} variants\nGoing to filter 4\n")
                    if len(self.filter4(result)) > 0:
                        result = self.filter4(result)
                        print(f"Filter 4")
                        print(f"Remaining {len(result)} variants\nGoing to filter 5\n")
                        if len(self.filter5(result)) > 0:
                            result = self.filter5(result)
                            print(f"Filter 5")
                            print(f"Remaining {len(result)} variants\nGoing to filter 6\n")
                            if len(self.filter6(result)) > 0:
                                result = self.filter6(result)
                                print(f"Filter 6")
                                print(f"Remaining {len(result)} variants\nGoing to filter 7\n")
                                if len(self.filter7(result)) > 0:
                                    result = self.filter7(result)
                                    print(f"Filter 7")
                                    print(f"Remaining {len(result)} variants\nGoing to filter 8\n")
                                    if len(self.filter8(result)) > 0:
                                        result = self.filter8(result)
                                        print(f"Filter 8")
                                        print(f"Remaining {len(result)} variants\nGoing to filter 9\n")
                                        if len(self.filter9(result)) > 0:
                                            result = self.filter9(result)
                                            status = "PASS"
                                            print(f"Filter 9")
                                            print(f"Remaining {len(result)} variants\nFiltering Done!\n")
                                        else:
                                            result = self.filter9(result)
                                            status = "PASS"
                                            print(f"Filter 9")
                                            print(f"Remaining {len(result)} variants\nStop filtering at filter 9")
                                    else:
                                        result = self.filter8(result)
                                        status = "PASS"
                                        print(f"Filter 8")
                                        print(f"Remaining {len(result)} variants\nStop filtering at filter 8")
                                else:
                                    result = self.filter7(result)
                                    status = "PASS"
                                    print(f"Filter 7")
                                    print(f"Remaining {len(result)} variants\nStop filtering at filter 7")
                            else:
                                result = self.filter6(result)
                                status = "PASS"
                                print(f"Filter 6")
                                print(f"Remaining {len(result)} variants\nStop filtering at filter 6")
                        else:
                            result = self.filter5(result)
                            status = "PASS"
                            print(f"Filter 5")
                            print(f"Remaining {len(result)} variants\nStop filtering at filter 5")
                    else:
                        result = self.filter4(result)
                        status = "PASS"
                        print(f"Filter 4")
                        print(f"Remaining {len(result)} variants\nStop filtering at filter 4")
                else:
                    result = self.filter3(result)
                    status = "PASS"
                    print(f"Filter 3")
                    print(f"Remaining {len(result)} variants\nStop filtering at filter 3")
            else:
                result = self.filter2(result)
                status = "PASS"
                print(f"Filter 2")
                print(f"Remaining {len(result)} variants\nStop filtering at filter 2")
        else:
            result = self.filter_1(new_snv)
            status = "PASS"
            print(f"Filter 1")
            print(f"Remaining {len(result)} variants\nStop filtering at filter 1")
        
        return result, status
