# Variant_filtering
Filter variants for WES analysis result

# Usage
python Variant_report.py 1O922112301_01A (xlsx file)

# Variant_report.py module import
## Desplay progress bar and filtering info
from New_SNV_filter import SNV_filtering
from New_INDEL_filter import INDEL_filtering
## Not to desplay progress bar and save log file 
from SNV_filter import SNV_filtering
from INDEL_filter import INDEL_filtering
sys.stdout = open(self.result_dir + case_id + '_filter_log.txt', 'w')
