# Important:
##### Header of toxin dataset fasta file should be formated as follow:
- ids/name separated by '_'
- the last name should be the toxin code
- accepted toxin codes:
"3FTx", "5NUCL", "ACES", "BDEF", "BPP", "CNP", "CRISP", "CTL", "CVF", "CYS", "DIESTER",
                         "DIPEP", "DIS", "FA5V", "FAXV", "HYAL", "IPLA2", "KUNZ", "KUWAP", "LAO", "NGF", "OHAN",
                         "PLA2", "PLB", "SBPM", "SRTX", "SVLIPA", "SVMI", "SVMMP", "SVMP", "SVSP", "VEGF", "WAP"

##### * The final result is totally dependent of the toxin dataset provided as input.

##### Example:
```
python3 toxin_hunter.py /path/assembly/file.fasta /path/toxin/dataset/file.fasta 90 Sample_name
```
