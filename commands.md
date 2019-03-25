
# Set up new test directory with raw data files.
```shell
$ mkdir test
$ cp biomart_data/mart_export.txt targetscan_60_data/miR_Family_Info.txt test
$ cd test
```

# Format miRNA list for TargetScan
```shell
$ awk '{ print $4,$2,$3 }' OFS='\t' miR_Family_Info.txt | sort | tail -n +2 > miR_input.txt
```
## remove duplicates
```shell
$ sort -u -k 2 miR_input.txt | sort -k 1 > miR_input2.txt
```

# Clean BioMart data
```shell
# One line per gene 
$ awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' mart_export.txt | \
# convert to fasta
bioawk -c fastx '{print ">"$name,$seq}' OFS='\t' | \
# Replace pipe between fields with tab 
gsed 's/|/\t/g' | \
# Remove missing sequences 
grep -v "Sequence unavailable" | \
# Remove extrachomosomal scaffolds & write to file
grep -v "K[N,Z]" > biomart_UTRs.txt
```
## in one line:
```shell
$ awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' mart_export.txt | bioawk -c fastx '{print ">"$name,$seq}' OFS='\t' | gsed 's/|/\t/g' | grep -v "Sequence unavailable" | grep -v "K[N,Z]" | awk 'length($6) > 6 {print}' > biomart_UTRs.txt
```

# keep-longest.py
Removes duplicate UTRs and only keeps the longest one (https://gist.github.com/mkweskin/8869358).
```shell
# Generate fasta file for python script 
$ cut -f1,6 biomart_UTRs.txt | awk '{print $1, "\n", $2}' OFS='' > UTRs.fasta

$ python ../keep-longest.py UTRs.fasta > kl-output.fasta
```

# Format for fix-positions.R
```shell
# Convert fasta file of unique entries back to tsv (1 line per UTR)
$ bioawk -c fastx '{print ">"$name,$seq}' OFS='\t' kl-output.fasta > unique_UTRs.txt

# Join list of unique sequences to genomic position information
$ join -1 2 -2 6 -o1.1,2.1,2.2,2.3,2.4,2.5,2.6 -t $'\t' <(sort -k2 unique_UTRs.txt) <(sort -k6 biomart_UTRs.txt) | awk '$1 == $2 {print}' > out.txt

## Deal with duplicates
# Removed the 29 duplicate lines in out.txt manually. These pairs of lines were identical except that the genomic positions were reversed. To do this programmatically: if strand == 1, keep line where field 5 (before the ';') is lesser. If strand == 2, keep line where field 5 (before the ';') is greater. File containing manually selected duplicate UTRs is mar16/duplicates2.txt.
# Identify lines in out.txt but not in unique_UTRs
$ comm -2 -3 <(cut -f1 out.txt | sort) <(cut -f1 unique_UTRs.txt | sort) > onlyInOut.txt
# File of 29 duplicate pairs (58 lines):
$ grep -f onlyInOut.txt out.txt > duplicates.txt
# Output with the 58 duplicate lines removed:
$ grep -v -f onlyInOut.txt out.txt > out-dups-removed.txt
# Add the correct 29 lines back to the output file
$ cat duplicates2.txt >> out-dups-removed.txt

# Clean file for R script
$ cut -f2-7 out-dups-removed.txt | sed 's/>//g' > fix_pos-input.txt
```

## Expected line counts at this point:
```shell
$ wc -l *
   56912 UTRs.fasta
   28456 biomart_UTRs.txt
      58 duplicates.txt
      29 duplicates2.txt
   41102 kl-output.fasta
  538996 mart_export.txt
     215 miR_input.txt
     126 miR_input2.txt
      29 onlyInOut.txt
   20551 out-dups-removed.txt
   20580 out.txt
   20551 unique_UTRs.txt
  727235 total
```

## fix-positions.R

```R
R ../R/fix-positions.R fix_pos-input.txt
```

## Validate fix-positions.R output
```shell
# Check how many new lines were added:
$ wc -l fix_pos-*
   20551 fix_pos-input.txt
   22659 fix_pos-output.txt
# 2108 new lines were added
# Get counts for each number of semicolons in a single line:
$ cut -f4 fix_pos-input.txt | grep ";" | sed -E 's/[0-9]+//g' | sort | uniq -c
 934 ;
 104 ;;
  59 ;;;
  43 ;;;;
  29 ;;;;;
  13 ;;;;;;
   8 ;;;;;;;
  12 ;;;;;;;;
   3 ;;;;;;;;;
   1 ;;;;;;;;;;
   3 ;;;;;;;;;;;
   2 ;;;;;;;;;;;;
   1 ;;;;;;;;;;;;;;
   2 ;;;;;;;;;;;;;;;;
   1 ;;;;;;;;;;;;;;;;;;
   1 ;;;;;;;;;;;;;;;;;;;;;;;;
   1 ;;;;;;;;;;;;;;;;;;;;;;;;;;
   1 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
```
Num ; | Count | Num UTRs: ((Num ; + 1) * count)
--- | --- | ---
1	| 934 |	1868
2	| 104 |	312
3	| 59	| 236
4	| 43	| 215
5	| 29	| 174
6	| 13	| 91
7	| 8	| 64
8	| 12	| 108
9	| 3	| 30
10	| 1	| 11
11	| 3	| 36
12	| 2	| 26
14	| 1	| 15
15	| 2	| 32
18	| 1	| 19
24	| 1	| 25
26	| 1	| 27
34	| 1	| 35
**Total** | **1218**	| **3324**

diff=2106

Problem: 2 extra lines?? (2108)

# Format UTR data for TargetScan
```shell
$ awk '{ print $1, "7955", $6}' OFS='\t' fix_pos-output.txt > ts_input.txt
```

# Run TargetScan 6.0
```shell
$ perl ../targetscan_60_script/targetscan_60.pl miR_input2.txt ts_input.txt ts_output.txt
```

## Join TargetScan output with genomic position data
```shell
$ join -j1 -o1.1,1.2,1.4,1.5,1.9,2.2,2.3,2.4,2.5 -t $'\t' <(sort -k1 ts_output.txt) <(sort -k1 fix_pos-output.txt) > enrichment_input.txt
```

## Generate background file for GREAT
```shell
$ awk '{print $2, $4, $5, $1, 0, $3}' OFS='\t' fix_pos-output.txt > bg.bed6
```

# enrichment.R
Perform functional enrichment using GREAT
```R
R enrichment.R::ts.to.great() miR_input2.txt enrichment_input2.txt
```

# Visualize with heatmap.R
