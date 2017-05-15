#!/bin/sh


#trimm sequences
gunzip -c ZL097_S1_R1_001.fastq.gz | fastx_trimmer -l 32 -Q33 -o ZL097_1trimmed.fastq
gunzip -c ZL097_S1_R2_001.fastq.gz | fastx_trimmer -l 20 -Q33 -o ZL097_2trimmed.fastq
mkdir barcodesplitterqual

#qualityfilter data
paste -d '' ZL097_1trimmed.fastq ZL097_2trimmed.fastq | fastq_quality_filter -q 30 -p 90 -Q33 | awk "NR%4==2" | nl |awk '{print ">" $1 "\n" $2}'|fastx_barcode_splitter.pl --bcfile barcodefile.txt --prefix ~/ZL097/barcodesplitterqual/ --eol

cd barcodesplitterqual

for i in {97..132} ; do
#filter out reads with Ns, cut off indexes and unique datafiles
awk "NR%2==0" BC${i} | grep -v N | cut -b 1-44 |sort | uniq -c | sort -nr > ZL097processedBC${i}.txt
#split output files into two files per index, one that is containing the counts of each unique sequnce, the other the unique sequences themselves.
awk '{print $1}' ZL097processedBC${i}.txt > ZL097BC${i}_counts.txt
awk '{print $2}' ZL097processedBC${i}.txt > ZL097_BC${i}seq.txt
done


mkdir thresholds
cd thresholds

#pick thresholds from matlab plots, based on a steep drop of the 32+12 read counts. avoids too much crap in the data themselves. note, bash seems to start counting at 0, so put a random number at the first position of this array, to get the order right.
threshold=2
#do a collpase around the varital tags
for i in {97..132}; do
#j=${threshold[$i]}
(j=${threshold}
echo $j

awk '$1< '$j' {print NR}' ../ZL097BC${i}_counts.txt | head -n 1 >t
#grep -nr ^${threshold[$i]}$  ../ZL097BC${i}_counts.txt -m 1 | cut -f1 -d":" > t
thresh=$(cat t)
echo $thresh

head ../ZL097_BC${i}seq.txt -n $thresh | cut -b 1-32 | sort | uniq -c | sort -nr > ZL097${i}quickout.txt;)


done
wait
mkdir indexes

#ok, correct whats left
for i in {97..132}
do
(echo $i
in=ZL097${i}quickout.txt
#split off real barcodes from spikeins for easier analysis

grep -v 'ATCAGTCA$' $in | grep '[CT][CT]$' > ZL097BC${i}_quickprocessed.txt
awk '{print $1}' ZL097BC${i}_quickprocessed.txt > ZL097${i}_counts.txt
awk '{print $2}' ZL097BC${i}_quickprocessed.txt > ZL097${i}_seq.txt


nl ZL097${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ZL097_BC${i}fasta2u.txt; bowtie-build -q ZL097_BC${i}fasta2u.txt indexes/BC${i}fasta2u; bowtie -v 3 -p 10 -f --best -a indexes/BC${i}fasta2u ZL097_BC${i}fasta2u.txt bowtiealignment${i}_2u.txt
awk '{print $1}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_1.txt;awk '{print $3}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_3.txt;)&
done

wait



# deal with spike in sequences

for i in {97..132}; do 
(echo $i; in=ZL097${i}quickout.txt 
grep 'ATCAGTCA$' $in > ZL097spikes${i}_quickprocessed.txt; 
awk '{print $1}' ZL097spikes${i}_quickprocessed.txt > ZL097spikes${i}_counts.txt; 
awk '{print $2}' ZL097spikes${i}_quickprocessed.txt > ZL097spikes${i}_seq.txt;  
nl ZL097spikes${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ZL097_spikes${i}fasta2u.txt; 
bowtie-build -q ZL097_spikes${i}fasta2u.txt indexes/spikes${i}fasta2u; bowtie -v 3 -p 10 -f --best -a indexes/spikes${i}fasta2u ZL097_spikes${i}fasta2u.txt bowtiealignmentspikes${i}_2u.txt; 
awk '{print $1}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_1.txt;
awk '{print $3}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_3.txt;)&
done
wait