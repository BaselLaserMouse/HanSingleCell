#!/bin/sh


#trimm sequences
gunzip -c ZL102_S1_R1_001.fastq.gz | fastx_trimmer -l 32 -Q33 -o ZL102_1trimmed.fastq
gunzip -c ZL102_S1_R2_001.fastq.gz | fastx_trimmer -l 20 -Q33 -o ZL102_2trimmed.fastq
mkdir barcodesplitterqual

#qualityfilter data
paste -d '' ZL102_1trimmed.fastq ZL102_2trimmed.fastq | fastq_quality_filter -q 30 -p 90 -Q33 | awk "NR%4==2" | nl |awk '{print ">" $1 "\n" $2}'|fastx_barcode_splitter.pl --bcfile barcodefile.txt --prefix ~/ZL102/barcodesplitterqual/ --eol

cd barcodesplitterqual

for i in {133..168} ; do
#filter out reads with Ns, cut off indexes and unique datafiles
awk "NR%2==0" BC${i} | grep -v N | cut -b 1-44 |sort | uniq -c | sort -nr > ZL102processedBC${i}.txt
#split output files into two files per index, one that is containing the counts of each unique sequnce, the other the unique sequences themselves.
awk '{print $1}' ZL102processedBC${i}.txt > ZL102BC${i}_counts.txt
awk '{print $2}' ZL102processedBC${i}.txt > ZL102_BC${i}seq.txt
done


mkdir thresholds
cd thresholds

#pick thresholds from matlab plots, based on a steep drop of the 32+12 read counts. avoids too much crap in the data themselves. note, bash seems to start counting at 0, so put a random number at the first position of this array, to get the order right.
threshold=4
#do a collpase around the varital tags
for i in {133..168}; do
#j=${threshold[$i]}
(j=${threshold}
echo $j

awk '$1< '$j' {print NR}' ../ZL102BC${i}_counts.txt | head -n 1 >t
#grep -nr ^${threshold[$i]}$  ../ZL102BC${i}_counts.txt -m 1 | cut -f1 -d":" > t
thresh=$(cat t)
echo $thresh

head ../ZL102_BC${i}seq.txt -n $thresh | cut -b 1-32 | sort | uniq -c | sort -nr > ZL102${i}quickout.txt;)


done
wait
mkdir indexes

#ok, correct whats left
for i in {133..168}
do
(echo $i
in=ZL102${i}quickout.txt
#split off real barcodes from spikeins for easier analysis

grep -v 'ATCAGTCA$' $in | grep '[CT][CT]$' > ZL102BC${i}_quickprocessed.txt
awk '{print $1}' ZL102BC${i}_quickprocessed.txt > ZL102${i}_counts.txt
awk '{print $2}' ZL102BC${i}_quickprocessed.txt > ZL102${i}_seq.txt


nl ZL102${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ZL102_BC${i}fasta2u.txt; bowtie-build -q ZL102_BC${i}fasta2u.txt indexes/BC${i}fasta2u; bowtie -v 3 -p 10 -f --best -a indexes/BC${i}fasta2u ZL102_BC${i}fasta2u.txt bowtiealignment${i}_2u.txt
awk '{print $1}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_1.txt;awk '{print $3}' bowtiealignment${i}_2u.txt > bowtie${i}_2u_3.txt;)&
done

wait



# deal with spike in sequences

for i in {133..168}; do 
(echo $i; in=ZL102${i}quickout.txt 
grep 'ATCAGTCA$' $in > ZL102spikes${i}_quickprocessed.txt; 
awk '{print $1}' ZL102spikes${i}_quickprocessed.txt > ZL102spikes${i}_counts.txt; 
awk '{print $2}' ZL102spikes${i}_quickprocessed.txt > ZL102spikes${i}_seq.txt;  
nl ZL102spikes${i}_seq.txt | awk '{print ">" $1 "\n" $2}' > ZL102_spikes${i}fasta2u.txt; 
bowtie-build -q ZL102_spikes${i}fasta2u.txt indexes/spikes${i}fasta2u; bowtie -v 3 -p 10 -f --best -a indexes/spikes${i}fasta2u ZL102_spikes${i}fasta2u.txt bowtiealignmentspikes${i}_2u.txt; 
awk '{print $1}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_1.txt;
awk '{print $3}' bowtiealignmentspikes${i}_2u.txt > bowtiespikes${i}_2u_3.txt;)&
done
wait