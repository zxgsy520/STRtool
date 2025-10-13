# STRtool
Short Tandem Repeats prediction result processing tool
## Using help
```
export PATH=/Work/pipeline/software/Base/trf/v4.09.1/:$PATH
trf plasmid.genomic.fasta 2 7 7 80 10 50 500 -f -d -h
python trf2gff.py plasmid.genomic.fasta.2.7.7.80.10.50.500.dat --repeat_num 2 --percent_matches 100 --outfa --outfa all.trf.fa >plasmid.TRF.gff3
python stat_length_gc.py all.trf.fa >all.trf.len_gc.tsv
Rscript plot_length_bar.R  txt.xls --prefix stat_plasmid_length --xlab "Length(bp)" --maxx 500
cut_plasmid_repeat.py plasmid.genomic.fasta --gff plasmid.77mer-TRF.gff3 >plasmid.rm_TR.fasta 2>log
```
