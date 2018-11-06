#!/bin/bash

mypp() {
       go=$1;
       of=$(mktemp "$pp"/tmpg2p.XXXXXX); #create a temporary file
       grep "$go" "$pp"/"$gaf" > "$of"; 
}
export -f mypp;

#GO term query
q1="GO:0009819
GO:0097207
GO:1902584
GO:2000070
GO:0042630
GO:0042631
GO:0009414";

t="At"; #taxon, At, Sor, Pop
pp=$(pwd);
gaf="$t""RawGAF.txt"; #the GAF2.1 file associating GO terms and gene names

#species specific variables, "reg" is empirically determined as the regex representing all possible gene names in the gff3 file.
#it is used to query the multiple columns that may contain a gene name in the GAF file.
#"csub" is a string attached to chromosome number in column 1 of the gff3
if     [ "$t" == "At" ];  then gff="TAIR10_GFF3_genes.gff";      reg="^AT[1-5,C,M]G[0-9][0-9][0-9][0-9][0-9]$";               csub="Chr";
  elif [ "$t" == "Sor" ]; then gff="Sbi1.4.gff3";                reg="^Sb[0-3][0-9][0-9,g][0-9][0-5,s][0-9][0-9][0-9][0-9]$"; csub="chromosome_";
  elif [ "$t" == "Pop" ]; then gff="Ptrichocarpa_156_gene.gff3"; reg="^POPTR_[0-3][0-9][0-9][0-9]s[0-9][0-9][0-9][0-9]$";     csub="scaffold_"
fi;
export gaf gff pp;

rm -f tmpg2p.*; #clean up any old temporary files
echo "$q1" | parallel --env gaf --env gff --env pp --env mypp mypp; #extract the association btw GO term and gene name from GAF

find "$pp" -name "tmpg2p*" -size 0 -exec rm {} \; #remove all results of size 0
q2=$(cat tmpg2p.* | cut -d$'\t' -f3,10,11 | sort -u | sed 's/|/'$'\t/g'); #gene name query, must include 3 columns that contain gene identifiers

#extract the gene names associated with GO terms in proper format for querying the gff
q3="";
>log.txt; #start a log file to record problems
while read -r ll;
  do d=$(echo "$ll" | tr "\t" "\n" | grep "$reg" | sort -u); #retrieve gene names that fit "reg" pattern
    if [[ $(echo "$d" | wc -l | awk '{print $1}') > 1 ]]; then echo "Warning. More than 1 gene name:"$'\n'"$d"$'\n'"in: $ll" >> log.txt; fi;
    if [ "$d" = "" ]; then echo "Warning. No genes were found for: $ll" >> log.txt; fi;
    q3+="$d ";
  done<<<"$q2";
q3=$(echo "$q3" | gsed 's/  */ /g' | sed 's/ $//' | tr " " "\n" | sort -u); #condense space runs to one space, eliminate trailing space, convert return to space from when multiple genes were found, eliminate repeated genes

#create a db of relevant lines from the gff3 file, e.g. all 'genes'.
s=$(awk -F$'\t' '$3 == "gene"' "$pp"/"$gff"); #extract all lines that specify genes

of="$t"go2pos.txt; #define the output file
>"$of";
for i in $q3;
  do b=$(grep "$i" <(echo "$s") | cut -d$'\t' -f1,4,5 | sed 's/'$csub'//g';) # extract the chromosome, gene start, gene end from gff3
    if [ "$b" = "" ]; then echo "Warning. Gene name $i not found in gff3 file." >> log.txt; continue; fi;
    chr=$(echo "$b" | cut -d$'\t' -f1); #chromosome
    gs=$(echo "$b" | cut -d$'\t' -f2); #gene start
    ge=$(echo "$b" | cut -d$'\t' -f3); #gene end
    echo "$i     $chr"."$gs":"$chr"."$ge"; #output gene locations to stdout in human readable format
    echo "$chr"."$gs":"$chr"."$ge" >> "$of"; #output gene locations in format for haplotypista input
  done;

#clean up
rm -f tmpg2p.*; #clean up any old temporary files




#SOME TOOLS FOR USING GO2POS"
#query the gff3 file for all possible gene names, determine a regex to represent it
q2=$(awk -F$'\t' '$3 == "gene"' "Ptrichocarpa_156_gene.gff3" | cut -d$'\t' -f9 | cut -d\; -f1 | sed 's/ID=//g');
for i in {1..16}; do echo -n "$i: "; echo "$q2" | cut -c$i | sort -u | tr "\n" ","; echo; done; #determine all possible characters to make up the regex

#pipe the results file to haplotypista
cd /Volumes/attica/Desktop/Attica-Pats\ Folder/Scripts/Haplotype\ analysis/haplotypista; 
head /Volumes/attica/Desktop/Attica-Pats\ Folder/Scripts/go2pos/Atgo2pos.txt | ./haplotypista -i AtExample.txt -o atout -b 1 4 -m ? -p 1 -v Atpopid.txt;
cat /Volumes/attica/Desktop/Attica-Pats\ Folder/Scripts/go2pos/Atgo2pos.txt | ./haplotypista -i AtExample.txt -o atout -b 1 4 -m ? -p 1 -v Atpopid.txt;

echo "5.26300478:5.26300478
5.26519242:5.26519242
5.26712971:5.26712971" | ./haplotypista -i AtExample.txt -o atout -l atlog.txt -b 1 4 -m ? -p 1 -v Atpopid.txt;


1.75000:1.1000000,14.8697509:14.8697509