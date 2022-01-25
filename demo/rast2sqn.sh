#! /bin/sh

# usage ./rast2sqn.sh template filename locus_tag_prefix product_id_prefix

PATH=$PATH:.
export PATH

x=`dirname "$2"`/`basename "$2" ".gb"`

gbf2tbl.pl "$x".gb

mv "$x".tbl "$x".tmp

tblfix.pl -product_to_codon_recognized < "$x".tmp | \
tblfix.pl -revcomp_codon_recognized | \
tblfix.pl -split_compound_product | \
tblfix.pl -product_to_tc_number | \
tblfix.pl -product_to_ec_number | \
tblfix.pl -clean_rrna_product | \
tblfix.pl -create_locus_tag $3 | \
tblfix.pl -create_protein_id $4 > "$x".tbl

rm "$x".tmp

tbl2asn -t $1 -i "$x".fsa -U -a s6 -j "[gcode=11]" -R -T -c f -V v -Z "$x".dsc

