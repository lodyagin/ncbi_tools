#! /bin/sh

# usage ./rastbatch.sh template directory locus_tag_prefix product_id_prefix

PATH=$PATH:.
export PATH

find $2 -name "*.gb" | while read file; do
  ./rast2sqn.sh $1 $file $3 $4
done
