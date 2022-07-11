#!/bin/bash
istr="/path/to"
ostr="$PWD"

cat Circall_source.sh > Circall.sh
chmod 777 Circall.sh
eval "sed -i -e 's#"$istr"#"$ostr"#g'  Circall.sh R/load_packages.R R/getPseudoCircRNAsequence.R R/doPairEndFiltering.R R/doSingleEndFiltering.R R/getFdr.R R/export.R"

#chmod -R 777 ./linux
mv Circall.sh linux/bin/
echo "Done!"

