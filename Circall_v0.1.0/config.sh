#!/bin/bash
istr="/path/to"
ostr="$PWD"

eval "sed -i -e 's#"$istr"#"$ostr"#g'  Circall.sh R/load_packages.R R/getPseudoCircRNAsequence.R R/doPairEndFiltering.R R/doSingleEndFiltering.R R/getFdr.R R/export.R"

#chmod -R 777 ./linux

echo "Done!"

