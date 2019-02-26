while getopts ":i:t:o:h" o; do
    case "${o}" in

        i) input=${OPTARG}
            ;;
	t) code=${OPTARG}
            ;;
        o) output=${OPTARG}
            ;;
        h) echo " rename sequences with: gene_R/F_platename.ab1
			-i input folder containing .ab1 files.
			-t .csv formatted file exported from the excel sheet.
			-o output folder.
				"
               exit
           ;;
         esac
     done

mkdir $output
tr -d $'\r' < $code > $output/tmp.lst
cp $input/*.ab1 $output
cd $output
for i in *.ab1; 
	do
	a=$(echo $i | awk -F "_" '{print $2}' | sed 's/\.ab1//' );
	for m in $a;
	do echo $m;
	export b=$(echo $m | grep -w -f /dev/stdin tmp.lst | awk -F ";" '{print $2"_"}');
#	echo $m $b
	mv $i $b$i;
     done;
  done
rm tmp.lst

for e in *.ab1; do mv $e ${e%_*}".ab1"; done
