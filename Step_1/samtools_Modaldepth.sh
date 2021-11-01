file=$1

samtools stats $file > ${file}_temp

number=$(grep ^COV ${file}_temp | cut -f 2- | awk -v max=0 '{if($3>max){want=$2; max=$3}}END{print want} ' -)

rm ${file}_temp

echo $file $number > ${file}_modal_depth.txt
