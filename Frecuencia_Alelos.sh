#!/bin/bash

#echo -n 'Cual es la extension de los archivos? (csv,txt,dat,etc): '
#read extension

#consolida todos los archivos de entrada en una sola lista

for f in *.vcf
do

echo "procesando archivo $f ... "

grep -v "#" $f | awk '{ print $1 "_" $2 }' >> Raw.txt

done

#genera una nueva lista no redundante

#echo "empezando a consolidar"

sort Raw.txt | uniq -d > Duplicados.txt

sort Raw.txt | uniq -u > Unicos.txt

#echo "ya hice las listas"

cat Duplicados.txt Unicos.txt > Lista_Consolidada_Temp.txt

sort Lista_Consolidada_Temp.txt > Lista_Consolidada.txt

#echo "las juntamos y a eliminar los temporales"

rm Duplicados.txt Unicos.txt Raw.txt Lista_Consolidada_Temp.txt

#Crea listas considerando los criterios de clasificacion de variantes

for g in *.vcf
do

grep -v "#" $g | awk '{print $1 "\t" $2 "\t" $10}' | sed 's/:/\t/g' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' | sed 's/,/\t/g' |

awk 'FNR==NR { div = $4 / $6 ; div2 = $5 / $6 ; 

        if ($6 < 12) print $0 "\t" $6 "\tNA";

                else if ( $4 < $5 && div < 0.3 ) print $0 "\t" div "\tNA" ;
        
                        else if ( $4 < $5 && div >= 0.3 ) print $0 "\t" div "\tHetero"

                                else if ( $4 > $5 && div2 < 0.3 ) print $0 "\t" div2 "\tNA" ;

                                        else if ( $4 > $5 && div2 >= 0.3 ) print $0 "\t" div2 "\tHetero"

                                                else if ($4 == $5 && div == div2 ) print $0 "\t" div "\tHetero" ;

                                                        else if ($4 == $6 ) print $0 "\t" div "\tHomoRef" ;

                                                                else if ($5 == $6 ) print $0 "\t" div2 "\tHomoAlt" ;

                                                                         
}' > $g.nr.temp


awk '{print $1 "_" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' $g.nr.temp | sort | uniq -d >> $g.nr
awk '{print $1 "_" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' $g.nr.temp | sort | uniq -u >> $g.nr

rm $g.nr.temp

awk '
# process first file; use gene name and row number as indices for 2-dimensional array
FNR==NR { start[$1]
          next }

# process second file
{ for (x in start)
      # if gene name ($1) shows up as first dimension index in our array ...
      if ( $1==x )
         { print x "\t" $7
                 }
       
} ' Lista_Consolidada.txt $g.nr > $g.out.temp

rm $g.nr

cat $g.out.temp | 

awk 'FNR==NR {if ( $2 == "HomoAlt" ) print $1 "\t" 1;
	
         else if ( $2 == "Hetero" )  print $1 "\t" 2;

	 else if ( $2 == "NA" )  print $1 "\t" "NA";

         else		             print $1 "\t" 0; }

' > $g.out.temp2

sort $g.out.temp2 > $g.out

rm $g.out.temp $g.out.temp2

done

#Crea la tabla de salida

echo "Variant" > Tabla_temp.txt

awk '{print $1}' Lista_Consolidada.txt >> Tabla_temp.txt

for h in *.out
do

echo "$h" >> $h.casi

awk '{print $2}' $h >> $h.casi

#rm $h

done

paste Tabla_temp.txt *.casi > Tabla_Comparaciones.txt

rm *.casi Tabla_temp.txt

