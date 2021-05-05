for FILE in *.dot; 
   do
      filename=`echo $FILE | cut -f 1 -d '.'`
      dot -Teps $FILE -o $filename.eps
      epstopdf $filename.eps
      rm -f $filename.eps
done
