bib="bib/biblio.bib"
echo "bibliography file: $bib"
for mdname in "$@"; do
    echo "converting $mdname"
    pdfname=`echo ${mdname%%.*md}`
    pandoc $mdname --bibliography $bib --filter pandoc-eqnos --filter pandoc-fignos -N -o pdf/$pdfname.pdf
done
