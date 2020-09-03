phylocom_manual.pdf: phylocom_manual.tex phylocom_manual.bib ecology.bst
	pdflatex phylocom_manual.tex
	bibtex phylocom_manual 
	pdflatex phylocom_manual.tex

clean:
	rm -f *.aux &> /dev/null
	rm -f *.log &> /dev/null
	rm -f *.out &> /dev/null
	rm -f *.bbl &> /dev/null
	rm -f *.blg &> /dev/null
	rm -f *.idx &> /dev/null
	rm -f *.ilg &> /dev/null
	rm -f *.ind &> /dev/null
	rm -f *.toc &> /dev/null

