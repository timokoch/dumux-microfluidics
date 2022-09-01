name = main

all:
	pdflatex $(name)

final:
	pdflatex $(name)
	bibtex $(name)
	pdflatex $(name)
	bibtex $(name)
	pdflatex $(name)

clean:
	rm -f *.bbl  *.log *.aux  *.blg  *.out *.spl *.fdb_latexmk *.pag *.fls *.gz
