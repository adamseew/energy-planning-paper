all: sim_eq.pdf iros-2021.pdf

clean:
	$(RM) *~ *.pdf *.dvi *.log *.aux *.bbl *.blg *.toc *.lol *.loa *.lox \
		*.lot *.out *.lg *.tmp *.xref *.lof .*.swp *.synctex.gz

sim_eq.pdf: sim_eq.tex
	pdflatex sim_eq
	pdflatex sim_eq
	pdflatex sim_eq

iros-2021.pdf: iros-2021.tex
	pdflatex iros-2021
	bibtex iros-2021
	pdflatex iros-2021
	pdflatex iros-2021


