all: sim_eq.pdf

clean:
	$(RM) *~ *.pdf *.dvi *.log *.aux *.bbl *.blg *.toc *.lol *.loa *.lox \
		*.lot *.out *.lg *.tmp *.xref *.lof .*.swp *.synctex.gz

sim_eq.pdf: sim_eq.tex
	pdflatex sim_eq
	bibtex sim_eq
	pdflatex sim_eq
	pdflatex sim_eq

