all: sim_eq.pdf icra-2021.pdf

clean:
	$(RM) *~ *.pdf *.dvi *.log *.aux *.bbl *.blg *.toc *.lol *.loa *.lox \
		*.lot *.out *.lg *.tmp *.xref *.lof .*.swp *.synctex.gz

sim_eq.pdf: sim_eq.tex
	pdflatex sim_eq
	pdflatex sim_eq
	pdflatex sim_eq

icra-2021.pdf: icra-2021.tex
	pdflatex icra-2021
	bibtex icra-2021
	pdflatex icra-2021
	pdflatex icra-2021


