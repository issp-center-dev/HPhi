
TEX_COMMAND = platex
BIBTEX_COMMAND = bibtex
DVIPDF_COMMAND = dvipdfmx

TEX_FILES = \
abst01_jp.tex  \
ackn01_jp.tex  \
ap01_jp.tex  \
bib01_jp.tex  \
chap01_jp.tex  \
chap02_jp.tex  \
chap02_1_jp.tex  \
chap03_jp.tex  \
chap03_2_jp.tex  \
chap04_jp.tex  \
chap04_1_jp.tex  \
chap05_jp.tex \
chap06_jp.tex \
userguide_jp.tex

userguide_jp.pdf:userguide_jp.dvi
	dvipdfmx userguide_jp.dvi

userguide_jp.dvi:userguide_jp.tex $(TEX_FILES)
	platex userguide_jp.tex
	bibtex userguide_jp.aux
	platex userguide_jp.tex; platex userguide_jp.tex; platex userguide_jp.tex

clean:
	rm -f *.aux *.bbl *.blg *.dvi *.log *.pdf *.toc 

