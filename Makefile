PACKAGE = $(notdir $(CURDIR))
VERSION = 1.0.0
SRC = DESCRIPTION README.Rmd $(wildcard R/*.R) $(wildcard vignettes/*.Rmd)
R = R
RARGS = --no-init-file
RSCRIPT = Rscript


.PHONY: build .build install .install check .check clean src


all: clean build install 

# Need submake to deal with Emacs
check checkascran build install test: 
	$(MAKE) EMACS="" .$@


.build: $(PACKAGE)_${VERSION}.tar.gz README.md


$(PACKAGE)_$(VERSION).tar.gz: $(SRC) NAMESPACE man
	$(R) $(RARGS) CMD build $(CURDIR)


src: $(PACKAGE)_$(VERSION)_src.tar.gz


$(PACKAGE)_$(VERSION)_src.tar.gz: $(SRC) NAMESPACE man
	git archive --format=tar.gz -o $@ -v HEAD


.install: $(PACKAGE)_$(VERSION).tar.gz
	$(R) $(RARGS) CMD INSTALL $<


.test:
	$(RSCRIPT) $(RARGS) -e "devtools::test()"


.check:
	$(RSCRIPT) $(RARGS) -e "devtools::check()"


.checkascran:
	$(RSCRIPT) $(RARGS) -e "devtools::check(cran=TRUE)"


NAMESPACE man: $(SRC)
	$(RSCRIPT) $(RARGS) -e "devtools::document()"


README.md: README.Rmd
	$(RSCRIPT) $(RARGS) -e "knitr::knit('$<')"


QFASA.pdf: $(SRC)
	$(R) $(RARGS) CMD Rd2pdf --output=$@ .


clean:
	$(RM) $(PACKAGE)_*.tar.gz
	$(RM) NAMESPACE
	$(RM) -rf man
	$(RM) README.md
	$(RM) *.pdf


