PACKAGE = $(notdir $(CURDIR))
VERSION = 0.0.0.9000
SRC = DESCRIPTION $(wildcard R/*.R) $(wildcard vignettes/*.Rmd)
R = R
RARGS = --no-init-file
RSCRIPT = Rscript


.PHONY: build .build install .install check .check clean


all: clean build install 

# Need submake to deal with Emacs
check build install test: 
	$(MAKE) EMACS="" .$@


.build: $(PACKAGE)_${VERSION}.tar.gz


$(PACKAGE)_$(VERSION).tar.gz: $(SRC) NAMESPACE man
	$(R) $(RARGS) CMD build $(CURDIR)	


.install: $(PACKAGE)_$(VERSION).tar.gz
	$(R) $(RARGS) CMD INSTALL $<


.test:
	$(RSCRIPT) $(RARGS) -e "devtools::test()"


.check:
	$(RSCRIPT) $(RARGS) -e "devtools::check()"


NAMESPACE man: $(SRC)
	$(RSCRIPT) $(RARGS) -e "devtools::document()"


clean:
	$(RM) $(PACKAGE)_${VERSION}.tar.gz
	$(RM) NAMESPACE
	$(RM) -rf man
