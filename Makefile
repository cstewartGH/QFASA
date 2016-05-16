PACKAGE = $(notdir $(CURDIR))
VERSION = 0.0.0.9000
SRC = DESCRIPTION $(wildcard R/*.R) $(wildcard vignettes/*.Rmd)
R = R
RARGS = --no-init-file
RSCRIPT = Rscript


.PHONY: build install clean 

all: clean build install 


build: $(PACKAGE)_${VERSION}.tar.gz


$(PACKAGE)_$(VERSION).tar.gz: $(SRC) NAMESPACE man
	$(R) $(RARGS) CMD build $(CURDIR)	


install: $(PACKAGE)_$(VERSION).tar.gz
	$(R) $(RARGS) CMD INSTALL $<


clean:
	$(RM) $(PACKAGE)_${VERSION}.tar.gz
	$(RM) NAMESPACE
	$(RM) -rf man


NAMESPACE man: $(SRC)
	$(RSCRIPT) $(RARGS) -e "devtools::document()"
