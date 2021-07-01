all: ErrorModelSimpleEquant.so CommonDiet.so
	# Comment here preserves the prior tab
ErrorModelSimpleEquant.so: ErrorModelSimpleEquant.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('ErrorModelSimpleEquant.cpp','-O0 -g')"
CommonDiet.so: CommonDiet.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('CommonDiet.cpp','-O0 -g')"

clean:
	rm -rf *o

