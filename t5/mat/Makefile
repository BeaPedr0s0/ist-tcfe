
all: octave.log
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)

octave.log: lab5.m
	octave $< > $@

clean:
	@rm -f octave.log octave-workspace *.tex *.eps *.cir  *~ 
	#@rm -f data*
	#@rm -f *_tab.tex result_octave.txt

.PHONY: clean




