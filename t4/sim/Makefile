all: 1 2 

#ALINEA 1
1:ngspice-1.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	ps2pdf vo1f.ps vo1f.pdf
	ps2pdf vo2f.ps vo2f.pdf
	
	
ngspice-1.log: ../sim/t4.net
	ngspice -b $< -o $@
		
		
		
#ALINEA 2		
		
		
2:ngspice-2.log
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	
	
		
ngspice-2.log: ../sim/t42.net
	ngspice -b $< -o $@	
		
					
clean:
	@rm -f *.log *.tex *.pdf *.ps *~

.PHONY: all clean



