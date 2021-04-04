all: 1 2 3 4 5

#ALINEA 1
1:ngspice-1.log
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)

ngspice-1.log: ../sim/sim21.net
	ngspice -b $< -o $@
		
		
		
#ALINEA 2		
		
		
2:ngspice-2.log
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
		
		
ngspice-2.log: ../sim/sim22.net
	ngspice -b $< -o $@	
		
		
#ALINEA 3

3:ngspice-3.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	

	
ngspice-3.log: ../sim/sim23.net
	ngspice -b $< -o $@	
	
	
#ALINEA 4
4:ngspice-4.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	


ngspice-4.log: ../sim/sim24.net
	ngspice -b $< -o $@	
	
	
#ALINEA 5
5:ngspice-5.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	


ngspice-5.log: ../sim/sim25.net
	ngspice -b $< -o $@	
	
		
		
			
clean:
	@rm -f *.log *.tex *.pdf *.ps *~

.PHONY: all clean
