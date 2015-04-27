all: data paper.pdf

data:
	mkdir -p paper/figures/generated/
	cd data && python descr_stats.py ./ ../paper/figures/generated \
		| tee ../paper/results_def.tex

paper.pdf: paper/p.tex paper/results_def.tex
	$(MAKE) -C paper
	cp paper/p.pdf paper.pdf

clean:
	$(MAKE) -C paper clean
	rm -f paper.pdf
	rm -rf data/segmentation data/timeseries paper/figures/generated

.PHONY: data
