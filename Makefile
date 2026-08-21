.PHONY: setup data datasets validate experiments plots genomics table stage repro-small verify all

setup:
	./setup.sh

data datasets:
	.venv/bin/python datasets/manage.py download

validate:
	.venv/bin/python datasets/manage.py validate

experiments:
	.venv/bin/python theory_validation/run_experiments.py
	.venv/bin/python method_comparison/run_experiments.py

plots:
	.venv/bin/python scripts/reproduce_plots.py

genomics:
	./run_genomics.sh

table: genomics

stage:
	.venv/bin/python datasets/manage.py stage

repro-small:
	./scripts/repro_small.sh

verify:
	.venv/bin/python scripts/test_frequency_counter.py
	.venv/bin/python datasets/manage.py validate
	.venv/bin/python scripts/completion_check.py --mark-complete

all: experiments plots genomics stage repro-small verify
