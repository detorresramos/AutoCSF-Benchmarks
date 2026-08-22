.PHONY: setup datasets validate bound-validation synthetic-comparison genomics-comparison plots synthetic-clean stage repro-small verify reproduce all

setup:
	./setup.sh

datasets:
	.venv/bin/python datasets/manage.py download

validate:
	.venv/bin/python datasets/manage.py validate

bound-validation:
	.venv/bin/python experiments/bound_validation/run.py
	.venv/bin/python experiments/bound_validation/plot.py

synthetic-comparison:
	.venv/bin/python experiments/synthetic_comparison/run.py
	.venv/bin/python experiments/synthetic_comparison/plot.py

genomics-comparison:
	./experiments/genomics_comparison/run.sh

plots:
	.venv/bin/python scripts/reproduce_plots.py

synthetic-clean:
	./scripts/repro_synthetic_clean.sh

stage:
	.venv/bin/python datasets/manage.py stage

repro-small:
	./scripts/repro_small.sh

verify:
	.venv/bin/python scripts/test_frequency_counter.py
	.venv/bin/python datasets/manage.py validate
	.venv/bin/python scripts/completion_check.py --mark-complete

reproduce: bound-validation synthetic-comparison genomics-comparison verify

all: reproduce
