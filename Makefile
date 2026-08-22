PYTHON := .venv/bin/python

.PHONY: datasets bound-validation synthetic-comparison genomics-comparison \
        verify reproduce repro-clean publish-datasets

# Genomics inputs: ~500 MB processed, ~10 GB working space. Not needed by the
# two synthetic experiments.
datasets:
	$(PYTHON) datasets/manage.py download
	$(PYTHON) datasets/manage.py validate

bound-validation:
	$(PYTHON) experiments/bound_validation/run.py
	$(PYTHON) experiments/bound_validation/plot.py

synthetic-comparison:
	$(PYTHON) experiments/synthetic_comparison/run.py
	$(PYTHON) experiments/synthetic_comparison/plot.py

genomics-comparison:
	./experiments/genomics_comparison/run.sh

reproduce: bound-validation synthetic-comparison genomics-comparison verify

# Is what is in results/ complete and self-consistent?
verify:
	$(PYTHON) datasets/manage.py validate
	$(PYTHON) scripts/verify.py

# Does it still come out the same? Recomputes everything in a clean container
# that never sees results/, diffs the fresh numbers against the committed ones,
# and only promotes them if they match. SCOPE=synthetic (default) or SCOPE=small.
repro-clean:
	./scripts/repro_clean.sh

# Stage the genomics tables into an upload folder and push them to the
# Hugging Face dataset repo named in datasets/sources.json.
publish-datasets:
	$(PYTHON) datasets/manage.py stage
	$(PYTHON) datasets/manage.py publish
