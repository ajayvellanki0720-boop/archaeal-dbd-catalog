.PHONY: install verify panel-sensitivity independent-projection motif-holdout genus-blocked-motif lineage-blocked-motif quality-sensitivity taxonomic-prevalence family-separability publication-figures smoke all test

install:
	python -m pip install -e '.[test]'

verify:
	python -m archaeal_dbd.cli audit

panel-sensitivity:
	python -m archaeal_dbd.cli panel-sensitivity

independent-projection:
	python -m archaeal_dbd.cli independent-projection

motif-holdout:
	python -m archaeal_dbd.cli motif-holdout

genus-blocked-motif:
	python -m archaeal_dbd.cli genus-blocked-motif

quality-sensitivity:
	python -m archaeal_dbd.cli quality-sensitivity

taxonomic-prevalence:
	python -m archaeal_dbd.cli taxonomic-prevalence

family-separability:
	python -m archaeal_dbd.cli family-separability

lineage-blocked-motif:
	python -m archaeal_dbd.cli lineage-blocked-motif

publication-figures:
	python -m archaeal_dbd.cli publication-figures

smoke:
	python -m archaeal_dbd.smoke

all:
	python -m archaeal_dbd.cli all

test:
	pytest -q
