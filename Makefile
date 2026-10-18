.PHONY: install verify panel-sensitivity independent-projection motif-holdout all test

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

all:
	python -m archaeal_dbd.cli all

test:
	pytest -q
