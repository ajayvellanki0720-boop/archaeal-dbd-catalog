.PHONY: install verify panel-sensitivity all test

install:
	python -m pip install -e '.[test]'

verify:
	python -m archaeal_dbd.cli audit

panel-sensitivity:
	python -m archaeal_dbd.cli panel-sensitivity

all:
	python -m archaeal_dbd.cli all

test:
	pytest -q

