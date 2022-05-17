# ----- Name ------
PROJECTNAME := gcmpy

# ----- Directories -----

SRCPATH := $(CURDIR)/gcmpy
BUILDPATH := $(CURDIR)/build
DIST_DIR := $(CURDIR)/dist

# ----- Tools -----

PYTHON = python3
TWINE = twine

# ----- Commands -----

.PHONY: requirements
requirements:
	pip3 install -U -r $(CURDIR)/requirements.txt

.PHONY: lint
lint:
	flake8 $(SRCPATH) --tee --output $(BUILDPATH)/lint.txt

.PHONY: format
format:
	black $(SRCPATH)

.PHONY: pre-commit
pre-commit: lint format

.PHONY: install
install:
	python3 -m pip install --editable .

.PHONY: upload
upload: 
	$(PYTHON) setup.py sdist
	$(TWINE) upload dist/* 

.PHONY: clean-distribution
clean-distribution:
	rm -rf $(DIST_DIR)

.PHONY: clean-build
clean-build:
	rm -rf $(BUILDPATH)

.PHONY: clean-egg-pyc
clean-egg-pyc:
	rm -rf $(SRCPATH)/*.egg-info
	find $(SRCPATH) -name "*.pyc" -delete
	find $(SRCPATH) -type d -name __pycache__ -delete

.PHONY: clean-all
clean-all:
	clean-distribution
	clean-build
	clean-egg-pyc

.PHONY: run-unit-tests
discover-unit-tests:
	$(PYTHON) -m unittest discover -v -s test/ -p 'test_*.py'

.PHONY: coverage-report
coverage-report:
	coverage run --source gcmpy -m unittest discover && coverage report

# ----- Usage -----

.PHONY: help
help:
	@make usage

define HELP_MESSAGE
Available targets:
   
   make lint 		 			run flake8 linter and output to build/lint.txt
   make format       			run black formatter
   make install 	 			installs package using setup.py
   make pre-commit   			checks linting and formats code
   make requirements 			installs the requirements
   make upload					makes and uploads the distribution to Pypi
   make clean-distribution		removes dist folder and contents
   make clean-build  			removes build folder and contents
   make clean-egg-pyc 			removes egg and pycache files from source directory
   make clean-all				runs clean-build, clean-egg-pyc and clean distribution
   make run-unit-tests     		runs unit tests with unittest
   make coverage-report			runs a coverage report with unittest
endef
export HELP_MESSAGE

usage:
	@echo "$$HELP_MESSAGE"