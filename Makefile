PYTHONPATH = PYTHONPATH="/home/$$USER/ATB:/home/$$USER/ATB-Dependencies"

PYTHON = $(PYTHONPATH) python

predict_tags: tag_predictor.py
	$(PYTHON) $< --auto

test: moieties.py
	@- rm patterns/*.png
	$(PYTHON) $<

errors:
	/usr/local/python35/bin/pylint -E $$(find . -name '*.py')
.PHONY: errors

mypy:
	MYPYPATH="/home/$${USER}/ATB_ONE" /usr/local/python35/bin/mypy $$(find . -name '*.py' |  sed "s|^\./||")
.PHONY: mypy
