PYTHONPATH = PYTHONPATH="/home/$$USER/ATB:/home/$$USER/ATB-Dependencies"

PYTHON = $(PYTHONPATH) python

predict_tags: tag_predictor.py
	$(PYTHON) $< --auto

test: moieties.py
	@- rm patterns/*.png
	$(PYTHON) $<
