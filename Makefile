predict_tags: tag_predictor.py
	python $< --auto

test: moieties.py
	@- rm patterns/*.png
	python $<
