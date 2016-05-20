predict_tags: tag_predictor.py
	python $<

test: moieties.py
	@- rm patterns/*.png
	python $<
