.PHONY: setup pipeline dashboard test

PYTHON := $(shell if command -v python >/dev/null 2>&1; then echo python; elif command -v python3.11 >/dev/null 2>&1; then echo python3.11; else echo python3; fi)

setup:
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analyze.py

dashboard:
	$(PYTHON) -m streamlit run app.py --server.address 0.0.0.0 --server.port 8501

test:
	$(PYTHON) -m unittest discover -s tests -v
