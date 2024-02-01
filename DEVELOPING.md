# Development notes

Mostly notes to myself

# Unit Tests

```python -m unittest discover tests```

# CLI usage from repository root

```python -m pynrose.cli```

# Building docs
```commandline
cd docs
make clean
make html
firefox build/html/index.html
```

# Updating version

* 1.0.0 -> 2.0.0: ```hatch version major```
* 1.0.0 -> 1.1.0: ```hatch version minor```
* 1.0.0 -> 1.0.1: ```hatch version patch```

# Building dist and uploading

```commandline
hatch clean
hatch build
python -m twine upload dist/*
```

enter the literal value `__token__` as the username, and the API key as the password
