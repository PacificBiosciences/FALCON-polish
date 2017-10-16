WHEELHOUSE?=wheelhouse
PIP=pip wheel --wheel-dir ${WHEELHOUSE}

default:
pytest:
	python -c 'import falcon_polish; print falcon_polish'
	py.test ${MY_TEST_FLAGS} --junit-xml=test.xml --doctest-modules falcon_polish/ utest/
pylint:
	pylint --errors-only falcon_polish
autopep8:
	autopep8 --max-line-length=120 -ir -j0 falcon_polish/
wheel:
	which pip
	${PIP} --no-deps .
	ls -larth ${WHEELHOUSE}
clean:
	rm -f *.xml
	find . -name '*.pyc' | xargs rm
