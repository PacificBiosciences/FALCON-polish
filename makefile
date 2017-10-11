WHEELHOUSE?=wheelhouse
PIP=pip wheel --wheel-dir ${WHEELHOUSE}

default:
utest:
	nosetests -v utest/
wheel:
	which pip
	${PIP} --no-deps .
	ls -larth ${WHEELHOUSE}
.PHONY: utest
