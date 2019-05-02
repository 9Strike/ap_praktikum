first: sigval

sigval:
	make -C ext/sigval py_sigval
	cp ext/sigval/tmp/sigval.so sigval.so
