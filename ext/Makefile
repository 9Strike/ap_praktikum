cwd = $(shell pwd)
workspace_path = $(cwd)/..

module_name = sigval
ext_path = $(cwd)/$(module_name)
out_path = $(ext_path)/build

target: sigval.so
	cp $(shell find $(out_path) -name '*.so') $(workspace_path)/$(module_name).so
	rm -r $(out_path)

sigval.so: $(ext_path)/src/sigval.c
	cd $(ext_path); python3 setup.py build; cd $(cwd)

clean:
	rm $(workspace_path)/$(module_name).so