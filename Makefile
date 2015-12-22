build = build

all: $(build)/Makefile
	cd $(build) && make 

$(build)/Makefile: ./test/CMakeLists.txt
	mkdir -p $(build) && cd $(build) && cmake ../test -G"Unix Makefiles"

clean:
	$(MAKE) -C $(build) clean
	rm -rf $(build)

tests:
	octave -q ./test/test_vector_example.m
	octave -q ./test/test_eigen_vector_example.m
	octave -q ./test/test_eigen_matrix_example.m
	octave -q ./test/test_types.m

