build = build

all: $(build)/Makefile
	cd $(build) && make 

$(build)/Makefile: CMakeLists.txt
	mkdir -p $(build) && cd $(build) && cmake .. -G"Unix Makefiles"

clean:
	$(MAKE) -C $(build) clean
	rm -rf $(build)

tests:
	octave -q ./test/test_vector_example.m
	octave -q ./test/test_eigen_vector_example.m
	octave -q ./test/test_eigen_matrix_example.m
	octave -q ./test/test_types.m

