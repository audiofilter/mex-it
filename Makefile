build = build

all: $(build)/Makefile
	cd $(build) && make 

$(build)/Makefile: CMakeLists.txt
	mkdir -p $(build) && cd $(build) && cmake .. -G"Unix Makefiles"

clean:
	$(MAKE) -C $(build) clean
	rm -rf $(build)
