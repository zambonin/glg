CFLAGS = -Wall -pedantic -O3 -march=native -mtune=native -Iinclude \
		 -DCLK_SPEED=$(shell cat /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq)
LDFLAGS = -L. -lglg -lflint

COMB_OBJ = $(patsubst %.c,%.o,$(wildcard src/comb/*.c))
UTIL_OBJ = src/util/math.o src/util/fq_nmod_mat_extra.o
LIB_OBJ = $(COMB_OBJ) $(UTIL_OBJ)
LIB_TARGET = libglg.a

STRAT = bruhat nonsing rej schubert subgroup unrank
SRC = $(patsubst %,src/test/t-%.c,$(STRAT))
TARGETS = $(basename $(SRC))

define INNER =
src/test/$(2)$(3)-$(1): src/test/$(2)$(3)-$(1).o src/$(1).o $(LIB_TARGET) \
	src/test/$(3)stub.o src/test/$(2)stub.o

ifeq ($(3),c)
$(2)$(3)-$(1): src/test/$(2)$(3)-$(1) src/test/rndcnt.so
	@LD_PRELOAD=src/test/rndcnt.so ./$$< | sed 's/$$$$/, alg = $(1)/g'
else
$(2)$(3)-$(1): src/test/$(2)$(3)-$(1)
	@./$$<
	@echo "[ OK ] $(1)"
endif
endef

define OUTER =
src/test/$(2)$(3)-%.o: src/test/t-%.o
	cp $$< $$@

$$(foreach s,$$(STRAT),$$(eval $$(call INNER,$$(s),$(2),$(3))))

$(1): $$(addprefix $(2)$(3)-,$$(STRAT))

clean-$(2)$(3):
	$$(RM) $$(patsubst src/test/t-%,src/test/$(2)$(3)-%,$$(TARGETS))
endef

default:

$(LIB_TARGET): $(LIB_OBJ)
	$(AR) rcs $@ $^

src/test/rndcnt.so: src/test/rndcnt.c
	$(CC) -shared -fPIC -o $@ $<

src/util/enum: src/util/enum.o $(LIB_TARGET)

$(eval $(call OUTER,test,t,e))
$(eval $(call OUTER,test-c,t,c))
$(eval $(call OUTER,prof,p,e))
$(eval $(call OUTER,prof-c,p,c))

clean: clean-te clean-tc clean-pe clean-pc
	$(RM) $(TARGETS) $(wildcard src/*.o) $(wildcard src/**/*.o) \
		src/test/rndcnt.so src/util/enum $(LIB_TARGET)
