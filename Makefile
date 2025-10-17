CFLAGS = -Wall -pedantic -O3 -march=native -mtune=native -Iinclude
LDFLAGS = -lflint

STRAT = bruhat lu nonsing rej schubert subgroup unrank
SRC = $(patsubst %,src/test/t-%.c,$(STRAT))
OBJ = $(SRC:.c=.o)

COMB_SRC = $(wildcard src/comb/*.c)
UTIL_SRC = $(wildcard src/util/*.c)
COMB_OBJ = $(COMB_SRC:.c=.o)
UTIL_OBJ = $(UTIL_SRC:.c=.o)
TARGETS = $(basename $(SRC))

define INNER =
ifeq ($(1),unrank)
	EXTRA_OBJ_$(1) = $$(COMB_OBJ) $$(UTIL_OBJ)
endif
src/test/$(2)$(3)-$(1): src/test/$(2)$(3)-$(1).o src/$(1).o $$(EXTRA_OBJ_$(1)) \
	src/util/fq_nmod_mat_extra.o src/test/$(3)stub.o src/test/$(2)stub.o

ifeq ($(3),c)
$(2)$(3)-$(1): src/test/$(2)$(3)-$(1) src/test/rndcnt.so
	@LD_PRELOAD=src/test/rndcnt.so ./$$<
	@echo "[ OK ] $(1)"
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

src/test/rndcnt.so: src/test/rndcnt.c
	$(CC) -shared -fPIC -o $@ $<

$(eval $(call OUTER,test,t,e))
$(eval $(call OUTER,test-c,t,c))
$(eval $(call OUTER,prof,p,e))
$(eval $(call OUTER,prof-c,p,c))

clean: clean-te clean-tc clean-pe clean-pc
	$(RM) $(TARGETS) $(wildcard src/*.o) $(wildcard src/**/*.o) \
		src/test/rndcnt.so
