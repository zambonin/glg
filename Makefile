CFLAGS = -Wall -pedantic -O3 -march=native -mtune=native -Iinclude
LDFLAGS = -lflint

STRAT = lu nonsing rej subgroup unrank
SRC = $(patsubst %,src/test/t-%.c,$(STRAT))
OBJ = $(SRC:.c=.o)
TARGETS = $(basename $(SRC))

COMB_SRC = $(wildcard src/comb/*.c)
UTIL_SRC = $(wildcard src/util/*.c)
COMB_OBJ = $(COMB_SRC:.c=.o)
UTIL_OBJ = $(UTIL_SRC:.c=.o)

src/test/p-%.o: src/test/t-%.o
	cp $< $@

define COMP_TEST =
ifeq ($(1),unrank)
	EXTRA_OBJ_$(1) = $$(COMB_OBJ) $$(UTIL_OBJ)
endif
src/test/t-$(1): src/test/t-$(1).o src/$(1).o src/util/fq_nmod_mat_extra.o \
	$$(EXTRA_OBJ_$(1)) src/test/tstub.o
endef

define COMP_PROF =
ifeq ($(1),unrank)
	EXTRA_OBJ_$(1) = $$(COMB_OBJ) $$(UTIL_OBJ)
endif
src/test/p-$(1): src/test/p-$(1).o src/$(1).o src/util/fq_nmod_mat_extra.o \
	$$(EXTRA_OBJ_$(1)) src/test/pstub.o
endef

define RUN_TEST =
$(1): src/test/$(1)
	@./$$< >/dev/null
	@echo "[ OK ] $(1)"
endef

define RUN_PROF =
$(1): src/test/$(1)
	./$$<
	@echo "[ OK ] $(1)"
endef

test: $(addprefix t-,$(STRAT))

prof: $(addprefix p-,$(STRAT))

$(foreach s,$(STRAT),$(eval $(call COMP_TEST,$(s))))
$(foreach s,$(STRAT),$(eval $(call RUN_TEST,t-$(s))))
$(foreach s,$(STRAT),$(eval $(call COMP_PROF,$(s))))
$(foreach s,$(STRAT),$(eval $(call RUN_PROF,p-$(s))))

clean:
	$(RM) $(TARGETS) $(wildcard src/*.o) $(wildcard src/**/*.o) \
		$(patsubst src/test/t-%,src/test/p-%,$(TARGETS))
