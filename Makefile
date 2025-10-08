CFLAGS = -Wall -pedantic -O3 -march=native -mtune=native -Iinclude
LDFLAGS = -lflint

STRAT = unrank
SRC = $(patsubst %,src/test/t-%.c,$(STRAT))
OBJ = $(SRC:.c=.o)
TARGETS = $(basename $(SRC))

COMB_SRC = $(wildcard src/comb/*.c)
UTIL_SRC = $(wildcard src/util/*.c)
COMB_OBJ = $(COMB_SRC:.c=.o)
UTIL_OBJ = $(UTIL_SRC:.c=.o)

define COMP_TEST =
src/test/t-$(1): src/test/t-$(1).o src/$(1).o $$(COMB_OBJ) $$(UTIL_OBJ) \
	src/test/tstub.o
endef

define RUN_TEST =
$(1): src/test/$(1)
	@./$$< >/dev/null
	@echo "[ OK ] $(1)"
endef

test: $(addprefix t-,$(STRAT))

$(foreach s,$(STRAT),$(eval $(call COMP_TEST,$(s))))
$(foreach s,$(STRAT),$(eval $(call RUN_TEST,t-$(s))))

clean:
	$(RM) $(TARGETS) $(wildcard src/*.o) $(wildcard src/**/*.o)
