CFLAGS = -Wall -pedantic -O3 -march=native -mtune=native -Iinclude
LDFLAGS = -lflint

SRC = src/test/t-unrank.c
COMB_SRC = $(wildcard src/comb/*.c)
UTIL_SRC = $(wildcard src/util/*.c)

TARGET = $(basename $(SRC))
COMB_OBJ = $(COMB_SRC:.c=.o)
UTIL_OBJ = $(UTIL_SRC:.c=.o)

test: $(TARGET)
	@./$< >/dev/null
	@echo "[ OK ] $<"

$(TARGET): src/test/t-unrank.o $(COMB_OBJ) $(UTIL_OBJ)

clean:
	$(RM) $(TARGET) $(wildcard src/*.o) $(wildcard src/**/*.o)
