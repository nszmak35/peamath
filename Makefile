# peamath - plot math functions in the terminal

TARGET  = peamath
SRC     = src/peamath.c

CC      = gcc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99
LDLIBS  = -lm

PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $(SRC) $(LDFLAGS) $(LDLIBS)

install: $(TARGET)
	install -d $(DESTDIR)$(BINDIR)
	install -m 755 $(TARGET) $(DESTDIR)$(BINDIR)/$(TARGET)

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all install uninstall clean
