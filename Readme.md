# peamath

Plot graphs of functions right in your terminal. Works a bit like
[peaclock](https://github.com/octobanana/peaclock): the picture is drawn
in the middle of the terminal, `(0; 0)` is the centre of the screen, and
everything scales with the terminal size.

```
peamath -g "y=sin(x)"
peamath -g "x**2+y**2=1"
```

Graphs are drawn with Unicode braille characters (2x4 dots per cell), so
curves are smooth and a circle really looks like a circle.

## Build and install

Requirements: `gcc`, `make`, a UTF-8 terminal (every modern one).

```
make
sudo make install       # installs to /usr/local/bin
sudo make uninstall
make clean
```

Install somewhere else with `make install PREFIX=$HOME/.local`
(`DESTDIR` is supported too, for packagers).

## Usage

```
peamath -g EXPR [-g EXPR ...] [options]
peamath EXPR
```

| Option             | Meaning                                            |
|--------------------|----------------------------------------------------|
| `-g, --graph EXPR` | function or equation to plot (up to 8, repeatable) |
| `-z, --zoom N`     | initial zoom factor, default `1`                   |
| `-l, --legend`     | print the expressions in the top-left corner       |
| `-a, --ascii`      | plain ASCII instead of braille                     |
| `-v, --version`    | print the version                                  |
| `-h, --help`       | print help (also when run without arguments)       |

Keys: arrows or `h j k l` move the view, `+` / `-` zoom, `0` resets,
`q`, `Esc` or `Ctrl-C` quit.

If the output is not a terminal (a pipe or a file) peamath prints one frame
and exits, using `$COLUMNS` x `$LINES` as the size (80x24 by default).

### What can be plotted

* `y = f(x)`, e.g. `y=sin(x)`, `y=ln(x)`, `y=1/x`. The `y=` can be left out:
  `peamath "sin(x)"` is the same thing.
* any equation in `x` and `y`: `x**2+y**2=1`, `x=sin(y)`, `x**3-3*x=y**2`,
  `abs(x)+abs(y)=2`. An expression with `y` and no `=` means "`... = 0`".
* several graphs at once, each in its own colour:
  `peamath -g "y=ln(x)" -g "y=exp(x)" -l`

### Expressions

* operators: `+ - * / % ^ **` and parentheses; `^` and `**` are the same
  and are right-associative, so `2^3^2 = 512`
* implicit multiplication: `2x`, `3sin(x)`, `x(x+1)`, `2pi`
* variables: `x`, `y`
* constants: `pi` (or the letter π), `e`, `tau`
* functions:
  `sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh sec csc cot`,
  `exp ln log log2 log10 sqrt cbrt abs floor ceil round trunc sign fract`
  (`log` is the decimal logarithm, `ln` the natural one; `tg ctg arcsin
  arccos arctan arctg lg sgn` also work),
  `atan2(y, x) pow(a, b) min(a, b) max(a, b) mod(a, b) hypot(a, b)
  logb(base, x)`

Parentheses are always the safest way to write function arguments
(`sin x^2` is read as `sin(x^2)`, `sin 2x` as `sin(2)*x`). Remember to quote
the expression in the shell, otherwise `*` and `^` may be eaten by it.

## Notes

* Explicit graphs are sampled several times per dot column and connected, so
  steep functions and asymptotes (`tan`, `1/x`) look right.
* Equations are plotted by looking for sign changes of `lhs - rhs` on a grid
  of dots. This is fast and works for almost everything, but it cannot show
  an isolated point of a curve that never changes sign (such as `x^2+y^2=0`) and an equation with poles inside it
  can produce a spurious vertical line. Write such things as `y=f(x)`.
* The screen is redrawn only after a key press or a resize, so peamath does
  not use CPU while it just sits there.
* `NO_COLOR` is respected.
* Colours follow your terminal's colour scheme. On start peamath asks the
  terminal for its palette (OSC 4/10/11) and picks, from the scheme's own
  colours, the ones that stay readable on its background, so graphs look right
  on dark and light themes alike. Axes use a dimmed foreground colour (true
  colour if `COLORTERM` says it is supported). If the terminal does not answer
  (about 0.4 s wait at most), the built-in ANSI colours are used as before.
