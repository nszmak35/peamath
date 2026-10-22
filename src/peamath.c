/*
 * peamath - plot math functions in the terminal.
 * I LOVE MATH!!!
 */

#define _DEFAULT_SOURCE

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <poll.h>
#include <setjmp.h>
#include <signal.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/ioctl.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>

#define VERSION     "0.1"
#define MAXFUNCS    8
#define OVERSAMPLE  4

static void die(const char *msg)
{
	fprintf(stderr, "peamath: %s\n", msg);
	exit(1);
}

static void *xcalloc(size_t n, size_t sz)
{
	void *p = calloc(n ? n : 1, sz);
	if (!p)
		die("out of memory");
	return p;
}

static int clampi(int v, int lo, int hi)
{
	return v < lo ? lo : (v > hi ? hi : v);
}

static void write_all(int fd, const char *b, size_t n)
{
	while (n) {
		ssize_t w = write(fd, b, n);
		if (w < 0) {
			if (errno == EINTR)
				continue;
			return;
		}
		b += w;
		n -= (size_t)w;
	}
}

static char *obuf;
static size_t olen, ocap;

static void app_n(const char *s, size_t n)
{
	if (olen + n + 1 > ocap) {
		size_t nc = ocap ? ocap : 4096;
		while (nc < olen + n + 1)
			nc *= 2;
		obuf = realloc(obuf, nc);
		if (!obuf)
			die("out of memory");
		ocap = nc;
	}
	memcpy(obuf + olen, s, n);
	olen += n;
}

static void app_s(const char *s) { app_n(s, strlen(s)); }
static void app_c(char c)        { app_n(&c, 1); }

typedef double (*f1_t)(double);
typedef double (*f2_t)(double, double);

enum { K_X, K_Y, K_CONST, K_F1, K_F2 };

typedef struct {
	const char *name;
	int kind;
	double c;
	f1_t f1;
	f2_t f2;
} Sym;

static double f_sec(double x)   { return 1.0 / cos(x); }
static double f_csc(double x)   { return 1.0 / sin(x); }
static double f_cot(double x)   { return 1.0 / tan(x); }
static double f_ln(double x)    { return log(x); }
static double f_lg(double x)    { return log10(x); }
static double f_lb(double x)    { return log2(x); }
static double f_sign(double x)  { return x > 0 ? 1.0 : (x < 0 ? -1.0 : x); }
static double f_fract(double x) { return x - floor(x); }

static double g_pow(double a, double b)  { return pow(a, b); }
static double g_min(double a, double b)  { return a < b ? a : b; }
static double g_max(double a, double b)  { return a > b ? a : b; }
static double g_mod(double a, double b)  { return a - b * floor(a / b); }
static double g_logb(double b, double x) { return log(x) / log(b); }

#define F1(n, f) { n, K_F1, 0, f, NULL }
#define F2(n, f) { n, K_F2, 0, NULL, f }
#define CT(n, v) { n, K_CONST, v, NULL, NULL }

static const Sym syms[] = {
	{ "x", K_X, 0, NULL, NULL },
	{ "y", K_Y, 0, NULL, NULL },
	CT("pi", M_PI), CT("tau", 2 * M_PI), CT("e", M_E),

	F1("sin", sin), F1("cos", cos), F1("tan", tan), F1("tg", tan),
	F1("asin", asin), F1("acos", acos), F1("atan", atan),
	F1("arcsin", asin), F1("arccos", acos), F1("arctan", atan),
	F1("arctg", atan),
	F1("sinh", sinh), F1("cosh", cosh), F1("tanh", tanh),
	F1("asinh", asinh), F1("acosh", acosh), F1("atanh", atanh),
	F1("sec", f_sec), F1("csc", f_csc), F1("cot", f_cot), F1("ctg", f_cot),
	F1("exp", exp), F1("ln", f_ln), F1("log", f_lg), F1("lg", f_lg),
	F1("log2", f_lb), F1("log10", f_lg),
	F1("sqrt", sqrt), F1("cbrt", cbrt), F1("abs", fabs),
	F1("floor", floor), F1("ceil", ceil), F1("round", round),
	F1("trunc", trunc), F1("sign", f_sign), F1("sgn", f_sign),
	F1("fract", f_fract),

	F2("atan2", atan2), F2("pow", g_pow), F2("min", g_min),
	F2("max", g_max), F2("mod", g_mod), F2("hypot", hypot),
	F2("logb", g_logb),
};
#define NSYMS ((int)(sizeof syms / sizeof syms[0]))

static const Sym *find_sym(const char *name)
{
	for (int i = 0; i < NSYMS; i++)
		if (!strcmp(syms[i].name, name))
			return &syms[i];
	return NULL;
}

enum { N_NUM, N_X, N_Y, N_NEG, N_ADD, N_SUB, N_MUL, N_DIV, N_MOD, N_POW,
       N_F1, N_F2 };

typedef struct Node {
	int type;
	double v;
	f1_t f1;
	f2_t f2;
	struct Node *a, *b;
} Node;

static Node *mk(int type, Node *a, Node *b)
{
	Node *n = xcalloc(1, sizeof *n);
	n->type = type;
	n->a = a;
	n->b = b;
	return n;
}

static double eval(const Node *n, double x, double y)
{
	switch (n->type) {
	case N_NUM: return n->v;
	case N_X:   return x;
	case N_Y:   return y;
	case N_NEG: return -eval(n->a, x, y);
	case N_ADD: return eval(n->a, x, y) + eval(n->b, x, y);
	case N_SUB: return eval(n->a, x, y) - eval(n->b, x, y);
	case N_MUL: return eval(n->a, x, y) * eval(n->b, x, y);
	case N_DIV: return eval(n->a, x, y) / eval(n->b, x, y);
	case N_MOD: {
		double a = eval(n->a, x, y), b = eval(n->b, x, y);
		return a - b * floor(a / b);
	}
	case N_POW: return pow(eval(n->a, x, y), eval(n->b, x, y));
	case N_F1:  return n->f1(eval(n->a, x, y));
	case N_F2:  return n->f2(eval(n->a, x, y), eval(n->b, x, y));
	}
	return NAN;
}

static int has_var(const Node *n, int type)
{
	if (!n)
		return 0;
	if (n->type == type)
		return 1;
	return has_var(n->a, type) || has_var(n->b, type);
}

enum { T_END, T_NUM, T_SYM, T_PLUS, T_MINUS, T_MUL, T_DIV, T_MOD, T_POW,
       T_LP, T_RP, T_COMMA, T_EQ };

typedef struct {
	int type;
	int pos;
	double num;
	const Sym *sym;
} Tok;

typedef struct {
	const char *src;
	Tok *tok;
	int i;
	jmp_buf jb;
	char err[96];
	int errpos;
} Parser;

static void fail(Parser *p, int pos, const char *fmt, ...)
	__attribute__((noreturn));

static void fail(Parser *p, int pos, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vsnprintf(p->err, sizeof p->err, fmt, ap);
	va_end(ap);
	p->errpos = pos;
	longjmp(p->jb, 1);
}

static void lex(Parser *p)
{
	const char *s = p->src;
	size_t len = strlen(s), i = 0;
	int n = 0;

	p->tok = xcalloc(len + 2, sizeof(Tok));

	while (i < len) {
		unsigned char ch = (unsigned char)s[i];
		Tok t = { T_END, (int)i, 0, NULL };

		if (isspace(ch)) {
			i++;
			continue;
		}

		if (isdigit(ch) || (ch == '.' && isdigit((unsigned char)s[i + 1]))) {
			size_t j = i, k;
			char buf[64];

			while (isdigit((unsigned char)s[j]))
				j++;
			if (s[j] == '.') {
				j++;
				while (isdigit((unsigned char)s[j]))
					j++;
			}
			if (s[j] == 'e' || s[j] == 'E') {
				k = j + 1;
				if (s[k] == '+' || s[k] == '-')
					k++;
				if (isdigit((unsigned char)s[k])) {
					while (isdigit((unsigned char)s[k]))
						k++;
					j = k;
				}
			}
			if (j - i >= sizeof buf)
				fail(p, (int)i, "number is too long");
			memcpy(buf, s + i, j - i);
			buf[j - i] = 0;
			t.type = T_NUM;
			t.num = strtod(buf, NULL);
			i = j;
		} else if (isalpha(ch) || (ch == 0xCF && (unsigned char)s[i + 1] == 0x80)) {
			const Sym *best = NULL;
			size_t bl = 0;

			if (ch == 0xCF) {
				best = find_sym("pi");
				bl = 2;
			} else {
				for (int k = 0; k < NSYMS; k++) {
					size_t l = strlen(syms[k].name);
					if (l > bl && strncasecmp(s + i, syms[k].name, l) == 0) {
						best = &syms[k];
						bl = l;
					}
				}
			}
			if (!best) {
				char w[24];
				size_t j = i, l = 0;
				while (isalpha((unsigned char)s[j]) && l < sizeof w - 1)
					w[l++] = s[j++];
				w[l] = 0;
				fail(p, (int)i, "unknown name '%s'", w);
			}
			t.type = T_SYM;
			t.sym = best;
			i += bl;
		} else {
			switch (ch) {
			case '+': t.type = T_PLUS;  i++; break;
			case '-': t.type = T_MINUS; i++; break;
			case '/': t.type = T_DIV;   i++; break;
			case '%': t.type = T_MOD;   i++; break;
			case '^': t.type = T_POW;   i++; break;
			case '(': t.type = T_LP;    i++; break;
			case ')': t.type = T_RP;    i++; break;
			case ',': t.type = T_COMMA; i++; break;
			case '*':
				if (s[i + 1] == '*') {
					t.type = T_POW;
					i += 2;
				} else {
					t.type = T_MUL;
					i++;
				}
				break;
			case '=':
				t.type = T_EQ;
				i += (s[i + 1] == '=') ? 2 : 1;
				break;
			default:
				if (ch >= 0x80)
					fail(p, (int)i, "unexpected non-ASCII character");
				fail(p, (int)i, "unexpected character '%c'", ch);
			}
		}
		p->tok[n++] = t;
	}
	p->tok[n].type = T_END;
	p->tok[n].pos = (int)len;
}

static Tok *peek(Parser *p) { return &p->tok[p->i]; }

static void expect(Parser *p, int type, const char *what)
{
	Tok *t = peek(p);
	if (t->type != type)
		fail(p, t->pos, "expected %s", what);
	p->i++;
}

static Node *parse_expr(Parser *p);
static Node *parse_power(Parser *p);

static Node *parse_primary(Parser *p)
{
	Tok *t = peek(p);
	Node *n, *a, *b;

	switch (t->type) {
	case T_NUM:
		p->i++;
		n = mk(N_NUM, NULL, NULL);
		n->v = t->num;
		return n;
	case T_LP:
		p->i++;
		n = parse_expr(p);
		expect(p, T_RP, "')'");
		return n;
	case T_SYM: {
		const Sym *s = t->sym;
		p->i++;
		switch (s->kind) {
		case K_X:
			return mk(N_X, NULL, NULL);
		case K_Y:
			return mk(N_Y, NULL, NULL);
		case K_CONST:
			n = mk(N_NUM, NULL, NULL);
			n->v = s->c;
			return n;
		case K_F1:
			if (peek(p)->type == T_LP) {
				p->i++;
				a = parse_expr(p);
				if (peek(p)->type == T_COMMA)
					fail(p, peek(p)->pos, "%s takes one argument", s->name);
				expect(p, T_RP, "')'");
			} else {
				a = parse_power(p);
			}
			n = mk(N_F1, a, NULL);
			n->f1 = s->f1;
			return n;
		case K_F2:
			if (peek(p)->type != T_LP)
				fail(p, peek(p)->pos, "%s needs two arguments: %s(a, b)",
				     s->name, s->name);
			p->i++;
			a = parse_expr(p);
			expect(p, T_COMMA, "','");
			b = parse_expr(p);
			expect(p, T_RP, "')'");
			n = mk(N_F2, a, b);
			n->f2 = s->f2;
			return n;
		}
		break;
	}
	case T_END:
		fail(p, t->pos, "unexpected end of expression");
	case T_RP:
		fail(p, t->pos, "unmatched ')'");
	default:
		fail(p, t->pos, "unexpected token");
	}
	return NULL;
}

static Node *parse_unary(Parser *p);

static Node *parse_power(Parser *p)
{
	Node *b = parse_primary(p);
	if (peek(p)->type == T_POW) {
		p->i++;
		return mk(N_POW, b, parse_unary(p));
	}
	return b;
}

static Node *parse_unary(Parser *p)
{
	int t = peek(p)->type;
	if (t == T_MINUS) {
		p->i++;
		return mk(N_NEG, parse_unary(p), NULL);
	}
	if (t == T_PLUS) {
		p->i++;
		return parse_unary(p);
	}
	return parse_power(p);
}

static Node *parse_term(Parser *p)
{
	Node *l = parse_unary(p);
	for (;;) {
		int t = peek(p)->type;
		if (t == T_MUL) {
			p->i++;
			l = mk(N_MUL, l, parse_unary(p));
		} else if (t == T_DIV) {
			p->i++;
			l = mk(N_DIV, l, parse_unary(p));
		} else if (t == T_MOD) {
			p->i++;
			l = mk(N_MOD, l, parse_unary(p));
		} else if (t == T_NUM || t == T_SYM || t == T_LP) {
			l = mk(N_MUL, l, parse_unary(p));
		} else {
			break;
		}
	}
	return l;
}

static Node *parse_expr(Parser *p)
{
	Node *l = parse_term(p);
	for (;;) {
		int t = peek(p)->type;
		if (t == T_PLUS) {
			p->i++;
			l = mk(N_ADD, l, parse_term(p));
		} else if (t == T_MINUS) {
			p->i++;
			l = mk(N_SUB, l, parse_term(p));
		} else {
			break;
		}
	}
	return l;
}

static int parse_equation(const char *src, Node **lhs, Node **rhs,
                          char *err, size_t en, int *errpos)
{
	static Parser P;
	Parser *p = &P;

	memset(p, 0, sizeof *p);
	p->src = src;
	if (setjmp(p->jb)) {
		snprintf(err, en, "%s", p->err);
		*errpos = p->errpos;
		return -1;
	}
	lex(p);
	if (peek(p)->type == T_END)
		fail(p, 0, "empty expression");
	*lhs = parse_expr(p);
	*rhs = NULL;
	if (peek(p)->type == T_EQ) {
		p->i++;
		*rhs = parse_expr(p);
	}
	if (peek(p)->type != T_END) {
		if (peek(p)->type == T_RP)
			fail(p, peek(p)->pos, "unmatched ')'");
		fail(p, peek(p)->pos, "unexpected token");
	}
	return 0;
}

typedef struct {
	char *text;
	int implicit;
	Node *e;
} Func;

static Func funcs[MAXFUNCS];
static int nfuncs;

static double zoom0 = 1.0, zoom = 1.0, cx, cy;
static int opt_legend, opt_ascii, usecolor;

static int cols, rows, W, H;
static uint8_t *gm, *am, *gc, *lc;
static char *lb;
static size_t cells_alloc;

static const uint8_t BIT[4][2] = {
	{ 0x01, 0x08 }, { 0x02, 0x10 }, { 0x04, 0x20 }, { 0x40, 0x80 }
};

static void canvas_resize(int c, int r)
{
	size_t n;

	if (c < 1) c = 1;
	if (r < 1) r = 1;
	cols = c; rows = r; W = c * 2; H = r * 4;
	n = (size_t)c * (size_t)r;
	if (n > cells_alloc) {
		free(gm); free(am); free(gc); free(lc); free(lb);
		gm = xcalloc(n, 1); am = xcalloc(n, 1);
		gc = xcalloc(n, 1); lc = xcalloc(n, 1);
		lb = xcalloc(n, 1);
		cells_alloc = n;
	}
}

static void gdot(int px, int py, int fi)
{
	size_t cell;
	if (px < 0 || py < 0 || px >= W || py >= H)
		return;
	cell = (size_t)(py >> 2) * cols + (px >> 1);
	gm[cell] |= BIT[py & 3][px & 1];
	gc[cell] = (uint8_t)(fi + 1);
}

static void adot(int px, int py)
{
	if (px < 0 || py < 0 || px >= W || py >= H)
		return;
	am[(size_t)(py >> 2) * cols + (px >> 1)] |= BIT[py & 3][px & 1];
}

static void put_text(int r, int c0, const char *s, int col)
{
	if (r < 0 || r >= rows)
		return;
	for (; *s; s++, c0++) {
		if (c0 < 0 || c0 >= cols)
			continue;
		lb[(size_t)r * cols + c0] = *s;
		lc[(size_t)r * cols + c0] = (uint8_t)col;
	}
}

static double ppu_now(void)
{
	double b = fmin(H / 10.0, W / 16.0);
	if (b < 1.0)
		b = 1.0;
	return b * zoom;
}

static void draw_axes(double ppu)
{
	double fx = W / 2.0 - cx * ppu;
	double fy = H / 2.0 + cy * ppu;
	int xin = fx >= 0 && fx < W, yin = fy >= 0 && fy < H;
	int ax = xin ? (int)floor(fx) : 0, ay = yin ? (int)floor(fy) : 0;
	int ar = fy < 0 ? -1 : (fy >= H ? rows : (int)fy / 4);
	int ac = fx < 0 ? -1 : (fx >= W ? cols : (int)fx / 2);
	double raw = 16.0 / ppu, mag = pow(10.0, floor(log10(raw)));
	double m = raw / mag;
	double step = (m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10) * mag;
	char b[40];
	int k, k0, k1;

	if (yin)
		for (int px = 0; px < W; px++)
			adot(px, ay);
	if (xin)
		for (int py = 0; py < H; py++)
			adot(ax, py);

	k0 = (int)ceil((cx - W / 2.0 / ppu) / step);
	k1 = (int)floor((cx + W / 2.0 / ppu) / step);
	for (k = k0; k <= k1; k++) {
		double v;
		int px, l;
		if (k == 0)
			continue;
		v = k * step;
		px = (int)floor(W / 2.0 + (v - cx) * ppu);
		if (yin) {
			adot(px, ay - 1);
			adot(px, ay + 1);
		}
		snprintf(b, sizeof b, "%.6g", v);
		l = (int)strlen(b);
		put_text(clampi(ar + 1, 0, rows - 1), px / 2 - l / 2, b, 0);
	}

	k0 = (int)ceil((cy - H / 2.0 / ppu) / step);
	k1 = (int)floor((cy + H / 2.0 / ppu) / step);
	for (k = k0; k <= k1; k++) {
		double v;
		int py, l, col;
		if (k == 0)
			continue;
		v = k * step;
		py = (int)floor(H / 2.0 - (v - cy) * ppu);
		if (xin) {
			adot(ax - 1, py);
			adot(ax + 1, py);
		}
		snprintf(b, sizeof b, "%.6g", v);
		l = (int)strlen(b);
		col = ac + 1;
		if (col + l > cols)
			col = ac - l;
		put_text(clampi(py / 4, 0, rows - 1), col, b, 0);
	}

	if (xin && yin)
		put_text(clampi(ar + 1, 0, rows - 1), ac + 1, "0", 0);
}

static void seg(double x0, double y0, double x1, double y1, int fi)
{
	double dx, dy, hh = (double)H;
	int o0 = y0 < 0 || y0 >= hh, o1 = y1 < 0 || y1 >= hh, n;

	if (o0 && o1)
		return;
	if (o0 || o1) {
		double yb = o0 ? y0 : y1;
		double edge = yb < 0 ? 0.0 : hh - 1e-9;
		double t = (edge - y0) / (y1 - y0);
		double xe = x0 + t * (x1 - x0);
		if (o0) { x0 = xe; y0 = edge; }
		else    { x1 = xe; y1 = edge; }
	}
	dx = x1 - x0;
	dy = y1 - y0;
	n = (int)ceil(fmax(fabs(dx), fabs(dy))) + 1;
	for (int i = 0; i <= n; i++) {
		double t = (double)i / n;
		gdot((int)floor(x0 + dx * t), (int)floor(y0 + dy * t), fi);
	}
}

static void plot_explicit(const Func *f, int fi, double ppu)
{
	int n = W * OVERSAMPLE, have = 0;
	double px0 = 0, py0 = 0;

	for (int k = 0; k <= n; k++) {
		double px = (double)k / OVERSAMPLE;
		double x = cx + (px - W / 2.0) / ppu;
		double y = eval(f->e, x, 0.0);
		double py;

		if (!isfinite(y)) {
			have = 0;
			continue;
		}
		py = H / 2.0 - (y - cy) * ppu;
		if (!isfinite(py)) {
			have = 0;
			continue;
		}
		if (have)
			seg(px0, py0, px, py, fi);
		px0 = px;
		py0 = py;
		have = 1;
	}
}

static int sign_change(double a, double b)
{
	return (a < 0 && b > 0) || (a > 0 && b < 0);
}

static void plot_implicit(const Func *f, int fi, double ppu)
{
	double *r0 = xcalloc((size_t)W, sizeof(double));
	double *r1 = xcalloc((size_t)W, sizeof(double));

	for (int j = 0; j < H; j++) {
		double y = cy + (H / 2.0 - j - 0.5) / ppu, *t;

		for (int i = 0; i < W; i++)
			r1[i] = eval(f->e, cx + (i + 0.5 - W / 2.0) / ppu, y);

		for (int i = 0; i < W; i++) {
			double v = r1[i];
			if (v == 0.0)
				gdot(i, j, fi);
			if (i + 1 < W && sign_change(v, r1[i + 1])) {
				if (fabs(v) < fabs(r1[i + 1])) gdot(i, j, fi);
				else                           gdot(i + 1, j, fi);
			}
			if (j > 0 && sign_change(r0[i], v)) {
				if (fabs(r0[i]) < fabs(v)) gdot(i, j - 1, fi);
				else                       gdot(i, j, fi);
			}
		}
		t = r0; r0 = r1; r1 = t;
	}
	free(r0);
	free(r1);
}

static void build_frame(void)
{
	size_t n = (size_t)cols * rows;
	double ppu = ppu_now();

	memset(gm, 0, n); memset(am, 0, n); memset(gc, 0, n);
	memset(lb, 0, n); memset(lc, 0, n);

	draw_axes(ppu);
	for (int i = 0; i < nfuncs; i++) {
		if (funcs[i].implicit)
			plot_implicit(&funcs[i], i, ppu);
		else
			plot_explicit(&funcs[i], i, ppu);
	}
	if (opt_legend) {
		for (int i = 0; i < nfuncs && i < rows; i++) {
			char b[256];
			snprintf(b, sizeof b, "%.*s", cols > 2 ? cols - 2 : 0, funcs[i].text);
			put_text(i, 1, b, i + 1);
		}
	}
}

static const char *gcol[MAXFUNCS] = {
	"\033[1;96m", "\033[1;93m", "\033[1;95m", "\033[1;92m",
	"\033[1;91m", "\033[1;94m", "\033[1;97m", "\033[1;33m"
};
static const char *axcol = "\033[90m";

typedef struct { int r, g, b; } RGB;

static RGB pal[16], term_fg, term_bg;
static int pal_ok[16], fg_ok, bg_ok;
static char gcol_buf[MAXFUNCS][32];
static char axcol_buf[48];

static int parse_rgb(const char *s, RGB *out)
{
	int v[3];

	if (strncmp(s, "rgb:", 4))
		return 0;
	s += 4;
	for (int i = 0; i < 3; i++) {
		char *e;
		unsigned long x, maxv;
		int nd;

		if (!isxdigit((unsigned char)*s))
			return 0;
		x = strtoul(s, &e, 16);
		nd = (int)(e - s);
		if (nd < 1 || nd > 4)
			return 0;
		maxv = (1ul << (4 * nd)) - 1;
		v[i] = (int)((x * 255 + maxv / 2) / maxv);
		s = e;
		if (i < 2) {
			if (*s != '/')
				return 0;
			s++;
		}
	}
	out->r = v[0]; out->g = v[1]; out->b = v[2];
	return 1;
}

static void parse_osc(const char *p)
{
	RGB c;

	if (!strncmp(p, "10;", 3)) {
		if (parse_rgb(p + 3, &c)) { term_fg = c; fg_ok = 1; }
	} else if (!strncmp(p, "11;", 3)) {
		if (parse_rgb(p + 3, &c)) { term_bg = c; bg_ok = 1; }
	} else if (!strncmp(p, "4;", 2)) {
		char *e;
		long n = strtol(p + 2, &e, 10);

		if (e != p + 2 && *e == ';' && n >= 0 && n < 16 &&
		    parse_rgb(e + 1, &c)) {
			pal[n] = c;
			pal_ok[n] = 1;
		}
	}
}

static int has_da1(const char *b, size_t n)
{
	for (size_t i = 0; i + 3 < n; i++) {
		if (b[i] != 033 || b[i + 1] != '[' || b[i + 2] != '?')
			continue;
		for (size_t j = i + 3; j < n; j++) {
			unsigned char ch = (unsigned char)b[j];
			if (ch >= 0x40 && ch <= 0x7e)
				return ch == 'c';
		}
	}
	return 0;
}

static void parse_replies(const char *b, size_t n)
{
	for (size_t i = 0; i + 2 < n; i++) {
		char tmp[96];
		size_t j, k = 0;

		if (b[i] != 033 || b[i + 1] != ']')
			continue;
		for (j = i + 2; j < n && b[j] != 007 && b[j] != 033; j++)
			if (k < sizeof tmp - 1)
				tmp[k++] = b[j];
		tmp[k] = 0;
		parse_osc(tmp);
		i = j;
	}
}

static double now_ms(void)
{
	struct timespec ts;

	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

static int query_terminal(void)
{
	static char buf[8192];
	size_t len = 0;
	struct termios o, t;
	double end;
	char q[512];
	int fd, ql = 0;
	const char *term = getenv("TERM");

	if (!term || !*term || !strcmp(term, "dumb"))
		return 0;
	fd = open("/dev/tty", O_RDWR | O_NOCTTY);
	if (fd < 0)
		return 0;
	if (tcgetattr(fd, &o) < 0) {
		close(fd);
		return 0;
	}
	t = o;
	t.c_lflag &= ~(tcflag_t)(ECHO | ICANON | ISIG);
	t.c_cc[VMIN] = 0;
	t.c_cc[VTIME] = 0;
	tcsetattr(fd, TCSANOW, &t);
	tcflush(fd, TCIFLUSH);

	ql += snprintf(q + ql, sizeof q - (size_t)ql, "\033]10;?\033\\\033]11;?\033\\");
	for (int i = 0; i < 16; i++)
		ql += snprintf(q + ql, sizeof q - (size_t)ql, "\033]4;%d;?\033\\", i);
	ql += snprintf(q + ql, sizeof q - (size_t)ql, "\033[c");
	write_all(fd, q, (size_t)ql);

	end = now_ms() + 400.0;
	for (;;) {
		struct pollfd pf = { fd, POLLIN, 0 };
		double left = end - now_ms();
		ssize_t n;

		if (left <= 0)
			break;
		if (poll(&pf, 1, (int)left + 1) <= 0)
			break;
		n = read(fd, buf + len, sizeof buf - len - 1);
		if (n <= 0)
			break;
		len += (size_t)n;
		if (has_da1(buf, len) || len >= sizeof buf - 1)
			break;
	}
	tcflush(fd, TCIFLUSH);
	tcsetattr(fd, TCSANOW, &o);
	close(fd);

	parse_replies(buf, len);
	if (!bg_ok)
		return 0;
	for (int i = 1; i <= 6; i++)
		if (!pal_ok[i] || !pal_ok[i + 8])
			return 0;
	if (!fg_ok && pal_ok[7]) {
		term_fg = pal[7];
		fg_ok = 1;
	}
	return fg_ok;
}

static double lin(int v)
{
	double c = v / 255.0;
	return c <= 0.04045 ? c / 12.92 : pow((c + 0.055) / 1.055, 2.4);
}

static double lum(RGB c)
{
	return 0.2126 * lin(c.r) + 0.7152 * lin(c.g) + 0.0722 * lin(c.b);
}

static double contrast(RGB a, RGB b)
{
	double la = lum(a), lb2 = lum(b), t;

	if (la < lb2) { t = la; la = lb2; lb2 = t; }
	return (la + 0.05) / (lb2 + 0.05);
}

static int pick_color(int hi, int lo)
{
	double ch = contrast(pal[hi], term_bg), cl = contrast(pal[lo], term_bg);

	if (ch >= 3.0)
		return hi;
	if (cl >= 3.0)
		return lo;
	return ch >= cl ? hi : lo;
}

static int sgr_fg(char *dst, size_t sz, int idx)
{
	return snprintf(dst, sz, "\033[%dm", idx < 8 ? 30 + idx : 90 + idx - 8);
}

static int truecolor_ok(void)
{
	const char *ct = getenv("COLORTERM");

	return ct && (strstr(ct, "truecolor") || strstr(ct, "24bit"));
}

static void apply_theme(void)
{
	static const int hue[6][2] = {
		{ 14, 6 }, { 11, 3 }, { 13, 5 }, { 10, 2 }, { 9, 1 }, { 12, 4 }
	};
	int sel[6];

	for (int i = 0; i < 6; i++) {
		sel[i] = pick_color(hue[i][0], hue[i][1]);
		sgr_fg(gcol_buf[i], sizeof gcol_buf[i], sel[i]);
		gcol[i] = gcol_buf[i];
	}

	snprintf(gcol_buf[6], sizeof gcol_buf[6], "\033[39m");
	gcol[6] = gcol_buf[6];

	{
		int y2 = sel[1] == 11 ? 3 : 11;

		if (contrast(pal[y2], term_bg) >= 3.0)
			sgr_fg(gcol_buf[7], sizeof gcol_buf[7], y2);
		else
			snprintf(gcol_buf[7], sizeof gcol_buf[7], "\033[1;39m");
		gcol[7] = gcol_buf[7];
	}

	if (truecolor_ok()) {
		RGB m = term_fg;

		for (int pct = 40; pct <= 80; pct += 10) {
			m.r = term_bg.r + (term_fg.r - term_bg.r) * pct / 100;
			m.g = term_bg.g + (term_fg.g - term_bg.g) * pct / 100;
			m.b = term_bg.b + (term_fg.b - term_bg.b) * pct / 100;
			if (contrast(m, term_bg) >= 2.5)
				break;
		}
		snprintf(axcol_buf, sizeof axcol_buf, "\033[38;2;%d;%d;%dm", m.r, m.g, m.b);
	} else if (contrast(pal[8], term_bg) >= 2.0) {
		snprintf(axcol_buf, sizeof axcol_buf, "\033[90m");
	} else {
		snprintf(axcol_buf, sizeof axcol_buf, "\033[2;39m");
	}
	axcol = axcol_buf;
}

static void app_braille(unsigned m)
{
	char b[3];
	b[0] = (char)0xE2;
	b[1] = (char)(0xA0 | (m >> 6));
	b[2] = (char)(0x80 | (m & 0x3F));
	app_n(b, 3);
}

static char ascii_char(unsigned m, int graph)
{
	if (graph)
		return '*';
	if (m == 0x09 || m == 0x12 || m == 0x24 || m == 0xC0)
		return '-';
	if ((m & 0x47) == 0x47 || (m & 0xB8) == 0xB8)
		return '|';
	return '+';
}

static void render(int interactive)
{
	const char *cur = NULL;

	olen = 0;
	if (interactive)
		app_s("\033[?2026h");

	for (int r = 0; r < rows; r++) {
		size_t rowstart;

		if (interactive) {
			char b[32];
			int n = snprintf(b, sizeof b, "\033[%d;1H", r + 1);
			app_n(b, (size_t)n);
		}
		rowstart = olen;

		for (int c = 0; c < cols; c++) {
			size_t i = (size_t)r * cols + c;
			unsigned m = 0;
			char ch = ' ';
			int graph = 0;
			const char *col = NULL;

			if (lb[i] && lc[i]) {
				ch = lb[i];
				col = gcol[(lc[i] - 1) % MAXFUNCS];
			} else if (gm[i]) {
				m = gm[i];
				graph = 1;
				col = gcol[(gc[i] - 1) % MAXFUNCS];
			} else if (lb[i]) {
				ch = lb[i];
				col = axcol;
			} else if (am[i]) {
				m = am[i];
				col = axcol;
			}

			if (usecolor && col && col != cur) {
				app_s(col);
				cur = col;
			}
			if (m)
				opt_ascii ? app_c(ascii_char(m, graph)) : app_braille(m);
			else
				app_c(ch);
		}

		if (!interactive) {
			while (olen > rowstart && obuf[olen - 1] == ' ')
				olen--;
			if (cur) {
				app_s("\033[0m");
				cur = NULL;
			}
			app_c('\n');
		}
	}
	if (interactive)
		app_s("\033[0m\033[?2026l");
}

static struct termios orig_tio;
static int raw_on, screen_on;
static volatile sig_atomic_t got_winch, got_stop;

static void term_restore(void)
{
	static const char s[] = "\033[0m\033[?7h\033[?25h\033[?1049l";
	if (screen_on) {
		write_all(1, s, sizeof s - 1);
		screen_on = 0;
	}
	if (raw_on) {
		tcsetattr(0, TCSAFLUSH, &orig_tio);
		raw_on = 0;
	}
}

static void term_setup(void)
{
	static const char s[] = "\033[?1049h\033[?25l\033[?7l\033[2J";
	struct termios t;

	if (tcgetattr(0, &orig_tio) < 0)
		die("cannot read terminal attributes");
	t = orig_tio;
	t.c_iflag &= ~(tcflag_t)(BRKINT | ICRNL | INPCK | ISTRIP | IXON);
	t.c_lflag &= ~(tcflag_t)(ECHO | ICANON | IEXTEN | ISIG);
	t.c_cc[VMIN] = 1;
	t.c_cc[VTIME] = 0;
	tcsetattr(0, TCSAFLUSH, &t);
	raw_on = 1;
	atexit(term_restore);
	write_all(1, s, sizeof s - 1);
	screen_on = 1;
}

static void on_sig(int s)
{
	if (s == SIGWINCH)
		got_winch = 1;
	else
		got_stop = 1;
}

static void get_size(int *c, int *r)
{
	struct winsize ws;
	const char *ec, *el;

	if (ioctl(1, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 0 && ws.ws_row > 0) {
		*c = ws.ws_col;
		*r = ws.ws_row;
		return;
	}
	ec = getenv("COLUMNS");
	el = getenv("LINES");
	if (ec && el && atoi(ec) > 0 && atoi(el) > 0) {
		*c = atoi(ec);
		*r = atoi(el);
		return;
	}
	if (ioctl(0, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 0 && ws.ws_row > 0) {
		*c = ws.ws_col;
		*r = ws.ws_row;
		return;
	}
	*c = 80;
	*r = 24;
}

static unsigned char ibuf[256];
static int ilen, ipos;

static int next_byte(int ms)
{
	struct pollfd pf = { 0, POLLIN, 0 };
	int r;
	ssize_t n;

	if (ipos < ilen)
		return ibuf[ipos++];
	r = poll(&pf, 1, ms);
	if (r < 0)
		return errno == EINTR ? -1 : -2;
	if (r == 0)
		return -1;
	if (pf.revents & (POLLERR | POLLNVAL))
		return -2;
	n = read(0, ibuf, sizeof ibuf);
	if (n < 0)
		return (errno == EINTR || errno == EAGAIN) ? -1 : -2;
	if (n == 0)
		return -2;
	ilen = (int)n;
	ipos = 0;
	return ibuf[ipos++];
}

static int read_escape(void)
{
	int c = next_byte(40), f;

	if (c < 0)
		return -1;
	if (c == 'O') {
		f = next_byte(40);
		return f < 0 ? 0 : f;
	}
	if (c == '[') {
		for (;;) {
			f = next_byte(40);
			if (f < 0)
				return 0;
			if (f >= 0x40 && f <= 0x7e)
				return f;
		}
	}
	return 0;
}

static void pan(double dx, double dy)
{
	double st = 12.0 / ppu_now();
	cx += dx * st;
	cy += dy * st;
}

static void zoom_by(double f)
{
	zoom *= f;
	if (zoom < 1e-4) zoom = 1e-4;
	if (zoom > 1e5)  zoom = 1e5;
}

static int run_interactive(void)
{
	struct sigaction sa;
	int c0, r0, dirty = 1, quit = 0;

	memset(&sa, 0, sizeof sa);
	sa.sa_handler = on_sig;
	sigemptyset(&sa.sa_mask);
	sigaction(SIGINT, &sa, NULL);
	sigaction(SIGTERM, &sa, NULL);
	sigaction(SIGHUP, &sa, NULL);
	sigaction(SIGQUIT, &sa, NULL);
	sigaction(SIGWINCH, &sa, NULL);

	term_setup();
	get_size(&c0, &r0);
	canvas_resize(c0, r0);

	while (!quit && !got_stop) {
		int c;

		if (got_winch) {
			got_winch = 0;
			get_size(&c0, &r0);
			canvas_resize(c0, r0);
			dirty = 1;
		}
		if (dirty) {
			build_frame();
			render(1);
			write_all(1, obuf, olen);
			dirty = 0;
		}

		c = next_byte(100);
		while (c != -1) {
			if (c == -2) {
				quit = 1;
				break;
			}
			switch (c) {
			case 3: case 'q': case 'Q':
				quit = 1;
				break;
			case 27: {
				int f = read_escape();
				if (f == -1) quit = 1;
				else if (f == 'A') { pan(0, 1);  dirty = 1; }
				else if (f == 'B') { pan(0, -1); dirty = 1; }
				else if (f == 'C') { pan(1, 0);  dirty = 1; }
				else if (f == 'D') { pan(-1, 0); dirty = 1; }
				break;
			}
			case 'k': pan(0, 1);  dirty = 1; break;
			case 'j': pan(0, -1); dirty = 1; break;
			case 'l': pan(1, 0);  dirty = 1; break;
			case 'h': pan(-1, 0); dirty = 1; break;
			case '+': case '=': zoom_by(1.25);       dirty = 1; break;
			case '-': case '_': zoom_by(1 / 1.25);   dirty = 1; break;
			case '0': case 'r':
				zoom = zoom0; cx = cy = 0;
				dirty = 1;
				break;
			default:
				break;
			}
			if (quit)
				break;
			c = next_byte(0);
		}
	}
	term_restore();
	return 0;
}

static void usage(FILE *f)
{
	fprintf(f,
"peamath " VERSION " - plot math functions in the terminal\n"
"\n"
"Usage:\n"
"  peamath -g EXPR [-g EXPR ...] [options]\n"
"  peamath EXPR\n"
"\n"
"Options:\n"
"  -g, --graph EXPR   function or equation to plot (repeatable, up to %d)\n"
"  -z, --zoom N       initial zoom factor (default 1)\n"
"  -l, --legend       print the expressions in the top-left corner\n"
"  -a, --ascii        plain ASCII instead of braille characters\n"
"  -v, --version      print version and exit\n"
"  -h, --help         print this help and exit\n"
"\n"
"Examples:\n"
"  peamath -g \"y=sin(x)\"\n"
"  peamath -g \"x**2+y**2=1\"\n"
"  peamath -g \"y=ln(x)\" -g \"y=exp(x)\" -l\n"
"  peamath -g \"y=tan(x)\" -z 2\n"
"  peamath -g \"x**3-3*x=y**2\"\n"
"\n"
"Expressions: + - * / %% ^ **  ( )  x y pi e tau  sin cos tan asin acos atan\n"
"  sinh cosh tanh ln log log2 log10 exp sqrt cbrt abs floor ceil round sign\n"
"  atan2(y,x) pow(a,b) min(a,b) max(a,b) mod(a,b) hypot(a,b) logb(base,x)\n"
"  Implicit multiplication works: 2x, 3sin(x), x(x+1).\n"
"\n"
"Keys: arrows or hjkl move, + / - zoom, 0 reset, q / Esc / Ctrl-C quit.\n",
	        MAXFUNCS);
}

static void bad_usage(const char *msg, const char *arg)
{
	fprintf(stderr, "peamath: %s%s%s\nTry 'peamath -h' for help.\n",
	        msg, arg ? " " : "", arg ? arg : "");
	exit(2);
}

int main(int argc, char **argv)
{
	const char *exprs[MAXFUNCS];
	int ne = 0;

	if (argc == 1) {
		usage(stdout);
		return 0;
	}

	for (int i = 1; i < argc; i++) {
		const char *a = argv[i], *val = NULL;
		int is_graph = 0, is_zoom = 0;

		if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
			usage(stdout);
			return 0;
		} else if (!strcmp(a, "-v") || !strcmp(a, "--version")) {
			puts("peamath " VERSION);
			return 0;
		} else if (!strcmp(a, "-l") || !strcmp(a, "--legend")) {
			opt_legend = 1;
		} else if (!strcmp(a, "-a") || !strcmp(a, "--ascii")) {
			opt_ascii = 1;
		} else if (!strcmp(a, "-g") || !strcmp(a, "--graph")) {
			if (++i >= argc)
				bad_usage("option requires an argument:", a);
			val = argv[i];
			is_graph = 1;
		} else if (!strncmp(a, "--graph=", 8)) {
			val = a + 8;
			is_graph = 1;
		} else if (!strcmp(a, "-z") || !strcmp(a, "--zoom")) {
			if (++i >= argc)
				bad_usage("option requires an argument:", a);
			val = argv[i];
			is_zoom = 1;
		} else if (!strncmp(a, "--zoom=", 7)) {
			val = a + 7;
			is_zoom = 1;
		} else if (a[0] == '-' && a[1]) {
			bad_usage("unknown option", a);
		} else {
			val = a;
			is_graph = 1;
		}

		if (is_graph) {
			if (ne >= MAXFUNCS)
				bad_usage("too many graphs, the limit is", "8");
			exprs[ne++] = val;
		} else if (is_zoom) {
			char *end;
			double z = strtod(val, &end);
			if (*end || !(z > 0) || !isfinite(z))
				bad_usage("invalid zoom factor", val);
			zoom0 = zoom = z;
		}
	}

	if (ne == 0)
		bad_usage("nothing to plot, use -g EXPR", NULL);

	for (int k = 0; k < ne; k++) {
		Node *l = NULL, *r = NULL;
		char err[128];
		int ep = 0;
		Func *f;

		if (parse_equation(exprs[k], &l, &r, err, sizeof err, &ep) < 0) {
			fprintf(stderr, "peamath: %s\n  %s\n  %*s^\n", err, exprs[k], ep, "");
			return 1;
		}
		f = &funcs[nfuncs++];
		f->text = (char *)exprs[k];
		if (!r) {
			f->e = l;
			f->implicit = has_var(l, N_Y);
		} else if (l->type == N_Y && !has_var(r, N_Y)) {
			f->e = r;
		} else if (r->type == N_Y && !has_var(l, N_Y)) {
			f->e = l;
		} else {
			f->implicit = 1;
			f->e = mk(N_SUB, l, r);
		}
	}

	{
		const char *nc = getenv("NO_COLOR");
		usecolor = isatty(1) && !(nc && *nc);
	}

	if (usecolor && query_terminal())
		apply_theme();

	signal(SIGPIPE, SIG_IGN);

	if (isatty(0) && isatty(1))
		return run_interactive();

	{
		int c, r;
		get_size(&c, &r);
		canvas_resize(c, r);
		build_frame();
		render(0);
		write_all(1, obuf, olen);
	}
	return 0;
}
