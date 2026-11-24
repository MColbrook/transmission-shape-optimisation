// Run from the repository root: node certificates/check_transition_endpoints.mjs
// All signs and bounds below are evaluated with rational arithmetic.
import assert from 'node:assert/strict';

function gcd(a, b) {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b !== 0n) [a, b] = [b, a % b];
  return a;
}

function fraction(a, b = 1n) {
  assert(b > 0n);
  const g = gcd(a, b);
  return [a / g, b / g];
}

function decimal(value) {
  const [whole, digits = ''] = value.split('.');
  return fraction(BigInt(whole + digits), 10n ** BigInt(digits.length));
}

const add = ([a, b], [c, d]) => fraction(a * d + c * b, b * d);
const subtract = ([a, b], [c, d]) => fraction(a * d - c * b, b * d);
const multiply = ([a, b], [c, d]) => fraction(a * c, b * d);
const less = ([a, b], [c, d]) => a * d < c * b;
const square = (a) => multiply(a, a);
const zero = fraction(0n);
const one = fraction(1n);

function product(a, b) {
  const endpoints = [multiply(a[0], b[0]), multiply(a[0], b[1]),
    multiply(a[1], b[0]), multiply(a[1], b[1])];
  endpoints.sort((x, y) => less(x, y) ? -1 : less(y, x) ? 1 : 0);
  return [endpoints[0], endpoints[3]];
}

function series(z, order, modified) {
  assert(less(zero, z) && less(z, fraction(31n)));
  let term = one;
  let sum = one;
  const sign = modified ? 1n : -1n;
  for (let k = 1n; k <= 30n; k++) {
    term = multiply(term, multiply(z, fraction(sign, 4n * k * (k + order))));
    sum = add(sum, term);
  }
  const next = multiply(term, multiply(z, fraction(sign, 124n * (31n + order))));
  // For modified series the tail after term 31 has term ratios below 1/100.
  // The alternating series lie between the partial sums of orders 30 and 31.
  const interval = modified
    ? [sum, add(sum, multiply(next, fraction(100n, 99n)))]
    : [add(sum, next), sum];
  if (order === 0n) return interval;
  const factor = multiply(z, fraction(1n, 2n));
  return product(interval, [factor, factor]);
}

const j = (z) => series(z, 0n, false);
const h = (z) => series(z, 1n, false);
const i = (z) => series(z, 0n, true);
const k = (z) => series(z, 1n, true);

// J1 is positive on (0,2.41] by its first two alternating terms.
// Use J0 zero spacing above 2.4 and the J0/J1 interlacing to index the roots.
const alphaLo = decimal('5.78318');
const alphaHi = decimal('5.78319');
const betaLo = decimal('14.68196');
const betaHi = decimal('14.68198');
const zetaLo = decimal('30.47126');
const zetaHi = decimal('30.47127');
assert(less(zero, j(alphaLo)[0]) && less(j(alphaHi)[1], zero));
assert(less(zero, h(betaLo)[0]) && less(h(betaHi)[1], zero));
assert(less(j(zetaLo)[1], zero) && less(zero, j(zetaHi)[0]));
assert(less(square(decimal('2.4')), alphaLo));
assert(less(alphaHi, square(decimal('2.41'))));
assert(less(square(decimal('5.5')), zetaLo));
assert(less(zetaHi, square(decimal('5.53'))));
assert(less(alphaHi, betaLo) && less(betaHi, zetaLo));
assert(less(fraction(57n, 10n), alphaLo));
assert(less(fraction(169n), multiply(fraction(2n), multiply(alphaLo, betaLo))));

function diskDeterminant(n, x) {
  const y = multiply(n, x);
  const first = product(j(x), h(y));
  const second = product(h(x), j(y));
  return [subtract(first[0], second[1]), subtract(first[1], second[0])];
}

function radialRootBracket(contrast, lower, upper) {
  const n = decimal(contrast);
  const lo = decimal(lower);
  const hi = decimal(upper);
  assert(less(lo, hi));
  // The two intervals lie in the active disk branches.
  assert(less(hi, alphaLo));
  assert(less(betaHi, multiply(n, lo)));
  assert(less(multiply(n, hi), zetaLo));
  assert(less(zero, diskDeterminant(n, lo)[0]));
  assert(less(diskDeterminant(n, hi)[1], zero));
  return [n, lo, hi];
}

const [n, lo, hi] = radialRootBracket('5.281738', '5.1803707', '5.1803709');
assert(less(multiply(n, square(hi)), fraction(142n)));
const formUpper = add(subtract(fraction(6923409n, 1000000n),
  multiply(multiply(add(n, one), lo), fraction(2941129141n, 10000000000n))),
  multiply(multiply(n, square(hi)), fraction(7471409n, 400000000n)));
assert(less(formUpper, fraction(-328n, 100000000n)));
console.log('At contrast 5.281738 the active radial root lies in (5.1803707, 5.1803709).');
console.log(`Exact upper bound for the trial form divided by pi: ${formUpper.join('/')}.`);

function plateDeterminant(z) {
  const first = product(j(z), k(z));
  const second = product(h(z), i(z));
  return [add(first[0], second[0]), add(first[1], second[1])];
}

const gammaLo = decimal('10.21582622');
const gammaHi = decimal('10.21582624');
assert(less(alphaHi, gammaLo) && less(gammaHi, betaLo));
assert(less(decimal('10.2'), gammaLo));
assert(less(fraction(104n), square(gammaLo)) && less(square(gammaHi), fraction(105n)));
assert(less(zero, plateDeterminant(gammaLo)[0]));
assert(less(plateDeterminant(gammaHi)[1], zero));
const [nLo, xLo] = radialRootBracket('5.4254619', '4.38586394', '4.38586397');
const [nHi, , xHi] = radialRootBracket('5.4254621', '4.38586333', '4.38586337');
assert(less(square(gammaHi), multiply(nLo, square(xLo))));
assert(less(multiply(nHi, square(xHi)), square(gammaLo)));
console.log('The plate cutoff lies in (5.4254619, 5.4254621).');

// 2 (5.7)(10.2) (6/sqrt(5)) - 105 > 207.
assert(less(multiply(fraction(5n), square(fraction(312n))), square(decimal('697.68'))));
assert(104n * 207n > 142n ** 2n);
console.log('Checked t_D(5.281738) < 142 and the integer product margin.');
