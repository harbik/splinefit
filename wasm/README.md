# splinefit

B-spline curve fitting in WebAssembly. Fit smoothing, interpolating, or cardinal cubic splines to data and evaluate, integrate, or find roots -- all in the browser, with no server required.

Built from a pure-Rust translation of the classic [Dierckx FITPACK](http://www.netlib.org/dierckx/) library (the same engine behind SciPy's `splrep`/`splev`).

## Installation

```sh
npm install splinefit
```

Or use directly from a CDN:

```html
<script type="module">
  import init, { CubicSpline } from "https://esm.sh/splinefit";
  await init();
  // ...
</script>
```

## Quick start

```js
import init, { CubicSpline } from "splinefit";
await init();

// Sample data: sin(x) on [0, 2pi]
const x = Float64Array.from({ length: 50 }, (_, i) => i * 2 * Math.PI / 49);
const y = x.map(Math.sin);

// Fit a smoothing spline (rms = 0.05)
const spline = CubicSpline.smoothing(x, y, 0.05);

// Evaluate at 200 points
const xNew = Float64Array.from({ length: 200 }, (_, i) => i * 2 * Math.PI / 199);
const yFit = spline.evaluate(xNew); // Float64Array

console.log(yFit[0]);  // ~0.0
console.log(yFit[100]); // ~1.0
```

## API

### Constructors

All constructors take `Float64Array` inputs. `x` must be strictly increasing and have the same length as `y` (minimum 4 points).

#### `CubicSpline.smoothing(x, y, rms)`

Fit a smoothing spline. `rms` controls the trade-off between smoothness and closeness of fit: smaller values produce more knots and a tighter fit.

```js
const spline = CubicSpline.smoothing(x, y, 0.05);
```

#### `CubicSpline.interpolating(x, y)`

Fit an interpolating spline that passes exactly through every data point.

```js
const spline = CubicSpline.interpolating(x, y);
```

#### `CubicSpline.cardinal(x, y, dt)`

Fit a spline on a fixed equidistant knot grid with spacing `dt`.

```js
const spline = CubicSpline.cardinal(x, y, 0.5);
```

### Methods

#### `spline.evaluate(x) -> Float64Array`

Evaluate the spline at each value in `x`.

```js
const yFit = spline.evaluate(new Float64Array([0.5, 1.0, 1.5]));
```

#### `spline.integral(a, b) -> number`

Compute the definite integral of the spline over `[a, b]`.

```js
const area = spline.integral(0, Math.PI); // integral of sin(x) from 0 to pi ~ 2.0
```

#### `spline.roots() -> Float64Array`

Find all interior zeros of the spline, returned in ascending order. Zeros at the domain boundaries may not be included.

```js
const zeros = spline.roots(); // e.g. [3.14159...]
```

#### `spline.knots() -> Float64Array`

Return the knot vector.

```js
console.log(spline.knots()); // Float64Array [0, 0, 0, 0, ..., 6.28, 6.28, 6.28, 6.28]
```

#### `spline.coefficients() -> Float64Array`

Return the B-spline coefficients.

#### `spline.num_knots -> number`

Number of knots (getter property, no parentheses).

```js
console.log(spline.num_knots); // e.g. 12
```

## TypeScript

Type definitions (`.d.ts`) are included automatically. Your editor will provide full autocompletion and type checking.

## Examples

### Smoothing noisy sensor data

```js
const spline = CubicSpline.smoothing(timestamps, readings, 0.1);
const smooth = spline.evaluate(timestamps);
```

### Finding zero crossings

```js
const spline = CubicSpline.interpolating(x, y);
const crossings = spline.roots();
console.log(`Signal crosses zero at: ${Array.from(crossings).join(", ")}`);
```

### Computing area under a curve

```js
const spline = CubicSpline.interpolating(x, y);
const area = spline.integral(x[0], x[x.length - 1]);
console.log(`Total area: ${area}`);
```

## How it works

This package is compiled from the Rust crate [`splinefit`](https://crates.io/crates/splinefit) to WebAssembly using `wasm-pack`. It runs the same numerical algorithms as SciPy's `splrep`/`splev` (Paul Dierckx' FITPACK), translated line-by-line from the original Fortran into Rust.

No server, no network calls, no native dependencies -- the spline fitting runs entirely in your browser or Node.js process.

## License

Apache-2.0 OR MIT
