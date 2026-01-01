# Lera

**version: 0.1.0**

An optimized linear algebra library done in pure Rust.

- - -

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
lera = "0.1.0"
```

Or clone the repository and add it as a dependency:

```toml
[dependencies]
lera = { path = "path/to/lera" }
```

- - - 

## Features

- Fully specialized matrices/vector operations for common operations/sizes, fallbacks to generic implementation if not
  in said cases.
- Support for complex numbers.
- fully generic between real and complex f32 and f64.
- No dependencies.
- MIT licensed.

- - - 

## Specialization

The following operations are fully specialized:

- Matrix operations up to 4x4
    - Fully inlined
    - Fully unrolled in most cases
    - No branches in most cases
- Vector operations up to the fourth dimension
    - fully inlined
    - fully unrolled in most cases
    - No branches in most cases
- Generic fallbacks for everything else (still optimized)

- - - 

## Contributing

Please open an issue or submit a pull request.

- - - 

## License

MIT
