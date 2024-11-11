# Linalg Library

Original library: [https://gitlab.com/enrico_facca/linear_algebra](https://gitlab.com/enrico_facca/linear_algebra)

* Clone the repository
```
git clone git@github.com:ncrescenzio/linear_algebra.git
```

* Build the library
```bash
cd globals && mkdir build
cmake -S . -B build
cmake --build build
```

* Install the library
```
cmake --install build --prefix /path/to/install/dir
```

## Dependencies

* [`globals` library](https://github.com/ncrescenzio/globals)
* `lapack`
* `blas`
* Optional: `agmg`, `hsl`, `arpack`
