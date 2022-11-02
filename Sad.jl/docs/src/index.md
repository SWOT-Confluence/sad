```@meta
CurrentModule = Sad
```

# Sad

Documentation for [Sad](https://gitlab.com/kandread/Sad.jl).

## Running on Confluence

There are two ways to run Sad with Confluence files: directly through the `scripts/swot.jl` script, or with the `kandread/sad` Docker image.

In both cases, you have to provide two arguments: the NetCDF file containing the SWOT observations (e.g., `swotobs.nc`) and the NetCDF file with the SWORD data (e.g., `sword.nc`)

```bash
julia scripts/swot.jl swotobs.nc sword.nc
```

or 

```bash
docker run -v /path/to/netcdf/files:/data kandread/sad /data/swotobs.nc /data/sword.nc
```

Currently the output of ``A_0`` and ``n`` is written to `stdout`.


```@index
```

```@autodocs
Modules = [Sad]
```
