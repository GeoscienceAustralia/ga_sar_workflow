## Product packaging

Once the user has processed a stack, it's often desirable they can distribute that data to others in a standard format.
`gamma_insar` provides a packaging script which will export a stack's data as an (ODC)[https://github.com/opendatacube/datacube-core] dataset via (eodatasets)[https://github.com/opendatacube/eo-datasets] which allows users to distribute the stack as an Open Data Cube for users to easily analyse and produce higher level products upon.

The process is incredibly simple, as this example shows, the packaging of a backscatter for a specific Sentinel-1 track/frame:

```BASH
package \
    --track T133D --frame F35S \
    --input-dir /path/to/workflow_output_data \
    --pkgdir /path/to/track_frame_pkg_output
```

Additional options exist for handling various production scenarios (eg: how to handle existing data if you're re-processing an existing stack), how to package subsets of the data (eg: single polarisation), etc - please refer to `package --help` for further details.
