# Cross-indexing subfind and snapshot data

In the context of the GADGET+SubFind ecosystem, this simple code performs a reverse search so that starting from a list of IDs you could reconstruct to what fof and what fof’s sub-halo they belong to.

_In the following, `cidx` stands for the executable of the code; that is its default name which of course the user may want to change in the Makefile_.

### How it works in brief

1)  The first step is to build, for any given snapshot of interest, catalog files that contains the relevant informations about all the particles that belong to a `fof` halo, i.e. `{masked_id, generation, fof-id, subhalo-id}`, where

    -  `masked_id` is the particle’s `id` maskedfor the stellar-generation part (see the final note if you don’t know what we’re talking about);
    -  `generation` is the stellar generation of a given particles;
       obviously, only star particles and gas particles that underwent star formation have `masked_id != id` and a meaningful `generation`.
    -  `fof-id` is the id of the fof the particle belongs to;
    -  `subhalo-id` is the id of the subhalo the particle belongs to, if any.

    This catalogs are built separately for each type of particles from the subfind outputs and then contains _only_ particles that actually are thought to be part of a fof halo. Particles that are in a fof halo but _not_ in a subhalo are still present in the catalog with a pair `{fof-id, -1}`.

2)  Once the catalogs for the snapshot of interest are built, the user can search any ID from the snapshot in the same catalogs to retrieve the `{fof_id, subhalo_id}` pair relative to that ID.
    _Note that there is no need to mask the ID, the search is performed with the ID found in the snapshot_.
    The ID of interest can be passed to `cidx` either as a file containing the IDs or at command line as arguments.
    `cidx` returns a corresponding list of pairs `{fof_id, subhalo_id}` either as a file or as output to `stdout`; the `fof_id` and `subhalo_id` pair can have the following values:

    -  `{fof_id >=0, subhalo_id >=0 }` the particle belongs to a fof with index `fof_id` _and_ to a subhalo of index `subhalo_id` (remind: the subhalo ids are relative to the fof halo they belong).
    -  `{fof_id >=0, -1 }` the particle belongs to a fof with index `fof_id` but it _does not_ belong to a subhalo.
    -  `{-1, -1}` the particle does not belong to any fof (i.e. its ID has not been found in the catalog files).



### How to invoke `cidx`

The `cidx` command line has the form

```bash
path/to/cidx [ action ] [ input file options ]
```

where the `action` possible values are `-c` to create the catalogs, `-s` to search for a list of particles, `-h` to visualize a brief remind of all the options, and `-ask` to obtain some info on the current build of `cidx`.

| action | meaning      |
| ---- | ---- |
| `-c` | create catalogs      |
| `-s` | search for a list of particles     |
| `-h` | visualize help     |
|  `-ask` | visualize some info on the current build |

In order to both crate catalogs and search a list of particles, the user need to specify _(i)_ the number of the snapshot and subfind file(s) it wants to proceed for, _(ii)_ where the snapshot and subfind file(s) is located, _(iii)_ the basename of the snapshot and subfind file(s).
Note that the default values which correspond to common choices in the GADGET+SubFind ecosystem.

| what                              | option            | default value      | notes                                                        |
| --------------------------------- | ----------------- | ------------------ | ------------------------------------------------------------ |
| snapshot base name                | -snapf [name]     | snap_ + snapnum    |                                                              |
| subfind base name                 | -subf [name]      | sub_ + snapnum     |                                                              |
| snapshot number                   | -num [number]     | _none_             | the leading zeros must<br />be included (es. 002, 017, 090, …) |
| snapshot folder                   | -snapdir [folder] | snapdir_ + snapnum |                                                              |
| subfind folder                    | -subfdir [folder] | groups_ + snapnum  |                                                              |
| working directory                 | -wdir [folder]    | ./                 |                                                              |
| number of<br />stellar generation | -g [num]          | 4                  |                                                              |



### How to create the catalogs

Creating catalogs from a snapshot may be as simple as in the following examples

1.  creating catalogs for the snapshot `nnn` whose file(s) are named `snap_nnn` in the standard `snapdir_nnn/` location

    ```bash
    path/to/cidx -c -num nnn
    ```

    The snapshot could be split in multiple files `snap_nnn.0`, `snap_nnn.1`, …
    The subfind failes will be searched as `groups_nnn/sub_nnn.x`.

2.  creating catalogs for the snapshot `nnn` whose file(s) are named `snap_nnn` in the current working dir

    ```bash
    path/to/cidx -c -num nnn -snapdir ./ 
    ```

    The snapshot could be split in multiple files `snap_nnn.0`, `snap_nnn.1`, …

    The subfind failes will be searched as `groups_nnn/sub_nnn.x`.

3.  creating catalogs for the snapshot `nnn` whose files are named `mysnap_from_themoon.le_nnn` under the non-standard location `snapdata_nnn`, and for the subfind file(s) stored under `haloes_nnn/galaxies_nnn.x`:

    ```bash
    path/to/cidx -c -num nnn -snapdir snapdata_ -snapf mysnap_from_the_moon.le_ -subfdir haloes_ -subf galaxies_
    ```

    As before, the snapshot could be split in multiple files `snapdata_nnn.0`, `snapdata_nnn.1`, …

    After the run, in the current working directory there will be the files `nnn_catalog_type_x` where `x` ranges in `[0:5]`. Each file contains the data corresponding to the different particle types.
    The additional file `timings` contains information about the `cidx` run time.

### How to search for a list of particles

Searching for a list of particles can be done by:

1.  creating a file with the IDs of the particles;
    in this case, the user can inform `cidx` about his/her choice by the `-list` argument:

    ```bash
    psth/to/cidx -s -num nnn -list my_list_file [type]
    ```

    where `my_list_file` is the name of the list file (the path to the file must be included if different than ./) and `type` is a specifier that, if needed, informs `cidx` that all the ids in the file are of the given type.
    Not specifying `type` or specifying `-1` means that the IDs refer to particles of mixed types.
    Entering `type` `-2` amounts to specify that the file contains type information for each particle.
    _(see below “how to create a list of IDs”)_

2.  passing the (short) list of IDs directly through command line.
    TBD.



##### how to create a list of IDs

There are 2 types of list file: a file that contains only the IDs of the particles and a file that contains the pair `{id, type}` where `type` is the `[0:5]` type of each particle.
The first type of file, that which contains only the IDs, has the very simple binary format:
`id_size:int4`,`num_of_particles:int8`,`[id0, id1, id2, ...]:id_type`

where `id_size` is a 4bytes integer that specifies the size in bytes of the IDs stored in the file; the size may be either 4 or 8. `num_of_particles` is the number of particles that are in the file and it must be an int8.
The subsequent IDs are stored linearly, and `id_type` is either `uint4_t` or `uint8_t`, depending on the value of `id_size`.

<span style="color: red"> TBD: how to specify multiple-file list </span>

<span style="color: red"> TBD: how to create the second type of files </span>

Note that `cidx` is built so that you can compile with only a precise value for `id_type`, so either `uint4_t` or `uint8_t` (see below, “How to compile”). In case that the IDs are of a different type, the code warns you about recompiling and stops.

##### How to compile

`make` :D

_A more thorough explanation of the compilation option will follow asap_



